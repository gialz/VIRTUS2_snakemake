import os
import pandas as pd
import sys

# Configuration
configfile: "config.yaml"

fastq_dir = config.get("fastq_dir")
sample_df = pd.read_csv(config.get("sample_metadata"))

sample_df["fastq1"] = [os.path.join(fastq_dir, f) for f in sample_df["f1"]]
sample_df["fastq2"] = [os.path.join(fastq_dir, f) for f in sample_df["f2"]]

SAMPLES = sample_df.set_index("sample")[["fastq1", "fastq2"]].to_dict("index")
sample_names = list(SAMPLES.keys())


# Default parameters
NTHREADS = config.get("threads", 8)
KZ_THRESHOLD = config.get("kz_threshold", 0.1)
FILENAME_OUTPUT = config.get("filename_output", "VIRTUS.output.tsv")
OUTFILEPREFIX_HUMAN = config.get("outFileNamePrefix_human", "human")
FASTP_LENGTH = config.get("fastp_length", 40)

# Input files from config
GENOME_DIR_HUMAN = config.get("genomeDir_human", config.get("dir_name_star_human", "STAR_index_human/"))
GENOME_DIR_VIRUS = config.get("genomeDir_virus", config.get("dir_name_star_virus", "STAR_index_virus"))

# All output files
rule all:
    input:
        expand(["results/qc/{sample}_fastp_1.fq.gz", "results/qc/{sample}_fastp_2.fq.gz"], sample=sample_names),

        expand("results/human/{sample}.Log.final.out", sample=sample_names),

        #expand("results/virus_mapping/{sample}_virus.Aligned.sortedByCoord.out.bam", sample=sample_names),
        #expand("results/virus_mapping/{sample}_virus.Log.out", sample=sample_names),
        #expand("results/virus_mapping/{sample}_virus.Log.progress.out", sample=sample_names),
        #expand("results/virus_mapping/{sample}_virus.SJ.out.tab", sample=sample_names),
        #expand("results/virus_mapping/{sample}_virus.Log.final.out", sample=sample_names),

        #expand(["results/unmapped/{sample}_kz_1.fq.paired.fq", "results/unmapped/{sample}_kz_2.fq.paired.fq"], sample=sample_names),
        #expand(["results/unmapped/{sample}_kz_1.fq", "results/unmapped/{sample}_kz_2.fq"], sample=sample_names),

        #expand("results/unmapped/{sample}_unmapped.bam", sample=sample_names),
        #expand(["results/unmapped/{sample}_unmapped.1.fq", "results/unmapped/{sample}_unmapped.2.fq"], sample=sample_names),
        #expand("results/final/{sample}_" + FILENAME_OUTPUT, sample=sample_names),
        #expand("results/virus_mapping/{sample}.filtered.out.bam", sample=sample_names),
        #expand("results/coverage/{sample}_virus.coverage.txt", sample=sample_names),

#downloading rna sequencing data
rule download_sra:
    output:
        f1 = temp("data/{sample}_1.fastq"),
        f2 = temp("data/{sample}_2.fastq")
    params:
        sra_id = "{sample}",
        tempdir = "/scratch0/tmp_{sample}",
        sra_toolkit_path = config.get("sra_toolkit", "")
    threads: 4
    shell:
        """
        mkdir -p {params.tempdir}
        mkdir -p data/
        cd {params.tempdir}
        {params.sra_toolkit_path}prefetch --max-size 70000000000 --force {params.sra_id}
        rm -f ../../data/{params.sra_id}_*.fastq ../../data/{params.sra_id}_*.fastq.gz
        {params.sra_toolkit_path}fasterq-dump {params.sra_id} --split-files -O ../../data/ -t . > ../../logs/{params.sra_id}_fasterq_dump.log 2>&1
        cd -
        rm -rf {params.tempdir}
        """

# Step 1: Quality control with fastp
rule fastp_pe:
    input:
        sample = ["data/{sample}_1.fastq", "data/{sample}_2.fastq"]
    output:
        trimmed1 = "results/qc/{sample}_fastp_1.fq.gz",
        trimmed2 = "results/qc/{sample}_fastp_2.fq.gz",
        html = "results/qc/{sample}_fastp_report.html",
        json = "results/qc/{sample}_fastp_report.json",
    priority: 50
    log:
        "logs/{sample}_fastp_qc.log"
    params:
        length_param = "--length_required {}".format(FASTP_LENGTH)
    threads: 8
    singularity:
        "docker://quay.io/biocontainers/fastp:0.20.0--hdbcaa40_0"
    shell:
        """
        fastp --in1 {input.sample[0]} --in2 {input.sample[1]} \
        --out1 {output.trimmed1} \
        --out2 {output.trimmed2} \
        --thread {threads} \
        --trim_poly_x \
        --html {output.html} \
        --json {output.json} \
        {params.length_param}
        """

# Step 2: Align reads to human genome
rule star_mapping_pe_human:
    input:
        fq1 = "results/qc/{sample}_fastp_1.fq.gz",
        fq2 = "results/qc/{sample}_fastp_2.fq.gz",
        genomeDir = config.get("dir_name_star_human")
    output:
        aligned = temp("results/human/{sample}.Aligned.sortedByCoord.out.bam"),
        mappingstats = "results/human/{sample}.Log.final.out"
    params:
        nthreads = 4,  
        outSAMunmapped = "Within",
        outSAMtype = "BAM SortedByCoordinate",
        readFilesCommand = "zcat",
        outFileNamePrefix = "results/human/{sample}.",
        tmpDir = "/scratch0/tmp_{sample}_star",
    singularity:
        "docker://quay.io/biocontainers/star:2.7.10a--h9ee0642_0"
    shell:
        """
        STAR --runMode alignReads \
        --genomeDir {input.genomeDir} \
        --runThreadN {params.nthreads} \
        --outSAMunmapped {params.outSAMunmapped} \
        --outSAMtype {params.outSAMtype} \
        --outTmpDir {params.tmpDir} \
        --readFilesCommand {params.readFilesCommand} \
        --outFileNamePrefix {params.outFileNamePrefix} \
        --readFilesIn {input.fq1} {input.fq2}
        """
# #Step 3: Extract unmapped reads with samtools
# rule remove_multi: #should have a high priority
#     input:
#         "results/human/{sample}.Aligned.sortedByCoord.out.bam"
#     output:
#         bam = "results/unmapped/{sample}_unmapped.bam",
#         idx = "results/unmapped/{sample}_unmapped.bam.bai"
#     priority: 50
#     log:
#         "results/unmapped/{sample}_samtools_view.log"
#     params:
#         extra = "-f 4"
#     threads: 4
#     singularity:
#         "docker://yyasumizu/bam_filter_polyx:1.3"
#     shell:
#         """
#         mkdir -p results/unmapped
#         samtools view -@ {threads} {params.extra} {input} | \
#         grep -v "uT:A:3" | \
#         samtools view -@ {threads} -bS - > {output.bam} 2> {log}
#         samtools index {output.bam} 2>> {log}
#         """

# # Step 4: Convert unmapped BAM to FASTQ using bedtools bamtofastq
# rule bedtools_bamtofastq:
#     input:
#         "results/unmapped/{sample}_unmapped.bam"
#     output:
#         "results/unmapped/{sample}_unmapped.1.fq",
#         "results/unmapped/{sample}_unmapped.2.fq"
#     log:
#         "results/unmapped/{sample}.separate.log"
#     threads: 4
#     singularity:
#         "docker://quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0"
#     shell:
#         """
#         bedtools bamtofastq \
#         -i {input} \
#         -fq {output[0]} \
#         -fq2 {output[1]} \
#         2> {log}
#         """

# # Step 5: Filter low complexity sequences with kz_filter for read 1
# rule kz_filter_fq1:
#     input:
#         input_fq = "results/unmapped/{sample}_unmapped.1.fq"
#     output:
#         output = "results/unmapped/{sample}_kz_1.fq"
#     params:
#         threshold = KZ_THRESHOLD
#     log:
#         "results/unmapped/{sample}_kz_1.separate.log"
#     singularity:
#         "docker://yyasumizu/ko:0.1"
#     shell:
#         """
#         kz --filter --threshold {params.threshold} < {input.input_fq} > {output.output} 2> {log}
#         """

# # Step 6: Filter low complexity sequences with kz_filter for read 2
# rule kz_filter_fq2:
#     input:
#         input_fq = "results/unmapped/{sample}_unmapped.2.fq"
#     output:
#         output = "results/unmapped/{sample}_kz_2.fq"
#     params:
#         threshold = KZ_THRESHOLD
#     log:
#         "results/unmapped/{sample}_kz_2.separate.log"
#     singularity:
#         "docker://yyasumizu/ko:0.1"
#     shell:
#         """
#         kz --filter --threshold {params.threshold} < {input.input_fq} > {output.output} 2> {log}
#         """

# #Step 7: Re pair filtered reads with fastq
# rule fastq_pair:
#     input:
#         fq1 = "results/unmapped/{sample}_kz_1.fq",
#         fq2 = "results/unmapped/{sample}_kz_2.fq"
#     output:
#         fq1_paired = "results/unmapped/{sample}_kz_1.fq.paired.fq",
#         fq2_paired = "results/unmapped/{sample}_kz_2.fq.paired.fq"
#     log:
#         "results/unmapped/{sample}_fq.paired.log"
#     singularity:
#         "docker://quay.io/biocontainers/fastq-pair:1.0--he1b5a44_1"
#     shell:
#         """
#         fastq_pair {input.fq1} {input.fq2} 2> {log}
#         """

# # Step 8: Map filtered reads to virus genome with STAR
# rule star_mapping_pe_virus:
#     input:
#         fq1 = "results/unmapped/{sample}_kz_1.fq.paired.fq",
#         fq2 = "results/unmapped/{sample}_kz_2.fq.paired.fq",
#         genomeDir = GENOME_DIR_VIRUS
#     output:
#         aligned = "results/virus_mapping/{sample}_virus.Aligned.sortedByCoord.out.bam",
#         log_out = "results/virus_mapping/{sample}_virus.Log.out",
#         log_progress = "results/virus_mapping/{sample}_virus.Log.progress.out",
#         sj_out = "results/virus_mapping/{sample}_virus.SJ.out.tab",
#         mappingstats = "results/virus_mapping/{sample}_virus.Log.final.out"
#     params:
#         outFileNamePrefix = "results/virus_mapping/{sample}_virus.",
#         outSAMtype = "BAM SortedByCoordinate",
#         outSAMunmapped = "Within",
#         nthreads = 4,
#     log:
#         "results/virus_mapping/{sample}_fq.paired.log"
#     threads: 4
#     singularity:
#         "docker://quay.io/biocontainers/star:2.7.10a--h9ee0642_0"
#     shell:
#         """
#         STAR --runMode alignReads \
#         --genomeDir {input.genomeDir} \
#         --readFilesIn {input.fq1} {input.fq2} \
#         --runThreadN {params.nthreads} \
#         --outSAMunmapped {params.outSAMunmapped} \
#         --outFileNamePrefix {params.outFileNamePrefix} \
#         --outSAMtype {params.outSAMtype}
#         """

# # Step 9: Filter polyX sequences from virus-mapped BAM
# rule bam_filter_polyx:
#     input:
#         input_bam = "results/virus_mapping/{sample}_virus.Aligned.sortedByCoord.out.bam"
#     output:
#         output = "results/virus_mapping/{sample}.filtered.out.bam"
#     singularity:
#         "docker://yyasumizu/bam_filter_polyx:1.3"
#     shell:
#         """
#         samtools view -h {input.input_bam} | \
#         grep -v "AAAAAAAAAAAAAAAAAAAA" | grep -v "TTTTTTTTTTTTTTTTTTTT" | grep -v "TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG" | \
#         samtools view -bS > {output.output} 
#         """

# # Step 10: Calculate coverage with samtools
# rule samtools_coverage:
#     input:
#         input_bam = "results/virus_mapping/{sample}.filtered.out.bam"
#     output:
#         output = "results/coverage/{sample}_virus.coverage.txt"
#     log:
#         "results/coverage/{sample}.coverage.log"
#     singularity:
#         "docker://quay.io/biocontainers/samtools:1.15--h1170115_1"
#     shell:
#         """
#         samtools coverage {input.input_bam} > {output.output} 2> {log}
#         """

# # Step 11: Generate final summary report
# rule mk_summary_virus_count:
#     input:
#         input_STARLog = "results/human/{sample}.Log.final.out",
#         input_virus_cov = "results/coverage/{sample}_virus.coverage.txt"
#     output:
#         f_output = "results/final/{sample}_" + FILENAME_OUTPUT
#     params:
#         layout = "PE",
#     # singularity:
#     #     "docker://yyasumizu/mk_summary_virus_count:2.0"
#     run:
#         df_STARLog = pd.read_csv(input.input_STARLog, sep='\t', header=None, index_col=0)
#         num_reads = int(df_STARLog.loc['                   Uniquely mapped reads number |', 1]) + \
#                     int(df_STARLog.loc['        Number of reads mapped to multiple loci |', 1])
#         df_cov = pd.read_csv(input.input_virus_cov, sep='\t') 
#         df_cov.columns = ['virus', 'startpos', 'endpos', 'numreads', 'covbases', 'coverage', 'meandepth', 'meanbaseq', 'meanmapq']
#         if params.layout == 'PE':
#             df_cov['numreads'] /= 2
#         df_cov['rate_hit'] = df_cov['numreads'] / num_reads
#         df_cov = df_cov.sort_values(by='rate_hit', ascending=False)
#         df_cov = df_cov[df_cov.numreads > 0]
#         df_cov.to_csv(output.f_output, index=None, sep='\t')