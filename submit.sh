#!/bin/bash
#Submit to the cluster, give it a unique name
#$ -S /bin/bash

#$ -cwd
#$ -V
#$ -l h_vmem=1.9G,h_rt=20:00:00,tmem=1.9G
#$ -pe smp 2

# join stdout and stderr output
#$ -j y
#$ -R y

# source activate /SAN/vyplab/vyplab_reference_genomes/conda_envs/splicing_env/
# --cluster-sync "qsub -l tmem={cluster.tmem},h_vmem={cluster.h_vmem},h_rt={cluster.h_rt} -o submissions/{rule}.{wildcards}.{jobid}.sh.o {cluster.submission_string}" \

snakemake -s downgrade_test.smk \
--conda-prefix "/SAN/vyplab/alb_projects/" \
--use-conda \
--use-singularity \
--singularity-args "--bind /SAN/vyplab/:/SAN/vyplab/,/scratch0/:/scratch0/" \
--jobscript cluster_qsub.sh \
--cluster-config cluster.yaml \
--cluster-sync "qsub -l tmem={cluster.tmem},h_vmem={cluster.h_vmem},h_rt={cluster.h_rt} -o submissions/{rule}.{wildcards}.\$JOB_ID.sh.o {cluster.submission_string}" \
-j 8 \
--rerun-triggers mtime \
--nolock \
--rerun-incomplete --conda-frontend conda --latency-wait 60
