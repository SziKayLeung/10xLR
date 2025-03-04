#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p bioseq # partition to use
#SBATCH --time=24:00:00 # maximum walltime
#SBATCH -A Research_Project-BioProduction # research project to submit under
#SBATCH --nodes=1 # specify number of nodesq pq
#SBATCH --exclusive
#SBATCH --job-name=cr7_rna


. "/data1/miniconda3/etc/profile.d/conda.sh"
conda activate illumina_bcl

bcl2fastq --version

bcl2fastq --use-bases-mask=Y28N*,I10,I10N*,Y150N* \
            --create-fastq-for-index-reads \
            --minimum-trimmed-read-length=8 \
            --mask-short-adapter-reads=8 \
            --ignore-missing-positions \
            --ignore-missing-controls \
            --ignore-missing-filter \
            --ignore-missing-bcls \
            --runfolder-dir  .. \
            --output-dir=V0311-gex-demux \
            --sample-sheet=ss.bcl.gex.csv
