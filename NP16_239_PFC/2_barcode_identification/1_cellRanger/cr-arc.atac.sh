#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p bioseq # partition to use
#SBATCH --time=24:00:00 # maximum walltime
#SBATCH -A Research_Project-BioProduction # research project to submit under
#SBATCH --nodes=1 # specify number of nodesq pq
#SBATCH --exclusive
#SBATCH --job-name=cr-arc-atac

. "/data1/miniconda3/etc/profile.d/conda.sh"
conda activate illumina_bcl


CR_ARC=/data1/software/cellranger-arc-2.0.2/cellranger-arc

bcl2fastq --version
${CR_ARC} --version

${CR_ARC} mkfastq \
  --use-bases-mask=Y90N*,I8N2,Y24,Y100N* \
  --run=../ \
  --id=V0311-atac \
  --csv=ss.cr-arc.csv \
  --output-dir=V0311-atac-demux \
  --barcode-mismatches=0
