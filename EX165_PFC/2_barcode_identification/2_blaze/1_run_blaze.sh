#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=20:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=1_run_blaze.o
#SBATCH --error=1_run_blaze.e

module load Miniconda2
source activate nanopore

# raw long-read FASTQ basecalled
LRFASTQ=$(ls /lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/Gina/P0179_20250326_11623/*)
WKD_ROOT=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/Gina/analysis/EX165_PFC

# merge fastq files
echo ${LRFASTQ}
cat ${LRFASTQ} > ${WKD_ROOT}/basecalled.fastq.gz
gunzip ${WKD_ROOT}/basecalled.fastq.gz


# run BLAZE
# number of cells from short-read sequencing data
mkdir -p ${WKD_ROOT}/6_blaze ${WKD_ROOT}/6_blaze/1_blaze
cd ${WKD_ROOT}/6_blaze/1_blaze
whitelist=/lustre/projects/Research_Project-MRC148213/lsl693/software/cellranger-arc-2.0.2/lib/python/cellranger/barcodes/737K-arc-v1.txt 
blaze --force-cells 2588 --output-prefix EX165_PFC --threads=12 ${WKD_ROOT}/basecalled.fastq --full-bc-whitelist ${whitelist} &> blaze.log