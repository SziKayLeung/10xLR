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

# raw short-read FASTQ
SRFASTQ=/lustre/projects/Research_Project-MRC190311/SingleCell/BMIAMP/1_raw/project_11062/11062/V0311-gex-demux/11062

# raw long-read FASTQ basecalled
LRFASTQ=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/Gina/basecalled_sequence

# run BLAZE
# number of cells from short-read sequencing data
cd /lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/Gina/analysis/6_blaze/1_blaze
blaze --force-cells 1378 --output-prefix NP16_239_PFC --threads=12 ${LRFASTQ} &> blaze.log
