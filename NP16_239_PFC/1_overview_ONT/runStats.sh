#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=50:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion


##-------------------------------------------------------------------------

echo Job started on:
date -u

# source config and function
module load Miniconda2/4.3.21
source activate lrp

# load config file provided on command line when submitting job
# Check if a config file was provided on the command line
if [ -z "$1" ]; then
    echo "Error: No config file provided."
    exit 1
fi
config=$(realpath "$1")
echo "Loading config file for project: ${config}" 
source ${config}
source ${SCRIPT_ROOT}/processing/01_source_functions.sh

# create stats output for QC downstream
PrefixOriginal=${ALL_SAMPLES_NAMES[0]}
seqkit stats -a ${WKD_ROOT}/1_basecalled/original/${PrefixOriginal}_merged.fastq > ${WKD_ROOT}/1b_demultiplex_merged/${NAME}_readstats.txt
