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
#SBATCH --output=1b_identify_noTSO-%A_%a.o
#SBATCH --error=1b_identify_noTSO-%A_%a.e

# 28/08/2025: identify reads with no TSO sequence after running BLAZE and cutadapt (to subset for ONT to run via their pipeline)

## ----------------------------------------------------------

module load Miniconda2
source activate lrp

rootDir=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/Gina/analysis/EX167_Hypo/
basecalledReads=${rootDir}/basecalled.fastq

## ----------------------------------------------------------

cd ${rootDir}/6_blaze/1_blaze

# subset the sequences that have not been trimmed by cutadapt i.e. do not have TSO sequences
grep "^@" EX167_HYPmatched_reads_untrimmed.fastq | cut -f 1 | awk -F'#' '{print $2}' - | cut -d'_' -f1 > EX167_HYPmatched_reads_untrimmed_reads.txt

# extract the read IDs from the basecalled reads
seqtk subseq ${basecalledReads} EX167_HYPmatched_reads_untrimmed_reads.txt > ${rootDir}/basecalled_noTSO.fastq

# subset a 10th of the reads
seqtk sample -s100 ${rootDir}/basecalled_noTSO.fastq 0.1 > ${rootDir}/basecalled_noTSO_subsetted.fastq
