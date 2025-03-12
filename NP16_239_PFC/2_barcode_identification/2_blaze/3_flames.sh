#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=5:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --mem=200G # specify bytes of memory to reserve
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=s.k.leung@exeter.ac.uk # email address

module load Miniconda2
source activate lrp

export PATH=$PATH:/lustre/projects/Research_Project-MRC148213/lsl693/software/FLAMES/python/
GENOME_FASTA=/lustre/projects/Research_Project-MRC148213/lsl693/references/human/hg38.fa
GFF=/lustre/projects/Research_Project-MRC148213/lsl693/references/human/gencode.v40.annotation.gff3
inputDir=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/Gina/analysis/6_blaze

sc_long_pipeline.py --gff3 ${GFF} --infq ${inputDir}/1_blaze/NP16_239_PFCmatched_reads.fastq -b ${inputDir}/2_minimap2/NP16_239_PFCmatched_reads_filtered_sorted.bam -o  ${inputDir}/3_flames/ 	--genomefa ${GENOME_FASTA} --minimap2_dir /lustre/projects/Research_Project-MRC190311/software/minimap2

