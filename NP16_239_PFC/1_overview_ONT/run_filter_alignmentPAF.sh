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


filter_alignment(){

    cd $2/PAF
    picard FilterSamReads I=$2/$1.bam O=$2/$1.filtered.bam READ_LIST_FILE=$2/PAF/$1_filteredreads.txt FILTER=includeReadList &> $2/PAF/$1.picard.log
    samtools bam2fq $2/$1.filtered.bam| seqtk seq -A > $2/$1.filtered.fa
    samtools sort -O bam -o "$2/$1.filtered.sorted.bam" "$2/$1.filtered.bam"
    

}

parts=(00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20)
for i in ${parts[@]}; do 
  sample=${ALL_SAMPLES_NAMES[0]}${i}
  echo $sample
  filter_alignment ${sample}_mapped ${WKD_ROOT}/5_cupcake/5_align
done
