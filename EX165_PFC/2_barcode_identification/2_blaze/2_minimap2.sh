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
#SBATCH --output=2_minimap2.o
#SBATCH --error=2_minimap2.e

module load Miniconda2
source activate lrp

run_minimap2(){

  name=$(basename $1 .fastq.gz)

  if [ -f $2/${name}_sorted.sam ]; then
    echo -e "Minimap2: \e[32mCompleted\e[0m"
  
  else
  
    echo "Aligning ${name} using Minimap2"
    
    # remove secondary and supplementary alignments, but keep duplicates (required for phasing)
    minimap2 -t 46 -ax splice --secondary=no -R "@RG\tID:${name}\tSM:${name}\tLB:lib1\tPL:ONT" ${GENOME_FASTA} $1 > $2/${name}.sam 2> $2/${name}_minimap2.log
    samtools view -h -F 2308 $2/${name}.sam > $2/${name}_filtered.sam

    # sort sam file 
    samtools sort -O SAM $2/${name}_filtered.sam > $2/${name}_filtered_sorted.sam  

    # convert to bam file
    samtools view -S -b $2/${name}_filtered.sam | samtools sort -o $2/${name}_filtered_sorted.bam
    samtools index $2/${name}_filtered_sorted.bam
  
    htsbox samview -pS $2/${name}_sorted.sam > $2/${name}.paf
    awk -F'\t' '{if ($6!="*") {print $0}}' $2/${name}.paf > $2/${name}.filtered.paf
    awk -F'\t' '{print $1,$6,$8+1,$2,$4-$3,($4-$3)/$2,$10,($10)/($4-$3),$5,$13,$15,$17}' $2/${name}.filtered.paf | sed -e s/"mm:i:"/""/g -e s/"in:i:"/""/g -e s/"dn:i:"/""/g | sed s/" "/"\t"/g > $2/${name}"_mappedstats.txt"
  
  fi

}

GENOME_FASTA=/lustre/projects/Research_Project-MRC148213/lsl693/references/human/hg38.fa
WKD_ROOT=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/Gina/analysis/EX165_PFC
run_minimap2 ${WKD_ROOT}/6_blaze/1_blaze/EX165_PFCmatched_reads.fastq.gz ${WKD_ROOT}/6_blaze/1_blaze/2_minimap