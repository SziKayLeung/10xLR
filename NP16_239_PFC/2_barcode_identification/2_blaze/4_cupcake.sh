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

BLAZE_READS=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/Gina/analysis/6_blaze/1_blaze/NP16_239_PFCmatched_reads.fastq
ALIGNED_BAM=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/Gina/analysis/6_blaze/2_minimap2/NP16_239_PFCmatched_reads_filtered_sorted.bam
OUTPUT_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/Gina/analysis/6_blaze/4_cupcake
GENOME_FASTA=/lustre/projects/Research_Project-MRC148213/lsl693/references/human/hg38.fa
GENOME_GTF=/lustre/projects/Research_Project-MRC148213/lsl693/references/human/gencode.v40.annotation.gtf

cd ${OUTPUT_DIR}
#gunzip /lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/Gina/analysis/6_blaze/1_blaze/NP16_239_PFCmatched_reads.fastq.gz
#export TMPDIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/Gina/analysis/6_blaze/4_cupcake

#pbmm2 align --preset ISOSEQ ${GENOME_FASTA} ${BLAZE_READS} NP16_239_PFC_mapped.bam --log-level DEBUG --log-file NP16_239_PFC_mapped.log
samtools sort NP16_239_PFC_mapped.bam -o NP16_239_PFC_mapped_sorted.bam
isoseq3 collapse NP16_239_PFC_mapped_sorted.bam NP16_239_PFC_collapsed.gff \
      --min-aln-coverage 0.85 --min-aln-identity 0.95 --do-not-collapse-extra-5exons \
      --log-level TRACE --log-file NP16_239_PFC_collapsed.log

awk '{split($1, a, "#"); print a[2], $2}' NP16_239_PFC_collapsed.read_stat.txt > NP16_239_PFC_collapsed_split.read_stat.txt
awk '{split($1, a, "_"); print a[1], $2}' NP16_239_PFC_collapsed_split.read_stat.txt > NP16_239_PFC_collapsed_split2.read_stat.txt

export SOFTDIR=/lustre/projects/Research_Project-MRC148213/lsl693/software
export SEQUENCE=$CUPCAKE/sequence
export PYTHONPATH=$PYTHONPATH:$SEQUENCE
export SQANTI3_DIR=${SOFTDIR}/SQANTI3
export CAGE_PEAK=$SQANTI3_DIR/data/ref_TSS_annotation/human.refTSS_v3.1.hg38.bed
export POLYA=$SQANTI3_DIR/data/polyA_motifs/mouse_and_human.polyA_motif.txt
export SQANTI_JSON=/lustre/projects/Research_Project-MRC190311/scripts/sequencing/longReadseq/SQANTI3-5.1/SQANTI3-5.1/utilities/filter/filter_adapted.json

python $SQANTI3_DIR/sqanti3_qc.py -v
name=NP16_239_PFC
python $SQANTI3_DIR/sqanti3_qc.py NP16_239_PFC_collapsed.gff ${GENOME_GTF} ${GENOME_FASTA} \
--CAGE_peak ${CAGE_PEAK} \
--polyA_motif_list ${POLYA} --skipORF \
--genename --isoAnnotLite --report skip -t 30 &> ${name}.sqanti.qc.log

echo "Processing Sample ${name} for SQANTI filter"
name=NP16_239_PFC_collapsed
python $SQANTI3_DIR/sqanti3_filter.py rules ${name}"_classification.txt" --gtf ${name}"_corrected.gtf" -j=${SQANTI_JSON} --skip_report &> ${name}.sqanti.filter.log
