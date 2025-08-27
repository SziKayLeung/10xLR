#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=144:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --mem=200G # specify bytes of memory to reserve
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=s.k.leung@exeter.ac.uk # email address

module load Miniconda2
source activate lrp

BLAZE_READS=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/Gina/analysis/EX167_Hypo/6_blaze/1_blaze/EX167_HYPmatched_reads_trimmed.fastq
ALIGNED_BAM=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/Gina/analysis/EX167_Hypo/6_blaze/2_minimap/EX167_HYPmatched_reads_trimmed.fastq_filtered_sorted.sam
OUTPUT_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/Gina/analysis/EX167_Hypo/6_blaze/4_cupcake
GENOME_FASTA=/lustre/projects/Research_Project-MRC148213/lsl693/references/human/hg38.fa
GENOME_GTF=/lustre/projects/Research_Project-MRC148213/lsl693/references/human/gencode.v40.annotation.gtf

mkdir -p ${OUTPUT_DIR}
cd ${OUTPUT_DIR}
gunzip /lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/Gina/analysis/6_blaze/1_blaze/EX167_Hypmatched_reads.fastq.gz
export TMPDIR=${OUTPUT_DIR}

pbmm2 align --preset ISOSEQ ${GENOME_FASTA} ${BLAZE_READS} EX167_Hyp_mapped.bam --log-level DEBUG --log-file EX167_Hyp_mapped.log
samtools sort EX167_Hyp_mapped.bam -o EX167_Hyp_mapped_sorted.bam
isoseq3 collapse EX167_Hyp_mapped_sorted.bam EX167_Hyp_collapsed.gff \
      --min-aln-coverage 0.85 --min-aln-identity 0.95 --do-not-collapse-extra-5exons \
      --log-level TRACE --log-file EX167_Hyp_collapsed.log

awk '{split($1, a, "#"); print a[2], $2}' EX167_Hyp_collapsed.read_stat.txt > EX167_Hyp_collapsed_split.read_stat.txt
awk '{split($1, a, "_"); print a[1], $2}' EX167_Hyp_collapsed_split.read_stat.txt > EX167_Hyp_collapsed_split2.read_stat.txt

export SOFTDIR=/lustre/projects/Research_Project-MRC148213/lsl693/software
export SEQUENCE=$CUPCAKE/sequence
export PYTHONPATH=$PYTHONPATH:$SEQUENCE
export SQANTI3_DIR=${SOFTDIR}/SQANTI3
export CAGE_PEAK=$SQANTI3_DIR/data/ref_TSS_annotation/human.refTSS_v3.1.hg38.bed
export POLYA=$SQANTI3_DIR/data/polyA_motifs/mouse_and_human.polyA_motif.txt
export SQANTI_JSON=/lustre/projects/Research_Project-MRC190311/scripts/sequencing/longReadseq/SQANTI3-5.1/SQANTI3-5.1/utilities/filter/filter_adapted.json

# sqanti3
source activate sqanti3
NAME=EX167_Hyp
python $SQANTI3_DIR/sqanti3_qc.py  --isoforms ${NAME}_collapsed.gff --refGTF $GENOME_GTF --refFasta $GENOME_FASTA --output ${NAME} --CAGE_peak $CAGE_PEAK --polyA_motif_list $POLYA --skipORF --report skip --dir ${OUTPUT_DIR} &> ${NAME}.sqanti.qc.log

SQANTI_JSON=/lustre/projects/Research_Project-MRC190311/scripts/sequencing/longReadseq/SQANTI3-5.1/SQANTI3-5.1/utilities/filter/filter_default.json
python $SQANTI3_DIR/sqanti3_filter.py rules --sqanti_class ${NAME}"_classification.txt" --filter_gtf ${NAME}"_corrected.gtf" --dir ${OUTPUT_DIR}  -j=${SQANTI_JSON} &> ${NAME}.sqanti.filter.log

