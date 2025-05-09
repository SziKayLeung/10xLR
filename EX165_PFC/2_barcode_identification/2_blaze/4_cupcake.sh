#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=60:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --mem=200G # specify bytes of memory to reserve
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=s.k.leung@exeter.ac.uk # email address

module load Miniconda2
source activate lrp

BLAZE_READS=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/Gina/analysis/EX165_PFC/6_blaze/1_blaze/EX165_PFCmatched_reads.fastq
ALIGNED_BAM=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/Gina/analysis/EX165_PFC/6_blaze/1_blaze/2_minimap/EX165_PFCmatched_reads_filtered_sorted.bam
OUTPUT_DIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/Gina/analysis/EX165_PFC/6_blaze/4_cupcake
GENOME_FASTA=/lustre/projects/Research_Project-MRC148213/lsl693/references/human/hg38.fa
GENOME_GTF=/lustre/projects/Research_Project-MRC148213/lsl693/references/human/gencode.v40.annotation.gtf

cd ${OUTPUT_DIR}
gunzip /lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/Gina/analysis/6_blaze/1_blaze/EX165_PFC_PFCmatched_reads.fastq.gz
export TMPDIR=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/Gina/analysis/EX165_PFC/6_blaze/4_cupcake

pbmm2 align --preset ISOSEQ ${GENOME_FASTA} ${BLAZE_READS} EX165_PFC_PFC_mapped.bam --log-level DEBUG --log-file EX165_PFC_PFC_mapped.log
samtools sort EX165_PFC_PFC_mapped.bam -o EX165_PFC_PFC_mapped_sorted.bam
isoseq3 collapse EX165_PFC_PFC_mapped_sorted.bam EX165_PFC_PFC_collapsed.gff \
      --min-aln-coverage 0.85 --min-aln-identity 0.95 --do-not-collapse-extra-5exons \
      --log-level TRACE --log-file EX165_PFC_PFC_collapsed.log

awk '{split($1, a, "#"); print a[2], $2}' EX165_PFC_PFC_collapsed.read_stat.txt > EX165_PFC_PFC_collapsed_split.read_stat.txt
awk '{split($1, a, "_"); print a[1], $2}' EX165_PFC_PFC_collapsed_split.read_stat.txt > EX165_PFC_PFC_collapsed_split2.read_stat.txt

export SOFTDIR=/lustre/projects/Research_Project-MRC148213/lsl693/software
export SEQUENCE=$CUPCAKE/sequence
export PYTHONPATH=$PYTHONPATH:$SEQUENCE
export SQANTI3_DIR=${SOFTDIR}/SQANTI3
export CAGE_PEAK=$SQANTI3_DIR/data/ref_TSS_annotation/human.refTSS_v3.1.hg38.bed
export POLYA=$SQANTI3_DIR/data/polyA_motifs/mouse_and_human.polyA_motif.txt
export SQANTI_JSON=/lustre/projects/Research_Project-MRC190311/scripts/sequencing/longReadseq/SQANTI3-5.1/SQANTI3-5.1/utilities/filter/filter_adapted.json

python $SQANTI3_DIR/sqanti3_qc.py -v
name=EX165_PFC_PFC
python $SQANTI3_DIR/sqanti3_qc.py EX165_PFC_PFC_collapsed.gff ${GENOME_GTF} ${GENOME_FASTA} \
--CAGE_peak ${CAGE_PEAK} \
--polyA_motif_list ${POLYA} --skipORF \
--genename --isoAnnotLite --report skip -t 30 &> ${name}.sqanti.qc.log

echo "Processing Sample ${name} for SQANTI filter"
name=EX165_PFC_PFC_collapsed
python $SQANTI3_DIR/sqanti3_filter.py rules ${name}"_classification.txt" --gtf ${name}"_corrected.gtf" -j=${SQANTI_JSON} --skip_report &> ${name}.sqanti.filter.log

##-------------------------------------------------------------------------

export LOGEN_ROOT=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen
export PATH=$PATH:${LOGEN_ROOT}/transcriptome_stats
export PATH=$PATH:${LOGEN_ROOT}/compare_datasets
export PATH=$PATH:${LOGEN_ROOT}/target_gene_annotation
export PATH=$PATH:${LOGEN_ROOT}/merge_characterise_dataset
export PATH=$PATH:${LOGEN_ROOT}/miscellaneous 
export PATH=$PATH:${LOGEN_ROOT}/assist_ont_processing
name=EX165_PFC_PFC

# LOGEN: subset cupcake classification file 
# merge cupcake classification file with abundance
# filter cupcake classification file with minimum number of reads and counts
subset_quantify_filter_tgenes.R \
--classfile ${OUTPUT_DIR}/${name}_collapsed_RulesFilter_result_classification.txt \
--expression ${OUTPUT_DIR}/EX165_PFC_PFC_collapsed.flnc_count.txt \
--filter --nsample=1 --nreads=10 --monoexonic=2 --target_genes=NA

# working variables
finalanno=${OUTPUT_DIR}/${name}_collapsed_RulesFilter_result_classification.counts_filtered.txt 
finaliso=${OUTPUT_DIR}/${name}_collapsed_RulesFilter_result_classification.filtered_isoforms.txt

# LOGEN: subset fasta and gtf using the finalised list of target gene isoforms
subset_fasta_gtf.py --gtf ${OUTPUT_DIR}/${name}_collapsed.filtered.gtf -i ${finaliso} -o counts_filtered
subset_fasta_gtf.py --fa ${OUTPUT_DIR}/${name}_collapsed_corrected.fasta -i ${finaliso} -o counts_filtered

# run_cpat <input_fasta> <output_name> <output_dir>
export HEXAMER=/lustre/projects/Research_Project-MRC148213/lsl693/references/CPAT/Human_Hexamer.tsv
export LOGITMODEL=/lustre/projects/Research_Project-MRC148213/lsl693/references/CPAT/Human_logitModel.RData
cpat.py --version
cpat.py -x ${HEXAMER} -d ${LOGITMODEL} -g ${OUTPUT_DIR}/${name}_collapsed_corrected_counts_filtered.fa --min-orf=50 --top-orf=50 -o ${name} 2> ${name}"_cpat.e"

# extract_best_orf <sample> <root_dir>
extract_fasta_bestorf.py --fa ${name}".ORF_seqs.fa" --orf ${name}".ORF_prob.best.tsv" --o_name ${name}"_bestORF" --o_dir ${OUTPUT_DIR} &> orfextract.log