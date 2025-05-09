#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
## Assess percentage of short coding transcripts in dataset
## --------------------------------

## ---------- packages -----------------

library("data.table")

LOGenDir <- "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/"
source(paste0(LOGenDir, "transcriptome_stats/read_sq_classification.R"))

## ---------- read files -----------------

dir <- "/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/Gina/analysis/EX165_PFC/6_blaze/4_cupcake/"

# sqanti
class.files <- SQANTI_class_preparation(paste0(dir, "EX165_PFC_PFC_collapsed_RulesFilter_result_classification.counts_filtered.txt"), "ns")

# ORFs
cpat <- read.table(paste0(dir, "EX165_PFC_PFC.ORF_prob.best.tsv"), header = T)
noORF <- read.table(paste0(dir, "EX165_PFC_PFC.no_ORF.txt"))

## ---------- subset short isoforms -----------------

# < 700bo
short_isoforms <- class.files[class.files$length < 700,]
nrow(short_isoforms)/nrow(class.files)

# subset cpat for just the short isoforms with ORF
cpat_short <- cpat %>% filter(seq_ID %in% short_isoforms)

# subset cat for the the short isoforms with ORF and coding (CPAT default human threshold)
cpat_short_coding <- cpat_short %>% filter(Coding_prob >= 0.364)

# check that all short isoforms covered = TRUE
setdiff(short_isoforms$isoforms, c(cpat_short$seq_ID, noORF$V1))


## ---------- stats -----------------

nrow(short_isoforms)
# 81339

nrow(cpat_short)
# 79519

nrow(cpat_short_coding)
# 31173

nrow(cpat_short_coding[cpat_short_coding$exons == 1,])
#487

table(cpat_short_coding$structural_category)
#          FSM           ISM           NIC           NNC Genic_Genomic
#         4449         14259           731         10489           241
#    Antisense        Fusion    Intergenic
#           83            84           823
           
cpat_short_coding <- merge(cpat_short_coding, short_isoforms[,c("isoform","length","exons","associated_transcript","structural_category"),], 
by.x = "seq_ID", by.y = "isoform")


## ---------- output to generate gtf -----------------

write.table(cpat_short_coding$seq_ID, paste0(dir,"cpat_short_coding_isoform.txt"), col.names = F, quote = F, row.names = F)
write.table(short_isoforms$isoform, paste0(dir,"all_short_coding_isoform.txt"), col.names = F, quote = F, row.names = F)

# bash command line
#grep -f cpat_short_coding_isoform.txt EX165_PFC_PFC_collapsed.filtered_counts_filtered.gtf > cpat_short_coding.gtf
#grep -f all_short_coding_isoform.txt EX165_PFC_PFC_collapsed.filtered_counts_filtered.gtf > all_short.gtf