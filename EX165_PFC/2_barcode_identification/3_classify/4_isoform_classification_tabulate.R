#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
## Stats from manually merging barcodes, cell types and long-read collapsed isoforms in single cell dataset
## work on ISCA
## --------------------------------


## ---------- packages -----------------

suppressWarnings({
    suppressMessages(library("data.table"))
    suppressMessages(library("dplyr"))
})


## ---------- directories ----------------- 

dir = "/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/Gina/analysis/EX165_PFC/6_blaze/"
outputdir = "/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/Gina/analysis/EX165_PFC/6_blaze/5_postClassification/"


## ---------- load data ----------------- 

# barcodesAnno = output from Seurat of each barcode and corresponding cell type classified from short-read RNA-Seq data
barcodesAnno = fread(paste0(outputdir, "barcodesCellType.csv"), data.table = F)
# barcodes = output from BLAZE identifying barcodes from long-read ONT data
barcodes = fread(paste0(dir,"1_blaze/EX165_PFCputative_bc.csv"), data.table = F)
# readstat = output from Iso-Seq collapse
readstat <- fread(paste0(dir, "4_cupcake/EX165_PFC_PFC_collapsed_split2.read_stat.txt"))
colnames(readstat) <- c("read_id","pbid")

pdf(paste0(outputdir, "hist_barcode_quality.pdf"))
hist(barcodes$putative_bc_min_q)
dev.off()


## ---------- Filtering long-read sequencing data by barcode quality -----------------  

#1. Keeping reads that only contains AAAAAA and TTTT post UMI 
#2. Keeping reads with minimum barcode quality Q > 12 (94% accuracy, default threshold)
#3. Geeping reads that have GATCT, CTAGA in sequence before barcode

# 1. Containing only AAA and TTT post UMI
filteredbarcodes <- barcodes %>% filter(post_umi_flanking %in% c("AAAA", "AAAAA", "TTTTT", "TTTT"))

# keeping record of the number of barcodes with other sequences flanking post umi
post_umi_groups <- barcodes %>% group_by(post_umi_flanking) %>% tally() %>% as.data.frame()
write.csv(post_umi_groups, paste0(dir, "post_umi_groups.csv"), row.names = F)

# 2. Minimum barcode quality > 12 
filteredbarcodesQual <- filteredbarcodes %>% filter(putative_bc_min_q >= 12)

# plotting a distribution of barcode quality
pdf(paste0(dir, "hist_barcode_quality.pdf"))
hist(barcodes$putative_bc_min_q)
dev.off()

# 3. GATCT, CTAGA in sequence before barcode
longReadBarcodesFiltered <- filteredbarcodesQual %>% filter(pre_bc_flanking %in% c("GATCT", "CTAGA"))

# plotting a distribution of polyT end sequences
# not quite sure what metric does
pdf(paste0(dir, "polyT_end.pdf"))
hist(longReadBarcodesFiltered$polyT_end)
dev.off()

# save output
write.csv(longReadBarcodesFiltered, paste0(outputdir, "EX165_PFCputative_bc_filtered.csv"), row.names = F, quote = F)


## ---------- Merging long-read sequencing data with filtered barcodes to the cell types classified from Seurat -----------------  

# faster to merge as data.table
longReadBarcodesFiltered <- as.data.table(longReadBarcodesFiltered)
barcodesAnno<- as.data.table(barcodesAnno)
mergedFinal <- merge(filteredbarcodesQual, barcodesAnno, by.x = "putative_bc", by.y = "barcode")


## ---------- Collapsing cell-type specific long-read to isoforms -----------------  

# Using iso-seq collapsed read stat file which lists which ONT reads are collapsed
mergedFinalIsoform <- merge(readstat,  mergedFinal, by = "read_id", all.y = T)
write.csv(mergedFinalIsoform, paste0(outputdir, "EX165_PFCputative_bc_filtered_annotated_isoform.csv"), quote = F, row.names = F)

# keeping only the reads that are collapsed 
mergedFinalIsoformFiltered <- mergedFinalIsoform[!is.na(mergedFinalIsoform$pbid),]

# simplfy output to ease downstream loading
mergedFinalIsoformFilteredSimplifed <- mergedFinalIsoformFiltered %>% select(cellType, pbid)
write.csv(mergedFinalIsoformFilteredSimplifed, paste0(outputdir, "EX165_PFCputative_bc_filtered_annotated_isoform_simplifed.csv"), quote = F, row.names = F)


## ---------- Merging all barcodes with readstat ----------------- 

# note the number of reads match the abundance output file from Iso-Seq collapse
AllBarcodes <- merge(readstat, barcodes, by.x = "V1", by.y = "read_id")
write.csv(AllBarcodes, paste0(outputdir, "preFilteredbarcodes_isoform_merged.csv"), quote = F, row.names = F)