# multiome data analysis script for individual samples

# Script to check expression and annotate individual multiome samples
# 1. Uploads clustered multiome data (single sample)
# 2. Use known marker gene expression patterns to manually annotate
# 3. Prelim basic plots to check if ATAC informative
#
# Arguments
# 1. file path to project directory (project parent directory, e.g. BMIAMP)
# 2. file path to results directory (created if doesn't exist)
# 3. sample id (a directory with this name is created within results directory)
#
# Outputs
# 1. Creates subdir per sample in results dir (and a figures subdir within there) for various graphs 


library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(magrittr)
library(ggplot2)
library(readxl)
library(patchwork)
library(stringr)

dir_out = "C:/Users/sl693/OneDrive - University of Exeter/ExeterPostDoc/1_Projects/10XSingleCell/EX165_PFC/CellClassification/"
dir_script = "C:/Users/sl693/OneDrive - University of Exeter/ExeterPostDoc/2_Scripts/10xLR/"

## -------------- Load sample data ----------------

load(file = paste0(dir_out, "clustered_data.RData"))

DefaultAssay(sc) = "RNA"

# import marker genes
marker_genes = read.csv(paste0(dir_script, "utils/consensus-marker-genes.csv"))

# set cluster resolution to use
Idents(sc) = "wsnn_res.0.05"

## ----------  plot known marker genes UMAP  --------------------

pdf(paste0(dir_out, "rna-umap-all-markers.pdf"), width = 14, height = 5)
ngenes = dim(marker_genes)[1]
nplots = ifelse(ngenes%%3 == 0, ngenes%/%3, ngenes%/%3+1)  # extra plot if not divisible by 3
for (i in 1:nplots){
  print(FeaturePlot(sc, features = as.vector(marker_genes$Marker[(3*i-2):(3*i)]),
                    order = T, min.cutoff = "q1", reduction = "umap.wnn", ncol=3))
}
dev.off()

# select genes
select_genes <- marker_genes[marker_genes$Select == "TRUE",]
pdf(paste0(dir_out, "rna-umap-markers.pdf"), width = 14, height = 5)
ngenes = dim(select_genes)[1]
nplots = ifelse(ngenes%%3 == 0, ngenes%/%3, ngenes%/%3+1)  # extra plot if not divisible by 3
for (i in 1:nplots){
  p  <- FeaturePlot(sc, features = as.vector(select_genes$Marker[(3*i-2):(3*i)]),
                    order = T, min.cutoff = "q1", reduction = "umap.wnn", ncol=3)
  p <- p + plot_annotation(title = select_genes$Celltype[3*i-2])
  print(p)
}
dev.off()


## ----------  manually annotate cell types  --------------------

sc = RenameIdents(sc,
                  "0" = "oligodendrocytes",
                  "1" = "astrocytes",
                  "2" = "OPCs",
                  "3" = "inhibitory neurons",
                  "4" = "excitatory neurons",
                  "5" = "microglia",
                  "6" = "excitatory neurons")

sc$cell_type = Idents(sc)


## ----------  plot expression levels  --------------------

# N.B. oligo genes don't feature in scale.data (i.e. top 2000 variable features), can scale all as below to plot
# sc_scale_all = ScaleData(sc, assay = "RNA", features = rownames(sc[["RNA"]]))

# dotplots
pdf(paste0(dir_out, "rna-dotplot-markers.pdf"), width = 12, height = 5)
DotPlot(sc, features = marker_genes$Marker) + RotatedAxis()
dev.off()

# heatmaps
# Set raster = F to avoid blurring in preview view of pdf. No need to downsample as quite a small dataset anyway.
pdf(paste0(dir_out, "rna-heatmap-markers.pdf"), width = 12, height = 12)
DoHeatmap(sc, features = marker_genes$Marker, size = 3, raster = F)
# DoHeatmap(subset(sc, downsample = 100), features = marker_genes$Marker, size = 3, raster = F)
dev.off()

pdf(paste0(dir_out, "rna-heatmap-topvargenes.pdf"), width = 12, height = 12)
DoHeatmap(sc, features = VariableFeatures(sc)[1:100], size = 3, raster = F)
dev.off()


## ----------  plot ATAC data  --------------------

# quick plots of ATAC data for marker genes to see if looks like peaks coincide with RNA
pdf(paste0(dir_out, "atac_coverage_markers.pdf"), width = 10, height = 6)
for (i in marker_genes$Marker){
print(CoveragePlot(sc, region = i, features = i, assay = "ATAC", expression.assay = "RNA", peaks = T))
}
dev.off()

save(sc, file = paste0(dir_out, "annotated_data.RData"))


## -------------- output barcodes as txt file for merging with long read ----------------

# output barcodes and cell type
barcodesCellTypes <- cbind(as.data.frame(sc$cell_type), as.data.frame(sc$wsnn_res.0.05))
barcodesCellTypes <- barcodesCellTypes %>% tibble::rownames_to_column(., var = "barcode")
colnames(barcodesCellTypes) <- c("barcode","cellType","clusterNum_res0.05")
barcodesCellTypes$barcode <- word(barcodesCellTypes$barcode,c(1), sep = fixed("-"))
write.csv(barcodesCellTypes, paste0(dir_out, "barcodesCellType.csv"), row.names = F)
