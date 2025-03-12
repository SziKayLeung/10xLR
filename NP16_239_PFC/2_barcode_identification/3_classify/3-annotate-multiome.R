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
# library(BSgenome.Hsapiens.UCSC.hg38)
library(magrittr)
library(ggplot2)
library(readxl)
library(patchwork)

args = commandArgs(trailingOnly = T)


## ----------  set dir and files  --------------------

# set dir
dir_proj = args[1]
# dir_proj = "/lustre/projects/Research_Project-MRC190311/SingleCell/BMIAMP/"
dir_out = args[2]
# dir_out = "/lustre/projects/Research_Project-MRC190311/SingleCell/BMIAMP/5_annotated/"
sample_id = args[3]
# sample_id = "11062_NP16_239_PFC"

dir_data = paste0(dir_proj, "4_clustered/", sample_id, "/")
# dir_data = "/lustre/projects/Research_Project-MRC190311/SingleCell/BMIAMP/4_clustered/11062_NP16_239_PFC/"


dir_id = paste0(dir_out, sample_id, "/")
dir_fig = paste0(dir_id, "figures/")
if(!file.exists(dir_out)){dir.create(dir_out)}
if(!file.exists(dir_id)){dir.create(dir_id)}
if(!file.exists(dir_fig)){dir.create(dir_fig)}



## -------------- Load sample data ----------------

load(file = paste0(dir_data, "clustered_data.RData"))

DefaultAssay(sc) = "RNA"

# import marker genes
marker_genes = read_excel(paste0(dir_proj, "0_metadata/consensus-marker-genes.xlsx"), sheet = "Markers")


# set cluster resolution to use
Idents(sc) = "wsnn_res.0.05"


## ----------  plot known marker genes UMAP  --------------------

pdf(paste0(dir_fig, "rna-umap-all-markers.pdf"), width = 14, height = 5)
ngenes = dim(marker_genes)[1]
nplots = ifelse(ngenes%%3 == 0, ngenes%/%3, ngenes%/%3+1)  # extra plot if not divisible by 3
for (i in 1:nplots){
  print(FeaturePlot(sc, features = as.vector(marker_genes$Marker[(3*i-2):(3*i)]),
                    order = T, min.cutoff = "q1", reduction = "umap.wnn", ncol=3))
}
dev.off()


## ----------  plot selected marker genes UMAP  --------------------

# import marker genes
marker_genes = read_excel(paste0(dir_proj, "0_metadata/consensus-marker-genes.xlsx"), sheet = "Select")

pdf(paste0(dir_fig, "rna-umap-markers.pdf"), width = 14, height = 5)
ngenes = dim(marker_genes)[1]
nplots = ifelse(ngenes%%3 == 0, ngenes%/%3, ngenes%/%3+1)  # extra plot if not divisible by 3
for (i in 1:nplots){
  print(FeaturePlot(sc, features = as.vector(marker_genes$Marker[(3*i-2):(3*i)]),
                    order = T, min.cutoff = "q1", reduction = "umap.wnn", ncol=3))
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
pdf(paste0(dir_fig, "rna-dotplot-markers.pdf"), width = 12, height = 5)
DotPlot(sc, features = marker_genes$Marker) + RotatedAxis()
dev.off()

# heatmaps
# Set raster = F to avoid blurring in preview view of pdf. No need to downsample as quite a small dataset anyway.
pdf(paste0(dir_fig, "rna-heatmap-markers.pdf"), width = 12, height = 12)
DoHeatmap(sc, features = marker_genes$Marker, size = 3, raster = F)
# DoHeatmap(subset(sc, downsample = 100), features = marker_genes$Marker, size = 3, raster = F)
dev.off()

pdf(paste0(dir_fig, "rna-heatmap-topvargenes.pdf"), width = 12, height = 12)
DoHeatmap(sc, features = VariableFeatures(sc)[1:100], size = 3, raster = F)
dev.off()


## ----------  plot ATAC data  --------------------

# quick plots of ATAC data for marker genes to see if looks like peaks coincide with RNA
pdf(paste0(dir_fig, "atac_coverage_markers.pdf"), width = 10, height = 6)
for (i in marker_genes$Marker){
print(CoveragePlot(sc, region = i, features = i, assay = "ATAC", expression.assay = "RNA", peaks = T))
}
dev.off()
