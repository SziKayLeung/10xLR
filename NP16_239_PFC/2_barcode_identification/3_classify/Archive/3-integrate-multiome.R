# Multiome sample integration and clustering

# Script to upload filtered RNA data for samples, normalise, integrate and cluster
# 1. Uploads merged and filtered Seurat RNA data (note, maintains layer per sample)
# 2. Standard preprocessing per sample (normalise, variable features, scaling, PCA)
# 3. Plots UMAPs & clusters for merged (non-integrated) data
# 4. Integrates by Seurat's RCPA method and Harmony and plots results
#
# Arguments
# 1. file path to config.R (relative to submission dir: multiome-ann on ISCA)
# 2. project_id (i.e. name of directory with filtered data: <dir_qc>/<project_id>)
# 3. name of filtered dataset (e.g. "seurat_filtered".RData in directory above, loaded object should be seurat_filt)
#
# Outputs
# 1. Creates output directory (<dir_integrated>/<project_id>)
# 2. Plots UMAPs for merged & integrated data
# 3. Saves integrated data

# -----------------------------------------------------------------------
# LOAD PACKAGES AND SET PATHS
# -----------------------------------------------------------------------

# clear R env
rm(list = ls())

# packages
library(Seurat)
# library(magrittr)
library(ggplot2)
library(harmony)

# load filepaths and functions
args <- commandArgs(trailingOnly = T)
source(args[1])
# source("/lustre/projects/Research_Project-MRC190311/SingleCell/Scripts/multiome-ann/config/config.R")

source(paste0(dir_script,"/scripts/R/custom_functions.R"))

# set params
project_id <- args[2]
# project_id <- "V0335_and_V0341-counts"

# create output & figures directory
dir_out <- paste0(dir_integrated, "/", project_id)
if(!file.exists(dir_integrated)){dir.create(dir_integrated)}
if(!file.exists(dir_out)){dir.create(dir_out)}


# -----------------------------------------------------------------------
# LOAD DATA
# -----------------------------------------------------------------------

# filtered Seurat RNA data
load(paste0(dir_qc, "/", project_id, "/filtered/", args[3], ".RData"))
# load(paste0(dir_qc, "/", project_id, "/filtered/seurat_filtered.RData"))


# -----------------------------------------------------------------------
# STANDARD WORKFLOW (NO INTEGRATION)
# -----------------------------------------------------------------------

# each sample is separate layer in merged Seurat object so these steps are applied sample-wise
sc <- NormalizeData(seurat_filt)
sc <- FindVariableFeatures(sc)
sc <- ScaleData(sc)
sc <- RunPCA(sc, reduction.name = "pca.merged")

# temp save at intermediate stage (remove later)
save(sc, file = paste0(dir_out, "/sc_scaled.RData"))
# load(paste0(dir_out, "/sc_scaled.RData"))

pdf(paste0(dir_out, "/pca-elbow-plot-merged.pdf"), width = 6, height = 6)
ElbowPlot(sc, ndims = 50, reduction = "pca.merged") + ylim(0,NA)
invisible(dev.off())

# UMAP and clustering (Louvain for now)
sc <- FindNeighbors(sc, reduction = "pca.merged", dims = 1:30, assay = "RNA", verbose = F)
sc <- RunUMAP(sc, reduction = "pca.merged", dims = 1:30, assay = "RNA", reduction.name = "umap.merged", verbose = F)
sc <- FindClusters(sc, resolution = 0.5, cluster.name = "merged_clusters_0.5")     # 32 
sc <- FindClusters(sc, resolution = 0.1, cluster.name = "merged_clusters_0.1")     # 15
sc <- FindClusters(sc, resolution = 0.05, cluster.name = "merged_clusters_0.05")    #13

# plot merged data
png(paste0(dir_out, "/umap-merged-1.png"), width = 1400, height = 1200)
DimPlot(sc, reduction = "umap.merged", group.by = "sampleID") + labs(title = "Merged [sample ID]")
invisible(dev.off())
png(paste0(dir_out, "/umap-merged-2.png"), width = 1400, height = 1200)
DimPlot(sc, reduction = "umap.merged", group.by = "merged_clusters_0.1") + labs(title = "Merged [cluster resolution = 0.1]")
invisible(dev.off())


# -----------------------------------------------------------------------
# INTEGRATE DATA (RPCA & HARMONY)
# -----------------------------------------------------------------------

# Seurat's faster and more conservative algorithm (RPCA) ~ 30 min (need options setting to cope with 1.8G object - below command sets to 2G)
options(future.globals.maxSize = 2000 * 1024^2)
sc <- IntegrateLayers(object = sc, method = RPCAIntegration, orig.reduction = "pca.merged", 
                      new.reduction = "integrated.rpca", verbose = F)

save(sc, file = paste0(dir_out, "/sc_rpca.RData"))
file.remove(paste0(dir_out, "/sc_scaled.RData"))

# plot RPCA integration results
sc <- FindNeighbors(sc, reduction = "integrated.rpca", dims = 1:30)
sc <- FindClusters(sc, resolution = 0.5, cluster.name = "rpca_clusters_0.5")      
sc <- RunUMAP(sc, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca", verbose = F)
png(paste0(dir_out, "/umap-rpca-1.png"), width = 1400, height = 1200)
DimPlot(sc, reduction = "umap.rpca", group.by = "sampleID") + labs(title = "RPCA [sample ID]")
invisible(dev.off())
png(paste0(dir_out, "/umap-rpca-2.png"), width = 1400, height = 1200)
DimPlot(sc, reduction = "umap.rpca", group.by = "rpca_clusters_0.5") + labs(title = "RPCA [cluster resolution = 0.1]")
invisible(dev.off())

# Harmony
# load(paste0(dir_out, "/sc_rpca.RData"))
sc <- IntegrateLayers(object = sc, method = HarmonyIntegration, orig.reduction = "pca.merged", 
                      new.reduction = "integrated.harmony")

# plot harmony integration results
sc <- FindNeighbors(sc, reduction = "integrated.harmony", dims = 1:30)
sc <- FindClusters(sc, resolution = 0.5, cluster.name = "harmony_clusters_0.5")
sc <- RunUMAP(sc, reduction = "integrated.harmony", dims = 1:30, reduction.name = "umap.harmony", verbose = F)
png(paste0(dir_out, "/umap-harmony-1.png"), width = 1400, height = 1200)
DimPlot(sc, reduction = "umap.harmony", group.by = "sampleID") + labs(title = "Harmony [sample ID]")
invisible(dev.off())
png(paste0(dir_out, "/umap-harmony-2.png"), width = 1400, height = 1200)
DimPlot(sc, reduction = "umap.harmony", group.by = "harmony_clusters_0.5") + labs(title = "Harmony [cluster resolution = 0.1]")
invisible(dev.off())

save(sc, file = paste0(dir_out, "/sc_integrated.RData"))
file.remove(paste0(dir_out, "/sc_rpca.RData"))

