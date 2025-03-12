# Multiome clustering and marker gene expression analysis

# Script to upload integrated RNA data
# 1. Uploads integrated Seurat RNA data (note, still has 1 layer per sample; use harmony integration)
# 2. Test Leiden clustering with different resolutions
# 3. Plot expression data for some example genes
#
# Arguments
# 1. file path to config.R (relative to submission dir: multiome-ann on ISCA)
# 2. project_id (i.e. name of directory with integrated data: <dir_integrated>/<project_id>)
#               (N.B. Should contain sc_integrated.RData to load in for this analysis, name ")
# 3. file path to marker gene info file
#
# Outputs
# 1. Creates output directory (<dir_clustered>/<project_id>)
# 2. Plots UMAPs for clustered data

# -----------------------------------------------------------------------
# LOAD PACKAGES AND SET PATHS
# -----------------------------------------------------------------------

# clear R env
rm(list = ls())

# packages
library(Seurat)
library(ggplot2)
library(reticulate)         # required for Leiden clustering

# set python version for reticulate (so uses conda env with required python packages incl. leidenalg)
use_condaenv("singlecell")

# for trial runs via srun only-----------------------------------------------
# source("/lustre/projects/Research_Project-MRC190311/SingleCell/Scripts/multiome-ann/config/config.R")
# project_id <- "V0335_and_V0341-counts"
# source(paste0(dir_script,"/scripts/R/custom_functions.R"))

# for final script via sbatch -----------------------------------------------
# load filepaths and functions
args <- commandArgs(trailingOnly = T)
source(args[1])
source(paste0(dir_script,"/scripts/R/custom_functions.R"))

# set params
project_id <- args[2]
cluster_resolution <- "harmony_leiden_0.1"
file_markergenes <- paste0(dir_meta, "/consensus-marker-genes.csv")

# for all -------------------------------------------------------------------
# create output & figures directory
dir_out <- paste0(dir_clustered, "/", project_id)
if(!file.exists(dir_clustered)){dir.create(dir_clustered)}
if(!file.exists(dir_out)){dir.create(dir_out)}


# # -----------------------------------------------------------------------
# # LOAD DATA & TEST CLUSTERING RESOLUTION
# # -----------------------------------------------------------------------
# 
# # load integrated RNA data
# load(paste0(dir_integrated, "/", project_id, "/sc_integrated.RData"))
# 
# # test clustering resolution & plot UMAPs
# sc <- FindNeighbors(sc, reduction = "integrated.harmony", dims = 1:30)    # ensure stored NN is using harmony data
# res_params <- c(0.01, 0.1, 0.15, 0.2)
# for (i in res_params){
#   sc = FindClusters(sc, resolution = i, algorithm = 4, method = "igraph", cluster.name = paste0("harmony_leiden_", i))
#   png(paste0(dir_out, "/umap-harmony-leiden-", i, ".png"), width = 1400, height = 1200)
#   print(DimPlot(sc, reduction = "umap.harmony", group.by = paste0("harmony_leiden_", i)) + 
#             labs(title = paste0("Harmony [leiden ", i, "]")))
#   invisible(dev.off())
# }
# 
# # save clustered data
# save(sc, file = paste0(dir_out, "/sc_clustered.RData"))


# -----------------------------------------------------------------------
# LOAD DATA & PLOT CLUSTER EXPRESSION
# -----------------------------------------------------------------------

load(file = paste0(dir_out, "/sc_clustered.RData"))

# set clustering resolution
Idents(sc) <- cluster_resolution
sc$clusters <- Idents(sc)

# use RNA slot for expression data and join layers for expression analysis
DefaultAssay(sc) = "RNA"
sc <- JoinLayers(sc)

save(sc, file = paste0(dir_out, "sc_joined.RData"))

# # import marker genes & subset for those detected
# marker_genes_all = read.csv(file_markergenes, stringsAsFactors = F)
# colnames(marker_genes_all) = c("celltype", "gene")
# marker_genes = marker_genes_all[marker_genes_all$gene %in% rownames(sc),]
# 
# 
# 
# # check % cells expression each marker and max average expression in a cluster
# marker_genes$pc.cells = sapply(marker_genes$gene, 
#                                function (x) {sum(LayerData(sc, assay = "RNA", layer = "data")[x,]>0) / dim(sc)[2]})
# marker_genes$mean.exp = sapply(marker_genes$gene, 
#                                function (x) {max(AverageExpression(sc, features = x, assays = "RNA")$RNA)})
# 
# # remove rarely expressed genes (<1% cells) & order by cell type
# # marker_genes = marker_genes[which(pc.cells > 0.01),]
# marker_genes = marker_genes[order(marker_genes$celltype),]


