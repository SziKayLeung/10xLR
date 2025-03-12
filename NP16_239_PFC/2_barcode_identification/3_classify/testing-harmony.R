# Script to test Harmony parameters


library(Seurat)
library(harmony)

source("config/config.R")
source(paste0(dir_script,"/scripts/R/functions_qc.R"))

# set params
args <- commandArgs(trailingOnly = T)
data_id <- args[1]
filter_id <- args[2]

# # for testing as interactive script only
source("/lustre/projects/Research_Project-MRC190311/SingleCell/Scripts/multiome-ann/config/config.R")
source(paste0(dir_script,"/scripts/R/functions_qc.R"))
data_id <- "cellranger"
filter_id <- "filter3"

# create output dirs
dir_in <- paste0(dir_integrated, "/", data_id, "/", filter_id)
dir_out <- paste0(dir_integrated, "/testing_harmony")
if(!file.exists(dir_out)){dir.create(dir_out, recursive = T)}

load data and select filtered cells
load(paste0(dir_qc, "/", data_id, "/seurat_merged.RData"))
load(file = paste0(dir_in, "/filtered_cell_ids_", filter_id, ".RData"))

# subset data
sc <- subset(seurat_merged, cells = cell_ids)
# 
# # # small sample
# # sc_small <- subset(seurat_merged, subset = nCount_RNA > 2000)
# # sc_small <- NormalizeData(sc_small)
# # save(sc_small, file = paste0(dir_out, "/sc_small_norm.RData"))
# # n1 <- 2000
# # sc_small <- FindVariableFeatures(sc_small, nfeatures = n1)
# # sc_small <- ScaleData(sc_small)
# # sc_small <- RunPCA(sc, reduction.name = "pca.merged")
# # save(sc_small, file = paste0(dir_out, "/sc_small_pca_", n1, ".RData"))
# # 
# # # each sample is separate layer in merged Seurat object so these steps are applied sample-wise
sc <- NormalizeData(sc)
# save(sc, file = paste0(dir_out, "/sc_norm.RData"))
n2 <- 4000
sc <- FindVariableFeatures(sc, nfeatures = n2)
sc <- ScaleData(sc)
sc <- RunPCA(sc, reduction.name = "pca.merged")
#
# save(sc, file = paste0(dir_out, "/sc_pca_", n2, ".RData"))
# load(file = paste0(dir_out, "/sc_pca_", n2, ".RData"))

png(paste0(dir_out, "/harmony_seurat.png"), height = 800, width = 600)
# sc <- RunHarmony(sc, "orig.ident", plot_convergence = TRUE, reduction = "pca.merged")
IntegrateLayers(object = sc, method = HarmonyIntegration,
                      orig.reduction = "pca.merged", new.reduction = "integrated.harmony",
                      plot_convergence = TRUE)
dev.off()
