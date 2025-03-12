# multiome data analysis script for individual samples

# Script to analyse pre-processed individual multiome samples
# 1. Uploads pre-processed multiome data (i.e. filtered by cell QC metrics)
# 2. Normalise and reduce dimensions (RNA data)
# 3. Latent semantic indexing (LSI) (ATAC data)
# 4. Combine modalities by WNN analysis and reduce dimensions (UMAP)
#
# Arguments
# 1. file path to data directory (this should contain preprocessed data as Seurat object)
# 2. file path to results directory (created if doesn't exist)
# 3. sample id (a directory with this name is created within results directory)
#
# Outputs
# 1. Creates subdir per sample in results dir (and a figures subdir within there) for various graphs / txt outputs
# 2. Saves clustered data as RData object in sample results subdir


library(Signac)
library(Seurat)
# library(EnsDb.Hsapiens.v86)
# library(BSgenome.Hsapiens.UCSC.hg38)
library(magrittr)
library(ggplot2)

args = commandArgs(trailingOnly = T)


## ----------  set dir and files  --------------------

# set dir
dir_data = args[1]
# dir_data = "/lustre/projects/Research_Project-MRC190311/SingleCell/BMIAMP/3_preprocessed/11062_NP16_239_PFC/"
dir_out = args[2]
# dir_out = "/lustre/projects/Research_Project-MRC190311/SingleCell/BMIAMP/4_clustered/"
sample_id = args[3]
# sample_id = "11062_NP16_239_PFC"


dir_id = paste0(dir_out, sample_id, "/")
dir_fig = paste0(dir_id, "figures/")
if(!file.exists(dir_out)){dir.create(dir_out)}
if(!file.exists(dir_id)){dir.create(dir_id)}
if(!file.exists(dir_fig)){dir.create(dir_fig)}



## -------------- Load processed sample data ----------------

load(file = paste0(dir_data, "processed_data.RData"))



## -------------- Gene expression data processing ----------------

# standard scRNA processing steps (normalisation & PCA)

DefaultAssay(sc) = "RNA"
sc = NormalizeData(sc, verbose = F) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = F) %>%
  ScaleData(verbose = F) %>%
  RunPCA(reduction.name = "pca.rna") %>%
  RunUMAP(reduction = "pca.rna", dims = 1:20, 
          reduction.name = "umap.rna", reduction.key = "rnaUMAP_", verbose = F)

# elbow plot for RNA PCA
pdf(paste0(dir_fig, "rna-pca-elbow-plot.pdf"), width = 6, height = 6)
ElbowPlot(sc, reduction = "pca.rna")
invisible(dev.off())
  

## -------------- ATAC data processing ----------------

# standard scATAC processing (latent semantic indexing, LSI)

DefaultAssay(sc) = "ATAC"
sc = RunTFIDF(sc) %>%
  FindTopFeatures(min.cutoff = 'q0') %>%
  RunSVD() %>%
  RunUMAP(reduction = "lsi", dims = 2:50, 
          reduction.name = "umap.atac", reduction.key = "atacUMAP_")

# check correlation between LSI components and sequencing depth (1st component captures seq depth so remove from downstream analysis)
pdf(paste0(dir_fig, "atac-lsi-plot.pdf"), width = 6, height = 4)
DepthCor(sc)
invisible(dev.off())



## -------------- Multimodal weighted nearest neighbour analysis ----------------

sc = FindMultiModalNeighbors(sc, reduction.list = list("pca.rna", "lsi"),
                             dims.list = list(1:20, 2:50), verbose = F) %>%
  RunUMAP(nn.name = "weighted.nn", reduction.name = "umap.wnn", 
          reduction.key = "wnnUMAP_")

for (i in c(0.05, 0.1, 0.25, 0.5, 0.8)){
  sc = FindClusters(sc, graph.name = "wsnn", algorithm = 1, verbose = F, resolution = i)
}

# plot relative modality weights
pdf(paste0(dir_fig, "wsnn-weights.pdf"), width = 6, height = 4)
VlnPlot(sc, features = "RNA.weight", group.by = "wsnn_res.0.5", sort = T, pt.size = 0.1) + NoLegend()
invisible(dev.off())

## ----------  plot cluster resolutions UMAP  --------------------

pdf(paste0(dir_fig, "umap-cluster-resolution.pdf"), width = 6, height = 6)
print(DimPlot(sc, reduction = "umap.wnn", group.by = "wsnn_res.0.05", label = T))
print(DimPlot(sc, reduction = "umap.wnn", group.by = "wsnn_res.0.1", label = T))
print(DimPlot(sc, reduction = "umap.wnn", group.by = "wsnn_res.0.25", label = T))
print(DimPlot(sc, reduction = "umap.wnn", group.by = "wsnn_res.0.5", label = T))
print(DimPlot(sc, reduction = "umap.wnn", group.by = "wsnn_res.0.8", label = T))
dev.off()

# set cluster resolution to use
Idents(sc) = "wsnn_res.0.05"

# plot UMAPs and clustering results
pdf(paste0(dir_fig, "wsnn-umap.pdf"), width = 15, height = 5)
p1 <- DimPlot(sc, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE, pt.size = 0.4) + ggtitle("RNA")
p2 <- DimPlot(sc, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE, pt.size = 0.4) + ggtitle("ATAC")
p3 <- DimPlot(sc, reduction = "umap.wnn", label = TRUE, label.size = 2.5, repel = TRUE, pt.size = 0.4) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
invisible(dev.off())



## -------------- Save clustered sample data ----------------

save(sc, file = paste0(dir_id, "clustered_data.RData"))
