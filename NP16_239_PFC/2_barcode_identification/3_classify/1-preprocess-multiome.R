# multiome pre-processing script for individual samples

# Script to pre-process individual multiome samples before merging/integrating for downstream analysis
# 1. Uploads CellRanger output matrices & ATAC fragments data
# 2. QC on cells (Signac)
#
# Arguments
# 1. file path to data directory (this should contain raw and filtered CellRanger output matrices as .h5)
# 2. file path to results directory (created if doesn't exist)
# 3. sample id (a directory with this name is created within results directory)
#
# Outputs
# 1. Creates subdir per sample in results dir (and a figures subdir within there) for various graphs / txt outputs
# 2. Saves preprocessed data for each sample as RData object in sample results subdir
#
# TODO
# - add initial ambient RNA removal step?
# - add scDblFinder step?
# - improve peak calling with MACS2?
#options(timeout = 1000)
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")


library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)


dir_data = "C:/Users/sl693/OneDrive - University of Exeter/ExeterPostDoc/1_Projects/10XSingleCell/pilot_V0311/CellRanger/"
dir_out = "C:/Users/sl693/OneDrive - University of Exeter/ExeterPostDoc/1_Projects/10XSingleCell/pilot_V0311/CellClassification/"

## ----------  load data  -----------------------------

print("------------------- Loading data to Seurat ---------------------------")
# load filtered counts matrix and ATAC fragments
counts = Read10X_h5(paste0(dir_data, "filtered_feature_bc_matrix.h5"))
fragpath = paste0(dir_data, "atac_fragments.tsv.gz")

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86, verbose = F)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))


## ----------  create Seurat object  ----------------

# create seurat object with RNA
sc = CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA",
  project = "V0311"
)

# create ATAC assay and add
sc[["ATAC"]] = CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)


## ----------  Quality control  ----------------

print("------------------- Calculate QC metrics ---------------------------")

# set short ID for plotting
sc@meta.data$id = "V0311"
Idents(sc) = sc@meta.data$id

# calculate per-cell QC metrics using ATAC data
DefaultAssay(sc) = "ATAC"
sc = NucleosomeSignal(sc, verbose = F)       # calculates strength of nucleosome signal per cell (147-294bp : <147 bp fragments)
sc = TSSEnrichment(sc, fast = F, verbose = F)

# calculate per-cell QC metrics using RNA data
DefaultAssay(sc) = "RNA"
sc$percent.mt = PercentageFeatureSet(sc, pattern = "^MT-")
sc$percent.ribosomal = PercentageFeatureSet(sc, pattern = "^RP[LS]")

# Visualize QC metrics
pdf(paste0(dir_out,'qc_pre_filter.pdf'), width = 10, height = 5)
VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribosomal"), ncol = 4)
VlnPlot(sc, features = c("nFeature_ATAC", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"), ncol = 4)
plot1 = FeatureScatter(sc, feature1 = "nCount_RNA", feature2 = "percent.mt") + NoLegend()
plot2 = FeatureScatter(sc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend()
plot1 + plot2
plot3 = FeatureScatter(sc, feature1 = "nCount_RNA", feature2 = "nCount_ATAC") + NoLegend()
plot4 = FeatureScatter(sc, feature1 = "nCount_ATAC", feature2 = "nFeature_ATAC") + NoLegend()
plot3 + plot4
invisible(dev.off())


# record various QC metrics
sink(paste0(dir_out, "qc_metrics.txt"))
cat("No of cells with RNA reads: ", length(sc$nCount_RNA), "\n")
cat("No of cells with RNA reads < 200bp: ", length(which(sc$nCount_RNA < 200)), "\n")
cat("Number of cells with mitochondrial RNA > 5%: ",length(which(sc$percent.mt > 5)), "\n")
cat("Number of cells with ribosomal RNA > 5%: ",length(which(sc$percent.ribosomal > 5)), "\n")
# filter cells by QC metrics (very lenient - removes 38 cells)
sc = subset(
  x = sc,
  subset =           
    nCount_RNA > 200 
)
cat("Filtering cells with RNA reads < 200bp ")
cat("Filtered dimensions:", dim(sc$RNA)[1], " x ", dim(sc$RNA)[2])
sink()


## -------------- Save processed sample data ----------------

save(sc, file = paste0(dir_out, "processed_data.RData"))