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


library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

args = commandArgs(trailingOnly = T)


## ----------  set dir and files  --------------------

# set dir
dir_data = args[1]
# dir_data = "/lustre/projects/Research_Project-MRC190311/SingleCell/BMIAMP/2_cellranger/11062_NP16_239_PFC/outs/"
dir_out = args[2]
# dir_out = "/lustre/projects/Research_Project-MRC190311/SingleCell/BMIAMP/3_preprocessed/"
sample_id = args[3]
# sample_id = "11062_NP16_239_PFC"
short_id = args[4]
# short_id = "PFC"

dir_id = paste0(dir_out, sample_id, "/")
dir_fig = paste0(dir_id, "figures/")
if(!file.exists(dir_out)){dir.create(dir_out)}
if(!file.exists(dir_id)){dir.create(dir_id)}
if(!file.exists(dir_fig)){dir.create(dir_fig)}

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
  project = short_id
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
sc@meta.data$id = short_id
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
pdf(paste0(dir_fig,'qc_pre_filter_vln.pdf'), width = 10, height = 5)
VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribosomal"), ncol = 4)
VlnPlot(sc, features = c("nFeature_ATAC", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"), ncol = 4)
invisible(dev.off())

pdf(paste0(dir_fig,'qc_pre_filter_plt.pdf'), width = 8, height = 4)
plot1 = FeatureScatter(sc, feature1 = "nCount_RNA", feature2 = "percent.mt") + NoLegend()
plot2 = FeatureScatter(sc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend()
plot1 + plot2
plot3 = FeatureScatter(sc, feature1 = "nCount_RNA", feature2 = "nCount_ATAC") + NoLegend()
plot4 = FeatureScatter(sc, feature1 = "nCount_ATAC", feature2 = "nFeature_ATAC") + NoLegend()
plot3 + plot4
invisible(dev.off())

sc$nucleosome_group = ifelse(sc$nucleosome_signal > 2, 'NS>2', 'NS<2')
sc$high.tss = ifelse(sc$TSS.enrichment > 3, 'TSS>3', 'TSS<3')

pdf(paste0(dir_fig, 'qc_atac.pdf'), width = 8, height = 4)
FragmentHistogram(sc, assay = "ATAC", group.by = 'nucleosome_group')
TSSPlot(sc, assay = "ATAC", group.by = 'high.tss') + NoLegend()
invisible(dev.off())


# record various QC metrics
sink(paste0(dir_id, "qc_metrics.txt"))
cat(paste0("Initial dimensions:\nRNA:\t\t", dim(sc$RNA)[1], " x ", dim(sc$RNA)[2]))
cat(paste0("\nATAC:\t\t", dim(sc$ATAC)[1], " x ", dim(sc$ATAC)[2]))
cat(paste0("\n\nNo. of cells with ...",
"\nnCount_RNA < 200:\t\t", length(which(sc$nCount_RNA < 200)),
"\nnCount_ATAC < 200:\t\t", length(which(sc$nCount_ATAC < 200)),
"\nTSS.enrichment < 1:\t\t", length(which(sc$TSS.enrichment < 1)),
"\nnucleosome_signal >2:\t\t", length(which(sc$nucleosome_signal > 2)),
"\npercent.mito > 5:\t\t", length(which(sc$percent.mt > 5)),
"\npercent.ribo > 5:\t\t", length(which(sc$percent.ribosomal > 5))
))
cat(paste0("\n\nNo. of features with ...",
           "\nRNA detected in 10+ cells:\t\t", length(which(rowSums(GetAssayData(sc, assay = "RNA", slot = "counts") != 0) >= 10)),
           "\nRNA detected in 50+ cells:\t\t", length(which(rowSums(GetAssayData(sc, assay = "RNA", slot = "counts") != 0) >= 50)),
           "\nRNA detected in 100+ cells:\t\t", length(which(rowSums(GetAssayData(sc, assay = "RNA", slot = "counts") != 0) >= 100)),
           "\nATAC detected in 10+ cells:\t\t", length(which(rowSums(GetAssayData(sc, assay = "ATAC", slot = "counts") != 0) >= 10)),
           "\nATAC detected in 50+ cells:\t\t", length(which(rowSums(GetAssayData(sc, assay = "ATAC", slot = "counts") != 0) >= 50)),
           "\nATAC detected in 100+ cells:\t\t", length(which(rowSums(GetAssayData(sc, assay = "ATAC", slot = "counts") != 0) >= 100))
))

# filter cells by QC metrics (very lenient - removes 38 cells)
sc = subset(
  x = sc,
  subset = 
    nCount_ATAC > 200 &              
    nCount_RNA > 200 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)

cat(paste0("\n\nFiltered dimensions:\nRNA:\t\t", dim(sc$RNA)[1], " x ", dim(sc$RNA)[2]))
cat(paste0("\nATAC:\t\t", dim(sc$ATAC)[1], " x ", dim(sc$ATAC)[2]))
sink()


# plot no. of features detected in X cells
pdf(paste0(dir_fig,'qc_features_plt.pdf'), width = 6, height = 6)
plot(c(1:dim(sc$RNA)[1]), sort(rowSums(GetAssayData(sc, assay = "RNA", slot = "counts"))), 
     main = "RNA features", xlab = "feature index", ylab = "number of cells", 
     pch = 16, cex = 0.5, log = "y", yaxt = "n", xlim = c(0,40000))
myTicks = axTicks(2)
axis(2, at = myTicks, labels = formatC(myTicks, format = "d"))
plot(c(1:dim(sc$ATAC)[1]), sort(rowSums(GetAssayData(sc, assay = "ATAC", slot = "counts"))), 
       main = "ATAC features", xlab = "feature index", ylab = "number of cells", pch = 16, cex = 0.5, log = "y", yaxt = "n")
myTicks = axTicks(2)
axis(2, at = myTicks, labels = formatC(myTicks, format = "d"))
invisible(dev.off())



## -------------- Save processed sample data ----------------

save(sc, file = paste0(dir_id, "processed_data.RData"))



