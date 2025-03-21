#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
## Stats from manually merging barcodes, cell types and long-read collapsed isoforms in single cell dataset
## --------------------------------


## ---------- packages -----------------

suppressWarnings({
  suppressMessages(library("data.table"))
  suppressMessages(library("dplyr"))
  suppressMessages(library("ggplot2"))
  suppressMessages(library("ggVennDiagram"))
  suppressMessages(library("cowplot"))
  suppressMessages(library("scales"))
  suppressMessages(library("ggbio"))
  suppressMessages(library("GenomicRanges"))
})

# other function scripts
LOGEN <- c("C:/Users/sl693/OneDrive - University of Exeter/ExeterPostDoc/2_Scripts/LOGen/")
source(paste0(LOGEN,"transcriptome_stats/read_sq_classification.R"))

# plot the type of RNA
plot_RNAtype <- function(class_file, geneType){
  
  dat <- merge(class_file, distinct(geneType[,c("GeneSymbol","Class")]), by.x = "associated_gene", by.y = "GeneSymbol", all.x = T)
  
  p <- dat %>% group_by(Class) %>% tally() %>% 
    ggplot(., aes(y = reorder(Class,n), x = n)) + geom_bar(stat = "identity") +
    scale_x_log10(labels = label_comma()) +
    labs(x = "number of transcripts", y = "Class") + 
    theme_classic() 
  
  output <- list(p, dat)
  names(output) <- c("plot", "classfile")
  return(output)
}

## ---------- load files -----------------

# RNA type of transcripts from GENCODE
geneType <- read.csv(paste0(LOGEN,"0_utils/references/gencode.v40.annotation.geneannotation_corrected.csv"), header = T)

dir <- "C:/Users/sl693/OneDrive - University of Exeter/ExeterPostDoc/1_Projects/10XSingleCell/pilot_V0311/postCellClassification/"

# SQANTI classification file
class.files.names <- c(paste0(dir, "NP16_239_PFC_collapsed_RulesFilter_result_classification.txt"))
class.files <- SQANTI_class_preparation(class.files.names, "ns")    
preFilteredRNAType <- plot_RNAtype(class.files, geneType)
preFilteredRNAType$plot

# only keep FSM mono-exonic, i.e. remove all other categories of mono-exonic
monoexonicDiscard <- class.files %>% filter(exons == 1) %>% filter(structural_category != "FSM")
class.files <- class.files %>% filter(!isoform %in% monoexonicDiscard$isoform)

# barcode output file from 4_isoform_classification_tabulate.R 
# note there will be duplicates but this is because it's a simplified dataframe from the original annotated_isoform which includes the ONT read id
# the duplicates are when there are multiple ONT reads collapsed to the same isoform => abundance
barcode <- fread(paste0(dir, "NP16_239_PFCputative_bc_filtered_annotated_isoform_simplifed.csv"))

# import marker genes
dir_script = "C:/Users/sl693/OneDrive - University of Exeter/ExeterPostDoc/2_Scripts/10xLR/"
marker_genes = read.csv(paste0(dir_script, "utils/consensus-marker-genes.csv"))


## ---------- group the barcode by cell type and isoform  -----------------

# tally the number of ONT reads per isoform and classified by cellType
# keeping only the ones kept in SQANTI
barcode_cellType <- barcode %>% filter(pbid %in% class.files$isoform) %>% group_by(pbid, cellType) %>% tally()

# data wrangle 
abundance <- tidyr::spread(barcode_cellType, "cellType","n") %>% tibble::column_to_rownames(., var = "pbid")
abundance[is.na(abundance)] <- 0
colnames(abundance) <- paste0(colnames(abundance),"_FL")

# merge with abundance
class.files <- merge(class.files, abundance, by.x = "isoform", by.y = 0)

# filter by isoforms with less than total 2 reads across all cell populations
class.files <- class.files %>% mutate(FL = rowSums(select(., contains("_FL")), na.rm = TRUE))
filtered.class.files <- class.files %>% filter(FL >= 5)

# known genes
class.files.annoGenes <- filtered.class.files[!grepl("novel", filtered.class.files$associated_gene),]
table(class.files.annoGenes$polyA_motif_found)

# novel genes
novelGenes <- filtered.class.files[grepl("novel", filtered.class.files$associated_gene),] 
novelGenes %>% group_by(structural_category) %>% tally()
table(novelGenes$polyA_motif_found)
table(novelGenes$all_canonical)
ggplot(novelGenes, aes(x = exons, y = length, colour = structural_category, shape = polyA_motif_found)) + geom_point()

RNAType <- plot_RNAtype(filtered.class.files, geneType)
RNAType$plot

merge(as.data.frame(table(preFilteredRNAType$classfile$Class)),
      as.data.frame(table(RNAType$classfile$Class)), by = "Var1") %>% 
  mutate(dropOffRate = (Freq.y - Freq.x)/ Freq.y) %>%
  arrange(dropOffRate)

## ---------- summary statistics ------------

message("Number of reads kept from filtering mono-exonic isoforms:", sum(class.files$FL))
message("Number of reads kept from expression filtering (minimum 5 FL reads):", sum(filtered.class.files$FL))
message("Number of reads to known genes:", sum(class.files.annoGenes$FL))
message("Number of reads to novel genes:", sum(novelGenes$FL))

message("Number of isoforms:", nrow(filtered.class.files))
message("Number of protein coding isoforms:", nrow(RNAType$classfile %>% filter(Class == "protein_coding")))
message("number of genes with protein coding isoforms: ", length(unique(RNAType$classfile %>% filter(Class == "protein_coding") %>% .[,c("associated_gene")])))

message("Number of genes:", length(unique(filtered.class.files$associated_gene)))
message("Number of known genes:", length(unique(class.files.annoGenes$associated_gene)))


## ---------- plot per gene  -----------------

plot_cellType_by_gene <- function(gene){
  
  df <- filtered.class.files[filtered.class.files$associated_gene == gene,] 
  
  p <- df %>% select(contains("FL"),structural_category, isoform) %>% 
    select(-"FL") %>% 
    reshape2::melt(variable.name = "cellType", value.name = "count", id = c("structural_category","isoform")) %>% 
    mutate(cellType = word(cellType, c(1), sep = fixed("_"))) %>%
    ggplot(., aes(x = cellType, y = count, colour = structural_category)) + geom_boxplot() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(x = "Cell Type", y = "Number of FL reads", title = gene) 
  
  return(p)
}


`Excitatory neurons` <- lapply(marker_genes[marker_genes$Celltype == "Excitatory neurons","Marker"], function(x) plot_cellType_by_gene(x))
Astrocytes <- lapply(marker_genes[marker_genes$Celltype == "Astrocytes","Marker"], function(x) plot_cellType_by_gene(x))
Microglia <- lapply(marker_genes[marker_genes$Celltype == "Microglia","Marker"], function(x) plot_cellType_by_gene(x))
`Inhibitory neurons` <- lapply(marker_genes[marker_genes$Celltype == "Inhibitory neurons","Marker"], function(x) plot_cellType_by_gene(x))
Oligodendrocytes <- lapply(marker_genes[marker_genes$Celltype == "Oligodendrocytes","Marker"], function(x) plot_cellType_by_gene(x))
OPCs <- lapply(marker_genes[marker_genes$Celltype == "OPCs","Marker"], function(x) plot_cellType_by_gene(x))
Pericytes <- lapply(marker_genes[marker_genes$Celltype == "Pericytes","Marker"], function(x) plot_cellType_by_gene(x))

titles <- lapply(unique(marker_genes$Celltype), function(x) ggdraw() + draw_label(x, fontface = 'bold', size = 14))
names(titles) <- unique(marker_genes$Celltype)
pCellTypeByGenes <- lapply(unique(marker_genes$Celltype), function(x) plot_grid(titles[[x]], plot_grid(plotlist = get(x)), ncol = 1, rel_heights = c(0.1, 1)))

pdf(paste0(dir, "transcriptExpressionLongReadCellType.pdf"), width = 12, height = 10)
pCellTypeByGenes
dev.off()

## ---------- venn diagram of transcripts and genes across cell types  -----------------

cellType_transcript <- apply(abundance, 2, function(x) rownames(abundance)[which(x != 0)])
oligosAll <- unique(cellType_transcript$oligodendrocytes_FL, cellType_transcript$OPCs_FL)
neurons <- unique(cellType_transcript$`inhibitory neurons_FL`, cellType_transcript$`excitatory neurons_FL`)

# venn diagram of the common number of transcripts per cell type
final_cellType_transcript <- list(oligosAll, neurons, cellType_transcript$microglia_FL, cellType_transcript$astrocytes_FL)
names(final_cellType_transcript) <- c("Oligodendrocytes","Neurons","Microglia","Astrocytes")
ggVennDiagram(final_cellType_transcript, label = "count")

# venn diagram of the annotated genes
cellType_gene <- lapply(final_cellType_transcript, function(x)
  unique(class.files.annoGenes[class.files.annoGenes$isoform %in% x, "associated_gene"]))
ggVennDiagram(cellType_gene, label = "count")


## ---------- distribution of transcripts across genome-----------------

### generate gtf files on ISCA and transfer to onedrive
#write.table(novelGenes$isoform,paste0(dir,"unique_novelGenes_ID.csv"), quote = F, row.names = F, col.names = F)
#grep -f unique_novelGenes_ID.csv pilotPFC_collapsed_corrected.gtf > pilotPFC_novelGenes.filtered.gtf

#proteinCoding <- RNAType$classfile %>% filter(Class == "protein_coding")
#write.table(proteinCoding$isoform,paste0(dir,"unique_proteincoding_ID.csv"), quote = F, row.names = F, col.names = F)
#grep -f unique_proteincoding_ID.csv pilotPFC_collapsed_corrected.gtf > pilotPFC_proteinCodingGenes.filtered.gtf


plot_distribution_across_genome <- function(gtf){
  
  # read gtf
  data <- read.table(paste0(dir, "/", gtf), sep="\t", header=FALSE, stringsAsFactors=FALSE)
  colnames(data) <- c("chr", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
  
  # subset to transcripts only and those annotated to standard chromosomes
  transcripts <- subset(data, type == "transcript")
  chromosomes <- paste0("chr", c(1:22, "X","Y"))
  transcripts2chr <- transcripts %>% filter(chr %in% chromosomes)
  
  # Convert to GenomicRanges object
  gr <- GRanges(seqnames = transcripts2chr$chr,
                ranges = IRanges(start = transcripts2chr$start, end = transcripts2chr$end),
                strand = transcripts2chr$strand)
  ptranscriptAcrossGenome <- autoplot(gr) + ggtitle("Transcript Distribution Across Genome")
  
  ptranscriptsAcrossChr <- lapply(chromosomes, function(x) {
    chr_transcripts <- transcripts %>% filter(chr == x)
    
    # Only plot if there are transcripts for this chromosome
    if (nrow(chr_transcripts) > 0) {
      plot_transcripts_across_genome(chr_transcripts, x)
    } else {
      return(NULL)  # Pass NULL if no transcripts
    }
  })
  
  output <- list(ptranscriptAcrossGenome,ptranscriptsAcrossChr)
  return(output)
}

novelGenesGtf <- plot_distribution_across_genome("pilotPFC_novelGenes.filtered.gtf")

pdf(paste0(dir, "novelGenesAcrossChr.pdf"), width = 10, height = 6)
novelGenesGtf[[1]]
novelGenesGtf[[2]]
dev.off()

