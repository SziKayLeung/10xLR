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
  library("scales")
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

# RNA type of transcripts from GENCODE
geneType <- read.csv(paste0(LOGEN,"0_utils/references/gencode.v40.annotation.geneannotation_corrected.csv"), header = T)


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

RNAType <- plot_RNAtype(class.files, geneType)
RNAType$plot

merge(as.data.frame(table(preFilteredRNAType$classfile$Class)),
      as.data.frame(table(RNAType$classfile$Class)), by = "Var1") %>% 
  mutate(dropOffRate = (Freq.y - Freq.x)/ Freq.y) %>%
  arrange(dropOffRate)


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

Excitatory <- lapply(marker_genes[marker_genes$Celltype == "Excitatory neurons","Marker"], function(x) plot_cellType_by_gene(x))
Astrocytes <- lapply(marker_genes[marker_genes$Celltype == "Astrocytes","Marker"], function(x) plot_cellType_by_gene(x))
Microglia <- lapply(marker_genes[marker_genes$Celltype == "Microglia","Marker"], function(x) plot_cellType_by_gene(x))
Inhibitory <- lapply(marker_genes[marker_genes$Celltype == "Inhibitory neurons","Marker"], function(x) plot_cellType_by_gene(x))
Oligodendrocytes <- lapply(marker_genes[marker_genes$Celltype == "Oligodendrocytes","Marker"], function(x) plot_cellType_by_gene(x))
plot_grid(plotlist = Excitatory)
plot_grid(plotlist = Astrocytes)
plot_grid(plotlist = Microglia)
plot_grid(plotlist = Inhibitory)
plot_grid(plotlist = Oligodendrocytes)

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
