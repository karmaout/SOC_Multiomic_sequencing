# Load necessary libraries
library(Seurat)
library(Signac)
library(EnsDb.Rnorvegicus.v79)
library(dplyr)
library(ggplot2)
library(patchwork)
library(readr)
library(GenomicRanges)
library(SummarizedExperiment)
library(ggrepel)
library(cowplot)
library(openxlsx)


# Load the aggregation data
aggr_data_pc <- read_csv("aggr.csv")

# Load the 10x hdf5 file
pc_data <- Read10X_h5("filtered_feature_bc_matrix.h5")

# Increase the maximum size allowed for future globals
options(future.globals.maxSize = 2 * 1024^3)  # 2 GB, adjust as needed

# Extract RNA and ATAC data
rna_counts <- pc_data$`Gene Expression`
atac_counts <- pc_data$Peaks

# Create Seurat object
PC <- CreateSeuratObject(counts = rna_counts)
PC[["percent.mt"]] <- PercentageFeatureSet(PC, pattern = "^Mt-")

# Add ATAC-seq data
granges <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
granges_use <- seqnames(granges) %in% standardChromosomes(granges)
atac_counts <- atac_counts[as.vector(granges_use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Rnorvegicus.v79)
seqlevelsStyle(annotations) <- "NCBI"
genome(annotations) <- "mRatBN7.2"

chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'mRatBN7.2',
  fragments = "atac_fragments.tsv.gz",
  min.cells = 10,
  annotation = annotations
)
PC[["ATAC"]] <- chrom_assay

# Extract numeric suffix from cell barcodes and create mapping
library_indices <- as.numeric(sapply(strsplit(colnames(PC), "-"), `[`, 2))
behavior_groups <- aggr_data_pc$behavior_group[library_indices]
library_ids <- aggr_data_pc$library_id[library_indices]

# Add metadata to Seurat object
metadata <- data.frame(behavior_group = behavior_groups, library_id = library_ids)
rownames(metadata) <- colnames(PC)
PC <- AddMetaData(PC, metadata = metadata)

# ATAC QC
DefaultAssay(PC) <- "ATAC"
PC <- NucleosomeSignal(PC)
PC <- TSSEnrichment(PC)
DensityScatter(PC, x = "nCount_ATAC", y = "TSS.enrichment", log_x = TRUE, quantiles = TRUE)

# Split Seurat object by library
library_list <- SplitObject(PC, split.by = "library_id")

# Perform RNA QC and generate violin plots
DefaultAssay(PC) <- "RNA"
plot_list_RNA <- lapply(seq_along(library_list), function(i) {
  VlnPlot(library_list[[i]], features = c("nCount_RNA", "percent.mt"), ncol = 2, log = TRUE, pt.size = 0) +
    NoLegend() +
    ggtitle(paste("Library", names(library_list)[i], "- Before Filtering"))
})
combined_QC_plot_RNA <- wrap_plots(plotlist = plot_list_RNA)
print(combined_QC_plot_RNA)

# Perform ATAC QC and generate violin plots
plot_list_ATAC <- lapply(seq_along(library_list), function(i) {
  VlnPlot(library_list[[i]], features = c("nCount_ATAC", "TSS.enrichment", "nucleosome_signal"), ncol = 3, log = TRUE, pt.size = 0) +
    NoLegend() +
    ggtitle(paste("Library", names(library_list)[i], "- Before Filtering"))
})
combined_QC_plot_ATAC <- wrap_plots(plotlist = plot_list_ATAC)
print(combined_QC_plot_ATAC)

# Mitochondria QC
VlnPlot(PC, features = c("nCount_ATAC", "nCount_RNA", "percent.mt"), ncol = 3, log = TRUE, pt.size = 0) + NoLegend()

# Filter low-quality cells
PC <- subset(
  x = PC,
  subset = nCount_ATAC < 100000 & 
    nCount_RNA < 100000 & 
    nCount_ATAC > 200 & 
    nCount_RNA > 200 & 
    percent.mt < 5
)
print(PC)

# RNA analysis
DefaultAssay(PC) <- "RNA"
PC <- SCTransform(PC, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
DefaultAssay(PC) <- "ATAC"
PC <- RunTFIDF(PC)
PC <- FindTopFeatures(PC, min.cutoff = 'q0')
PC <- RunSVD(PC)
PC <- RunUMAP(PC, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

# Calculate WNN graph  
PC <- FindMultiModalNeighbors(PC, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
PC <- RunUMAP(PC, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
PC <- FindClusters(PC, resolution = 0.1, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

# Save Seurat object
saveRDS(PC, file = "PC_neurons.rds")

# List of behavior groups and cell types
behavior_groups <- c("U_E", "P_E", "U_R", "P_R")
cell_types <- c("Inhibitory", "Excitatory")

# Assign each cell as cfos+ or cfos- based on group-specific thresholds
PC$cfos_status <- NA
for (bg in behavior_groups) {
  for (ct in cell_types) {
    subset_idx <- PC@meta.data$behavior_group == bg & PC@meta.data$cell_type == ct
    PC$cfos_status[subset_idx & PC@meta.data$IEG1 > 0] <- "cfos+"
    PC$cfos_status[subset_idx & PC@meta.data$IEG1 <= 0] <- "cfos-"
  }
}

# Define the desired order for the groups
desired_order <- c("U_E_Excitatory", "P_E_Excitatory", "U_E_Inhibitory", "P_E_Inhibitory",
                   "U_R_Excitatory", "P_R_Excitatory", "U_R_Inhibitory", "P_R_Inhibitory")

# Create a combined data frame for plotting
ieg_scores_df <- data.frame(
  IEG_Score = PC@meta.data$IEG1,
  Subset = paste(PC@meta.data$behavior_group, PC@meta.data$cell_type, sep = "_"),
  cfos_status = PC$cfos_status
)

# Filter out the "Other" cell type
ieg_scores_df <- subset(ieg_scores_df, !grepl("_Other$", Subset))

# Convert the Subset column to a factor with the specified order
ieg_scores_df$Subset <- factor(ieg_scores_df$Subset, levels = desired_order)

# Plot the violin plots for each subset
violin_plot <- ggplot(ieg_scores_df, aes(x = Subset, y = IEG_Score, fill = Subset)) +
  geom_violin() +
  geom_jitter(width = 0.2, alpha = 0.3) +
  labs(title = "IEG Expression Distribution in Different Subsets", x = "Subset", y = "IEG Score") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 24),
    axis.text.y = element_text(size = 24),
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    plot.title = element_text(size = 24),
    legend.text = element_text(size = 24),
    legend.title = element_text(size = 24)
  )

# Save the plot as a TIFF file
ggsave("PC_violin_plot.tiff", plot = violin_plot, width = 12, height = 8, dpi = 300, compression = "lzw")

# Print the plot to the console
print(violin_plot)

# Function to calculate the percentage of cfos+ neurons in each subset
calculate_cfos_percentage <- function(df) {
  total_cells <- nrow(df)
  cfos_positive_cells <- sum(df$cfos_status == "cfos+")
  percentage_cfos_positive <- (cfos_positive_cells / total_cells) * 100
  return(percentage_cfos_positive)
}

# Add library_id to ieg_scores_df based on matching cell identifiers
ieg_scores_df$library_id <- PC@meta.data[rownames(ieg_scores_df), "library_id"]

# Filter out NA values in cfos_status
ieg_scores_df <- ieg_scores_df[!is.na(ieg_scores_df$cfos_status), ]

# Split the data by Subset and Library combination
cfos_percentages_list <- split(ieg_scores_df, list(ieg_scores_df$Subset, ieg_scores_df$library_id))

# Remove empty subsets
cfos_percentages_list <- cfos_percentages_list[sapply(cfos_percentages_list, nrow) > 0]

# Calculate the cfos+ percentages for each group and library
cfos_percentages <- sapply(cfos_percentages_list, calculate_cfos_percentage)

# Convert the result into a data frame
cfos_percentages_raw_df <- data.frame(
  Subset_Library = names(cfos_percentages),
  Percentage = cfos_percentages
)

# Separate Subset and Library columns
cfos_percentages_raw_df <- transform(cfos_percentages_raw_df,
                                     Subset = sapply(strsplit(as.character(Subset_Library), "\."), `[`, 1),
                                     Library = sapply(strsplit(as.character(Subset_Library), "\."), `[`, 2))

# Save the raw data to a CSV file
write.csv(cfos_percentages_raw_df, "cfos_percentages_raw_data.csv", row.names = FALSE)

# Print the raw data to the console
print(cfos_percentages_raw_df)
