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

# Set base directory (use relative paths or custom paths)
base_dir <- "path/to/data/"  # Update with relative or custom path

# Load the aggregation file
aggr_file <- file.path(base_dir, "aggr.csv")
aggr_data <- read_csv(aggr_file)
print(head(aggr_data))

# Load the 10x hdf5 data (both RNA and ATAC)
h5_file <- file.path(base_dir, "filtered_feature_bc_matrix.h5")
BLA_data <- Read10X_h5(h5_file)

# Increase maximum allowed size for globals
options(future.globals.maxSize = 2 * 1024^3)  # 2 GB

# Extract RNA and ATAC data
rna_counts <- BLA_data$`Gene Expression`
atac_counts <- BLA_data$Peaks

# Create Seurat object for RNA data
BLA <- CreateSeuratObject(counts = rna_counts)
BLA[["percent.mt"]] <- PercentageFeatureSet(BLA, pattern = "^Mt-")

# Add ATAC-seq data
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Rnorvegicus.v79)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "mRatBN7.2"
frag.file <- file.path(base_dir, "atac_fragments.tsv.gz")
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'mRatBN7.2',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)
BLA[["ATAC"]] <- chrom_assay

# Add behavior groups and library IDs from aggregation file to metadata
library_indices <- as.numeric(sapply(strsplit(colnames(BLA), "-"), `[`, 2))
behavior_groups <- aggr_data$behavior_group[library_indices]
library_ids <- aggr_data$library_id[library_indices]
metadata <- data.frame(behavior_group = behavior_groups, library_id = library_ids)
rownames(metadata) <- colnames(BLA)
BLA <- AddMetaData(BLA, metadata = metadata)

# Perform ATAC QC
DefaultAssay(BLA) <- "ATAC"
BLA <- NucleosomeSignal(BLA)
BLA <- TSSEnrichment(BLA)
DensityScatter(BLA, x = "nCount_ATAC", y = "TSS.enrichment", log_x = TRUE, quantiles = TRUE)

# RNA QC - Violin plots for each library
library_list <- SplitObject(BLA, split.by = "library_id")
plot_list_RNA <- lapply(library_list, function(obj) {
  VlnPlot(obj, features = c("nCount_RNA", "percent.mt"), ncol = 2, log = TRUE, pt.size = 0) + NoLegend() +
    ggtitle(paste("Library", unique(obj$library_id), "- Before Filtering"))
})
combined_QC_plot_RNA <- wrap_plots(plotlist = plot_list_RNA)
print(combined_QC_plot_RNA)

# ATAC QC - Violin plots for each library
plot_list_ATAC <- lapply(library_list, function(obj) {
  VlnPlot(obj, features = c("nCount_ATAC", "TSS.enrichment", "nucleosome_signal"), ncol = 3, log = TRUE, pt.size = 0) + NoLegend() +
    ggtitle(paste("Library", unique(obj$library_id), "- Before Filtering"))
})
combined_QC_plot_ATAC <- wrap_plots(plotlist = plot_list_ATAC)
print(combined_QC_plot_ATAC)

# Mitochondrial QC
VlnPlot(BLA, features = c("nCount_ATAC", "nCount_RNA", "percent.mt"), ncol = 3, log = TRUE, pt.size = 0) + NoLegend()

# Filter out low-quality cells
BLA <- subset(
  x = BLA,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 100000 &
    nCount_ATAC > 200 &
    nCount_RNA > 200 &
    percent.mt < 5
)

# RNA analysis
DefaultAssay(BLA) <- "RNA"
BLA <- SCTransform(BLA, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
DefaultAssay(BLA) <- "ATAC"
BLA <- RunTFIDF(BLA) %>% FindTopFeatures(min.cutoff = 'q0') %>% RunSVD() %>% RunUMAP(reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

# Calculate WNN graph and UMAP embedding
BLA <- FindMultiModalNeighbors(BLA, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
BLA <- RunUMAP(BLA, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
BLA <- FindClusters(BLA, resolution = 0.05, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

# Plot WNN UMAP with cluster labels
p3 <- DimPlot(BLA, group.by = "seurat_clusters", reduction = "wnn.umap", label = TRUE, pt.size = 1) +
  theme(
    legend.position = "right",
    text = element_text(size = 24),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 24),
    legend.text = element_text(size = 24),
    legend.title = element_text(size = 24),
    plot.title = element_text(size = 24)
  )
ggsave("BLA_Cells_digits.tiff", plot = p3, width = 10, height = 10, dpi = 300, units = "in")

# Define and visualize marker genes
DefaultAssay(BLA) <- 'SCT'
Feature_plot <- FeaturePlot(BLA, features = c("Slc17a7", "Gad1", "Gad2", "Slc32a1", "Mog", "Cx3cr1", "Slco1a4", "Cspg4", "Slc1a2"), reduction = "wnn.umap")
ggsave("BLA_feature_plot.tiff", plot = Feature_plot, width = 12, height = 10, dpi = 300, units = "in")
print(Feature_plot)

# Find marker genes for each cluster and save to CSV
BLA_markers <- FindAllMarkers(BLA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(BLA_markers, file = "BLA_cluster_markers.csv", row.names = FALSE)

# Initial cluster annotations
cluster_annotations <- c(
  "0" = "Excitatory", "1" = "Oligodendrocyte", "2" = "Inhibitory", "3" = "Inhibitory",
  "4" = "Microglia", "5" = "Pericyte", "6" = "Inhibitory", "7" = "Inhibitory",
  "8" = "Astrocyte", "9" = "Inhibitory", "10" = "Excitatory", "11" = "Endothelial",
  "12" = "Mixed"
)
BLA$initial_annotation <- Idents(BLA) %>% as.character() %>% plyr::mapvalues(from = names(cluster_annotations), to = unlist(cluster_annotations))

# Refine annotations based on marker gene expression
rna_data <- GetAssayData(BLA, assay = "RNA", layer = "counts")
refined_annotations <- ifelse(
  BLA$initial_annotation == "Excitatory" & (rna_data["Slc17a7", ] > 0.5 | rna_data["Slc17a6", ] > 0.4), "Excitatory",
  ifelse(BLA$initial_annotation == "Inhibitory" & (rna_data["Gad1", ] > 1 | rna_data["Gad2", ] > 0.5 | rna_data["Slc32a1",] > 0.5), "Inhibitory",
         BLA$initial_annotation)
)
BLA$refined_annotation <- refined_annotations
BLA <- AddMetaData(BLA, metadata = refined_annotations, col.name = "cell_type")

# Save annotated plot
DimPlot(BLA, group.by = "cell_type", reduction = "wnn.umap", label = TRUE, repel = TRUE)
dim_plot2 <- DimPlot(BLA, group.by = "cell_type", reduction = "wnn.umap") +
  theme(
    legend.position = "right",
    text = element_text(size = 24),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 24),
    legend.text = element_text(size = 24),
    legend.title = element_text(size = 24),
    plot.title = element_text(size = 24)
  )
ggsave("BLA Cell-type annotation.tiff", plot = dim_plot2, width = 10, height = 10, dpi = 300, units = "in")

# Save Seurat object
saveRDS(BLA, file = file.path(base_dir, "BLA_neurons.rds"))

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
ggsave("BLA_violin_plot.tiff", plot = violin_plot, width = 12, height = 8, dpi = 300, compression = "lzw")

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
