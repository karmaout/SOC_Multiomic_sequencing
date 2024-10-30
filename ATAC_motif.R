# Load required libraries
library(Seurat)
library(dplyr)

# Load the Seurat object
REGION_neurons <- readRDS("REGION_neurons.rds") # Load without specific directory path

# Set default assay for SCT data extraction
DefaultAssay(REGION_neurons) <- "SCT"

# Subset the Seurat object to only include VGLUT1 neurons
PE_VGLUT1_neurons <- subset(REGION_neurons, subset = neuron_annotations == "VGLUT1" & behavior_group == "P_E")

# Find the markers between cfos+ and cfos- cells in the VGLUT1 subset
deg_results <- FindMarkers(
  object = PE_VGLUT1_neurons,
  ident.1 = "cfos+",
  ident.2 = "cfos-",
  group.by = "cfos_status"
)

# Extract relevant columns: Gene names and log fold changes
deg_results$Gene.Name <- rownames(deg_results)

# Read ATAC files
REGION_PE_VGLUT1_cfos_positive <- read.table('REGION_PE_VGLUT1_cfos_positive_atac.bed')
REGION_PE_VGLUT1_cfos_negative <- read.table('REGION_PE_VGLUT1_cfos_negative_atac.bed')

# Combine ATAC scores for samples
PE_score <- as.data.frame(cbind(REGION_PE_VGLUT1_cfos_positive$V1, REGION_PE_VGLUT1_cfos_positive$V2, REGION_PE_VGLUT1_cfos_positive$V3, REGION_PE_VGLUT1_cfos_positive$V5, REGION_PE_VGLUT1_cfos_negative$V5))

# Label data frame's header and rows
colnames(PE_score) <- c('chr', 'start', 'end', 'fos_pos_atac', 'fos_neg_atac')
rownames(PE_score) <- as.character(paste('rn7', PE_score$chr, as.character(as.numeric(PE_score$start) + 1), PE_score$end, '+', sep = '_'))

# Read motif file (RSAT output)
fos.rsat <- read.table("matrix-scan_results.txt", sep = '\t', header = TRUE)

# Select best motif per peak
fos.motif <- fos.rsat %>%
  group_by(seq_id) %>%
  filter(weight == max(weight)) %>%
  slice(1) %>%  # This selects only the first row when there are multiple rows with the same max weight
  select(seq_id, ft_name, strand, start, end, sequence, weight)

# Format as data frame and label rows
fos.motif <- as.data.frame(fos.motif)
rownames(fos.motif) <- as.character(fos.motif$seq_id)

# Add motif information with ATAC peak scores
PE_score$motif_weight <- fos.motif[rownames(PE_score), 'weight']
PE_score$sequence <- fos.motif[rownames(PE_score), 'sequence']
PE_score$strand <- fos.motif[rownames(PE_score), 'strand']

# Read annotation file (HOMER output)
annotation <- read.delim("coord-anno.txt", sep = '\t', header = TRUE)
rownames(annotation) <- as.character(paste('rn7', annotation$Chr, annotation$Start, annotation$End, '+', sep = '_'))

# Integrate gene annotation to motif information and ATAC peak scores
PE_score$Gene.Name <- annotation[rownames(PE_score), 'Gene.Name']
PE_score$Annotation <- annotation[rownames(PE_score), 'Annotation']
PE_score$Distance.to.TSS <- annotation[rownames(PE_score), 'Distance.to.TSS']

# Integrate fold change information from `FindMarkers` results
PE_score <- PE_score %>%
  left_join(deg_results %>% select(Gene.Name, avg_log2FC), by = "Gene.Name")

# Write file for plots in Python
write.table(x = PE_score, file = 'REGION_PE_anno-sorted4_foldchange.txt', quote = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE)

# Convert ATAC scores to numeric
PE_score$fos_pos_atac <- as.numeric(PE_score$fos_pos_atac)
PE_score$fos_neg_atac <- as.numeric(PE_score$fos_neg_atac)

# Select only rows where "promoter" is present in the "Annotation" column
promoter_data <- PE_score %>%
  filter(grepl("promoter", Annotation))

# Group by Gene.Name and keep the row with the maximum log-transformed ATAC accessibility ratio
result2 <- promoter_data %>%
  group_by(Gene.Name) %>%
  slice(which.max(log((fos_pos_atac + 1) / (fos_neg_atac + 1))))

# Plot the smooth scatter using SCT log fold change (avg_log2FC) and the maximum log-transformed ATAC ratio
smoothScatter(result2$avg_log2FC, log((result2$fos_pos_atac + 1) / (result2$fos_neg_atac + 1)), nrpoints = 100, cex = 0.1,
              main = "REGION PE RNA ATAC fold change density", xlab = "Log2 Fold Change (avg_log2FC)", ylab = "Max Log ATAC Ratio")


