# Load the Seurat object
base_dir <- "path/to/data/"  # Update with relative or custom path
REGION_neurons <- readRDS(file.path(base_dir, "REGION_neurons.rds"))

# List of behavior groups and library IDs
behavior_groups <- c("P_E", "P_R")
library_ids <- c("PE1", "PE2", "PE3", "PR1", "PR2", "PR3")
cfos_statuses <- c("cfos+", "cfos-")

# Subset neurons for each behavior group, cfos status, and library
subset_neurons <- list()
for (bg in behavior_groups) {
  for (status in cfos_statuses) {
    group_name <- paste("REGION", bg, "VGLUT1", status, "neurons", sep = "_")
    subset_neurons[[group_name]] <- subset(BLA_neurons, subset = neuron_annotations == "VGLUT1" & behavior_group == bg & cfos_status == status)
    
    # Subset further by library ID
    for (lib in library_ids) {
      if (grepl(bg, lib)) {  # Ensure library matches behavior group
        lib_group_name <- paste(group_name, lib, sep = "_")
        subset_neurons[[lib_group_name]] <- subset(subset_neurons[[group_name]], subset = library_id == lib)
      }
    }
  }
}

# Extract ATAC ChromatinAssay objects for specific groups
REGION_PE_VGLUT1_cfos_positive_atac <- subset_neurons[["REGION_P_E_VGLUT1_cfos+_neurons"]]@assays$ATAC
REGION_PE_VGLUT1_cfos_negative_atac <- subset_neurons[["REGION_P_E_VGLUT1_cfos-_neurons"]]@assays$ATAC
REGION_PR_VGLUT1_cfos_positive_atac <- subset_neurons[["REGION_P_R_VGLUT1_cfos+_neurons"]]@assays$ATAC
REGION_PR_VGLUT1_cfos_negative_atac <- subset_neurons[["REGION_P_R_VGLUT1_cfos-_neurons"]]@assays$ATAC

# Function to create BED file from a ChromatinAssay object
create_bed_file <- function(chromatin_assay, file_name) {
  # Access the peak ranges
  peaks <- granges(chromatin_assay)
  
  # Convert the counts matrix to a data frame
  counts_df <- as.data.frame(chromatin_assay@counts)
  
  # Get the chromosome, start, and end information
  chromosome <- as.character(seqnames(peaks))
  
  # Add "chr" prefix to each chromosome name
  chromosome <- paste0("chr", chromosome)
  
  start <- start(peaks) - 1  
  end <- end(peaks)
  
  # Compute total reads across all cells for this condition
  total_reads <- sum(chromatin_assay@counts)
  
  # Calculate normalized read counts as a score
  score <- (rowSums(counts_df) / total_reads) * 1e6
  
  # Prepare the BED-like data frame
  bed_df <- data.frame(
    chromosome = chromosome,
    start = start,
    end = end,
    name = ".",  
    score = score,
    strand = ".",  
    counts = rowSums(counts_df)  # Adding counts for reference
  )
  
  # Write to a BED file
  write.table(
    bed_df,
    file = file_name,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}

# Generate BED files for each ChromatinAssay object
create_bed_file(REGION_PE_VGLUT1_cfos_positive_atac, "REGION_PE_VGLUT1_cfos_positive_atac.bed")
create_bed_file(REGION_PE_VGLUT1_cfos_negative_atac, "REGION_PE_VGLUT1_cfos_negative_atac.bed")
create_bed_file(REGION_PR_VGLUT1_cfos_positive_atac, "REGION_PR_VGLUT1_cfos_positive_atac.bed")
create_bed_file(REGION_PR_VGLUT1_cfos_negative_atac, "REGION_PR_VGLUT1_cfos_negative_atac.bed")
