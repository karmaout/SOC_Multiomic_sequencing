# Load the Seurat object
base_dir <- "path/to/data/"  # Update with relative or custom path
BLA_neurons <- readRDS(file.path(base_dir, "PC_neurons.rds"))

# Standard preprocessing for neurons
DefaultAssay(PC_neurons) <- "SCT"
PC_neurons <- PC_neurons %>%
  RunPCA(npcs = 100) %>%
  FindNeighbors(dims = 1:60) %>%
  FindClusters(resolution = 0.1) %>%
  RunUMAP(dims = 1:60, reduction.name = 'wnn.umap.neurons')

# Plot neuron clusters
PC_neuron_cluster <- DimPlot(PC_neurons, group.by = "seurat_clusters", reduction = "wnn.umap.neurons", pt.size = 1) +
  theme(
    legend.position = "none",
    text = element_text(size = 24),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 24),
    legend.text = element_text(size = 24),
    legend.title = element_text(size = 24),
    plot.title = element_text(size = 24)
  )
ggsave("PC_neuron_annotation_no_label.tiff", plot = PC_neuron_cluster, width = 10, height = 10, dpi = 300, units = "in")

# Find markers and plot heatmap
markers <- FindAllMarkers(PC_neurons, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top_markers <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top_marker_genes <- unique(top_markers$gene)

heatmap_plot <- DoHeatmap(PC_neurons, features = top_marker_genes, group.by = "seurat_clusters") + NoLegend()
write.csv(markers, file = "PC_neuron_cluster_markers.csv", row.names = FALSE)
ggsave("PC_Neuron_heatmap_PC.tiff", plot = heatmap_plot, width = 10, height = 10, dpi = 300)
print(top_markers)

# Subset the neurons
inhibitory_cells <- WhichCells(PC_neurons, expression = cell_type == "Inhibitory")
excitatory_cells <- WhichCells(PC_neurons, expression = cell_type == "Excitatory")
UE_cells <- WhichCells(PC_neurons, expression = behavior_group == "U_E")
PE_cells <- WhichCells(PC_neurons, expression = behavior_group == "P_E")
UR_cells <- WhichCells(PC_neurons, expression = behavior_group == "U_R")
PR_cells <- WhichCells(PC_neurons, expression = behavior_group == "P_R")

# Fos distribution of PC neurons
p4 <- FeaturePlot(PC_neurons, features = "IEG1", cells = PE_cells, reduction = "wnn.umap.neurons", 
                  pt.size = 2, cols = c("grey", "red"), alpha = 1, order = TRUE) +
  ggtitle("Paired Encoding")
ggsave("PC_PE_fos_distribution.tiff", plot = p4, width = 10, height = 10, dpi = 300, units = "in")

p5 <- FeaturePlot(PC_neurons, features = "IEG1", cells = PR_cells, reduction = "wnn.umap.neurons",
                  pt.size = 2, cols = c("grey", "red"), alpha = 1, order = TRUE) +
  ggtitle("Paired Retrieval")
ggsave("PC_PR_fos_distribution.tiff", plot = p5, width = 10, height = 10, dpi = 300, units = "in")

p6 <- FeaturePlot(PC_neurons, features = "IEG1", cells = UE_cells, reduction = "wnn.umap.neurons",
                  pt.size = 2, cols = c("grey", "red"), alpha = 1, order = TRUE) +
  ggtitle("Unpaired Encoding")
ggsave("PC_UE_fos_distribution.tiff", plot = p6, width = 10, height = 10, dpi = 300, units = "in")

p7 <- FeaturePlot(PC_neurons, features = "IEG1", cells = UR_cells, reduction = "wnn.umap.neurons",
                  pt.size = 2, cols = c("grey", "red"), alpha = 1, order = TRUE) +
  ggtitle("Unpaired Retrieval")
ggsave("PC_UR_fos_distribution.tiff", plot = p7, width = 10, height = 10, dpi = 300, units = "in")

# Annotation and plotting
neuron_annotations <- c(
  "0" = "Pyramidal 1", 
  "1" = "IN 1", # Interneurons
  "2" = "IN 2", # Semilunar
  "3" = "Vip+",
  "4" = "SL", # Semilunar
  "5" = "Sst+",
  "6" = "Pyramidal 2", # Deep pyramidal
  "7" = "VGLUT2",
  "8" = "Mixed"
)
PC_neurons$neuron_annotations <- plyr::mapvalues(PC_neurons$seurat_clusters, from = names(neuron_annotations), to = neuron_annotations)

PC_neuron_cluster_annotated <- DimPlot(PC_neurons, group.by = "neuron_annotations", reduction = "wnn.umap.neurons", label = TRUE, pt.size = 1) +
  theme(
    legend.position = "right",
    text = element_text(size = 24),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 24),
    legend.text = element_text(size = 24),
    legend.title = element_text(size = 24),
    plot.title = element_text(size = 24)
  )
print(PC_neuron_cluster_annotated)
ggsave("Neurons_annotated_PC.tiff", plot = PC_neuron_cluster_annotated, width = 10, height = 10, dpi = 300)


# List of subclusters and behavior groups to be analyzed
subclusters <- c("Pyramidal 1", "Pyramidal 2", "SL", "Vip+", "Sst+","VGLUT2","IN 1","IN 2")
behavior_groups <- c("U_E", "P_E", "U_R", "P_R")

subcluster_behavior_list <- list()
volcano_plots <- list()
de_results_list <- list()

# Iterate through each subcluster and behavior group combination
for (subcluster in subclusters) {
  for (bg in behavior_groups) {
    subset_name <- paste(subcluster, bg, sep = "_")
    subset_obj <- subset(PC_neurons, subset = neuron_annotations == subcluster & behavior_group == bg)
    
    # Check the number of cells in each cfos group
    cfos_positive_cells <- sum(subset_obj$cfos_status == "cfos+")
    cfos_negative_cells <- sum(subset_obj$cfos_status == "cfos-")
    
    if (cfos_positive_cells < 3 || cfos_negative_cells < 3) {
      print(paste("Skipping", subset_name, "because it has fewer than 3 cells in one of the cfos groups"))
      next
    }
    
    # Store the subset in the list
    subcluster_behavior_list[[subset_name]] <- subset_obj
    
    # Set Idents to cfos_status for the subset
    Idents(subset_obj) <- subset_obj$cfos_status
    
    # Perform differential expression analysis and store the results
    de_results <- FindMarkers(subset_obj, ident.1 = "cfos+", ident.2 = "cfos-", logfc.threshold = 0.5)
    
    # Exclude Fos gene and genes starting with ENSRN or AABR
    de_results <- de_results[!grepl("^(Fos|AC|ENSRN|AABR)", rownames(de_results)), ]
    
    # Add gene names as a column
    de_results$gene <- rownames(de_results)
    
    # Store the differential expression results
    de_results_list[[subset_name]] <- de_results
    
    # Sort the results by p-value to get the top 20 genes
    top_genes <- de_results[order(de_results$p_val_adj), ][1:20, "gene"]
    
    # Adjust the y-axis limits based on the behavior group

    
    # Create and store the enhanced volcano plot
    plot_title <- paste(subcluster, bg, ": cfos+ vs. cfos-", sep = " ")
    volcano_plot <- EnhancedVolcano(
      de_results,
      lab = ifelse(de_results$gene %in% top_genes, de_results$gene, NA),  # Label only top 20 genes
      x = 'avg_log2FC',  # Use log2 fold change
      y = 'p_val_adj',
      title = plot_title,  # Corrected to use plot_title
      pCutoff = 0.05,
      FCcutoff = 0.5,
      pointSize = 4.0,
      labSize = 6.0,
      col = c('grey70', 'grey70', 'green', 'green'),
      legendLabels = c('NS', 'NS', 'Downregulation', 'Upregulation'),
      legendPosition = "right",
      drawConnectors = TRUE,
      widthConnectors = 0.5,
      arrowheads = FALSE,
      colConnectors = 'grey50',
      max.overlaps = Inf,
    )
    
    # Save individual volcano plot to a high-resolution file
    file_name <- paste("volcano_plot_", subset_name, ".tiff", sep = "")
    ggsave(file_name, plot = volcano_plot, width = 8, height = 6, dpi = 300)  # Save at 300 dpi
    print(volcano_plot)
    dev.off()
  }
}

# Create a new workbook
wb <- createWorkbook()

# Add each DE result as a sheet in the workbook
for (name in names(de_results_list)) {
  addWorksheet(wb, name)
  writeData(wb, sheet = name, de_results_list[[name]])
}

# Save the workbook to a file
saveWorkbook(wb, "PC_Subcluster_Differential_Expression_Results_PC(with Fos).xlsx", overwrite = TRUE)


# Scatter plot generation
perform_de_analysis <- function(seurat_obj, ident1, ident2) {
  cells_group1 <- WhichCells(seurat_obj, idents = ident1)
  cells_group2 <- WhichCells(seurat_obj, idents = ident2)
  
  if (length(cells_group1) < 3 | length(cells_group2) < 3) {
    message(paste("Skipping DE analysis for", ident1, "vs.", ident2, "due to insufficient cells."))
    return(NULL)
  } else {
    de_results <- FindMarkers(seurat_obj, ident.1 = ident1, ident.2 = ident2, logfc.threshold = 0.25)
    de_results <- de_results[!grepl("^(Fos|ENSRN|AABR)", rownames(de_results)), ]
    de_results$gene <- rownames(de_results)
    return(de_results)
  }
}

# Function to generate scatter plots comparing DE results
generate_scatter_plot <- function(de_results_1, de_results_2, group1_name, group2_name, output_file) {
  if (is.null(de_results_1) | is.null(de_results_2)) {
    message(paste("Skipping scatter plot for", group1_name, "vs", group2_name, "due to insufficient DE results."))
    return(NULL)
  }
  
  merged_de_results <- merge(de_results_1, de_results_2, by = "gene", suffixes = c("_1", "_2"))
  filtered_de_results <- merged_de_results %>% 
    filter(!grepl("^(Fos|ENSRN|AABR)", gene))
  
  filtered_de_results$significance <- "Not Significant"
  filtered_de_results$significance[filtered_de_results$p_val_adj_1 < 0.05 & filtered_de_results$p_val_adj_2 < 0.05] <- "Significant in Both"
  filtered_de_results$significance[filtered_de_results$p_val_adj_1 < 0.05 & filtered_de_results$p_val_adj_2 >= 0.05] <- paste("Significant in", group1_name)
  filtered_de_results$significance[filtered_de_results$p_val_adj_1 >= 0.05 & filtered_de_results$p_val_adj_2 < 0.05] <- paste("Significant in", group2_name)
  
  red_genes_to_label <- filtered_de_results %>%
    filter(significance == "Significant in Both")
  
  scatter_plot <- ggplot(filtered_de_results, aes(x = avg_log2FC_1, y = avg_log2FC_2)) +
    geom_point(data = filtered_de_results %>% filter(significance == "Not Significant"),
               color = "grey70", alpha = 0.6, size = 4) +
    geom_point(data = filtered_de_results %>% filter(significance != "Not Significant"),
               aes(color = significance), alpha = 0.8, size = 8) +
    geom_text_repel(data = red_genes_to_label, 
                    aes(label = gene), 
                    size = 5, 
                    max.overlaps = Inf, 
                    box.padding = 2, 
                    point.padding = 2, 
                    segment.color = 'grey50',
                    force = 3,
                    nudge_x = 0.2,
                    nudge_y = 0.2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 1.5) +
    geom_smooth(method = "loess", se = FALSE, color = "blue", linetype = "solid", linewidth = 2) +
    labs(title = paste("Fold Change Comparison:", group1_name, "vs", group2_name), 
         x = paste(group1_name, "log2 Fold Change (cfos+ vs. cfos-)"), 
         y = paste(group2_name, "log2 Fold Change (cfos+ vs. cfos-)")) +
    scale_color_manual(values = setNames(
      c("red", "blue", "green", "grey70"),
      c("Significant in Both", paste("Significant in", group1_name), paste("Significant in", group2_name), "Not Significant")
    )) +
    theme_minimal(base_size = 15) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 20, margin = margin(t = 10, b = 10)),
      axis.title.x = element_text(size = 20, face = "bold"),
      axis.title.y = element_text(size = 20, face = "bold"),
      axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 18),
      legend.position = "none",
      panel.background = element_rect(fill = "white", color = "white"),
      plot.background = element_rect(fill = "white", color = "white")
    ) +
    coord_fixed(ratio = 1, xlim = c(-2.5, 10), ylim = c(-2.5, 10))
  
  ggsave(output_file, plot = scatter_plot, width = 10, height = 10, device = "tiff", dpi = 300, bg = "white")
}

# Perform DE analysis and generate scatter plots
for (subcluster in subclusters) {
  for (bg_pair in list(c("U_E", "P_E"), c("U_R", "P_R"))) {
    group1_name <- paste(bg_pair[1], subcluster, sep = " ")
    group2_name <- paste(bg_pair[2], subcluster, sep = " ")
    
    de_results_1 <- perform_de_analysis(subcluster_groups[[paste(bg_pair[1], subcluster, "neurons", sep = "_")]], "cfos+", "cfos-")
    de_results_2 <- perform_de_analysis(subcluster_groups[[paste(bg_pair[2], subcluster, "neurons", sep = "_")]], "cfos+", "cfos-")
    
    output_file <- paste("scatter_plot_", gsub(" ", "_", group1_name), "_vs_", gsub(" ", "_", group2_name), ".tiff", sep = "")
    generate_scatter_plot(de_results_1, de_results_2, group1_name, group2_name, output_file)
  }
}
