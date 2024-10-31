# Load the Seurat object
base_dir <- "path/to/data/"  # Update with relative or custom path
BLA_neurons <- readRDS(file.path(base_dir, "BLA_neurons.rds"))


# Set the default assay to SCT and perform standard preprocessing
DefaultAssay(BLA_neurons) <- "SCT"
BLA_neurons <- BLA_neurons %>%
  RunPCA(npcs = 100) %>%
  FindNeighbors(dims = 1:60) %>%
  FindClusters(resolution = 0.1) %>%
  RunUMAP(dims = 1:60, reduction.name = 'wnn.umap.neurons')

# Plot and save UMAP of neurons by cluster
neuron_cluster <- DimPlot(BLA_neurons, group.by = "seurat_clusters", reduction = "wnn.umap.neurons", label = TRUE, pt.size = 1) +
  theme(
    legend.position = "right",
    text = element_text(size = 24),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 24),
    legend.text = element_text(size = 24),
    legend.title = element_text(size = 24),
    plot.title = element_text(size = 24)
  )
ggsave("BLA_neuron_annotation.tiff", plot = neuron_cluster, width = 10, height = 10, dpi = 300, units = "in")

# Find marker genes for each cluster and save results
markers <- FindAllMarkers(BLA_neurons, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, file = "BLA_neuron_cluster_markers.csv", row.names = FALSE)

# View top marker genes for each cluster
top_markers <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

# Define and plot heatmap for top marker genes
top_marker_genes <- unique(top_markers$gene)
heatmap_plot <- DoHeatmap(BLA_neurons, features = top_marker_genes, group.by = "seurat_clusters") + NoLegend()
ggsave("BLA_Neuron_heatmap.tiff", plot = heatmap_plot, width = 10, height = 10, dpi = 300)

# Plot Fos distribution across BLA cells
p9 <- FeaturePlot(BLA_neurons, features = "IEG1", reduction = "wnn.umap.neurons", pt.size = 2, cols = c("grey", "red"), alpha = 1, order = TRUE) +
  ggtitle("Fos distribution across BLA cells")
ggsave("Fos_distribution_across_BLA_cells.tiff", plot = p9, width = 10, height = 10, dpi = 300, units = "in")

# Subset neurons by cell type and behavior group
inhibitory_cells <- WhichCells(BLA_neurons, expression = cell_type == "Inhibitory")
excitatory_cells <- WhichCells(BLA_neurons, expression = cell_type == "Excitatory")
UE_cells <- WhichCells(BLA_neurons, expression = behavior_group == "U_E")
PE_cells <- WhichCells(BLA_neurons, expression = behavior_group == "P_E")
UR_cells <- WhichCells(BLA_neurons, expression = behavior_group == "U_R")
PR_cells <- WhichCells(BLA_neurons, expression = behavior_group == "P_R")

# Plot Fos distribution in different behavior groups and save plots
plot_behavior_group <- function(feature, cells, group_label, filename) {
  p <- FeaturePlot(BLA_neurons, features = feature, cells = cells, reduction = "wnn.umap.neurons",
                   pt.size = 2, cols = c("grey", "red"), alpha = 1, order = TRUE) +
    ggtitle(group_label)
  ggsave(filename, plot = p, width = 10, height = 10, dpi = 300, units = "in")
}

plot_behavior_group("IEG1", PE_cells, "Paired Encoding", "PE_Fos_distribution.tiff")
plot_behavior_group("IEG1", PR_cells, "Paired Retrieval", "PR_Fos_distribution.tiff")
plot_behavior_group("IEG1", UE_cells, "Unpaired Encoding", "UE_Fos_distribution.tiff")
plot_behavior_group("IEG1", UR_cells, "Unpaired Retrieval", "UR_Fos_distribution.tiff")

# Annotate clusters with defined neuron types
BLA_neuron_annotations <- c(
  "0" = "VGLUT1",
  "1" = "Sst+",
  "2" = "MSNs",
  "3" = "Npy+",
  "4" = "PV+ with Sst enriched",
  "5" = "Vip+",
  "6" = "Unknown",
  "7" = "VGLUT2",
  "8" = "Unknown"
)

BLA_neurons$neuron_annotations <- plyr::mapvalues(BLA_neurons$seurat_clusters, from = names(BLA_neuron_annotations), to = BLA_neuron_annotations)

# Verify annotations and plot UMAP with new annotations
table(BLA_neurons$neuron_annotations)
BLA_neuron_cluster_annotated <- DimPlot(BLA_neurons, group.by = "neuron_annotations", reduction = "wnn.umap.neurons", label = TRUE, pt.size = 1) +
  theme(
    legend.position = "right",
    text = element_text(size = 24),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 24),
    legend.text = element_text(size = 24),
    legend.title = element_text(size = 24),
    plot.title = element_text(size = 24)
  )
ggsave("BLA_neuron_Cell_type.tiff", plot = BLA_neuron_cluster_annotated, width = 10, height = 10, dpi = 300, units = "in")



# List of subclusters and behavior groups
subclusters <- c("VGLUT1", "Sst+", "MSNs", "Npy+", "PV+ with Sst enriched", "Vip+", "VGLUT2")
behavior_groups <- c("U_E", "P_E", "U_R", "P_R")

# Initialize lists to store subsets, volcano plots, and differential expression results
subcluster_behavior_list <- list()
volcano_plots <- list()
de_results_list <- list()

# Function to perform differential expression analysis and return cleaned DE results
perform_de_analysis <- function(seurat_obj, ident1, ident2, subset_name) {
  cfos_positive_cells <- sum(seurat_obj$cfos_status == ident1)
  cfos_negative_cells <- sum(seurat_obj$cfos_status == ident2)
  
  if (cfos_positive_cells < 3 || cfos_negative_cells < 3) {
    message(paste("Skipping", subset_name, "due to insufficient cells in one of the cfos groups."))
    return(NULL)
  }
  
  Idents(seurat_obj) <- seurat_obj$cfos_status
  de_results <- FindMarkers(seurat_obj, ident.1 = ident1, ident.2 = ident2, logfc.threshold = 0.5)
  de_results <- de_results[!grepl("^(Fos|AC|ENSRN|AABR)", rownames(de_results)), ]
  de_results$gene <- rownames(de_results)
  return(de_results)
}

# Function to generate volcano plot for given DE results
generate_volcano_plot <- function(de_results, plot_title, y_limit, output_file) {
  top_genes <- de_results[order(de_results$p_val_adj), ][1:20, "gene"]
  
  volcano_plot <- EnhancedVolcano(
    de_results,
    lab = ifelse(de_results$gene %in% top_genes, de_results$gene, NA),
    x = 'avg_log2FC',
    y = 'p_val_adj',
    title = plot_title,
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
    xlim = c(-5, 10),
    ylim = c(0, y_limit)
  )
  
  ggsave(output_file, plot = volcano_plot, width = 8, height = 6, dpi = 300)
}

# Perform DE analysis and generate volcano plots for each subcluster and behavior group combination
for (subcluster in subclusters) {
  for (bg in behavior_groups) {
    subset_name <- paste(subcluster, bg, sep = "_")
    subset_obj <- subset(BLA_neurons, subset = neuron_annotations == subcluster & behavior_group == bg)
    
    de_results <- perform_de_analysis(subset_obj, "cfos+", "cfos-", subset_name)
    if (is.null(de_results)) next
    
    # Store DE results
    de_results_list[[subset_name]] <- de_results
    
    # Generate volcano plot
    y_limit <- if (bg == "U_E") 20 else 15
    plot_title <- paste(subcluster, bg, ": cfos+ vs. cfos-", sep = " ")
    output_file <- paste("volcano_plot_", subset_name, ".tiff", sep = "")
    generate_volcano_plot(de_results, plot_title, y_limit, output_file)
  }
}

# Create a new workbook and save DE results to an Excel file
wb <- createWorkbook()
for (name in names(de_results_list)) {
  addWorksheet(wb, name)
  writeData(wb, sheet = name, de_results_list[[name]])
}
saveWorkbook(wb, "BLA_Subcluster_Differential_Expression_Results.xlsx", overwrite = TRUE)

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
