# want DEGs list and volcano plots for CON M vs CON F, CON M vs COLD M, CON F vs COLD F
# for other clusters besides cycling NSPCs (drawing file DEGs from "DEGs for ORA workflow")

library(Seurat)
library(dplyr)
library(patchwork)
library(clustree)
library(pheatmap)
library(EnhancedVolcano)
library(ggplot2)
library(ggrepel)

# CON M vs. CON F (other clusters besides Cycling NSPCs) ------------------------------------------

# Define cluster names
clusters <- c("Neuronal_Precursors_Ascl1", 
              "Glial_Precursors",
              "NPCs-Glia", 
              "NPCs-Neuronal", 
              "Neuronal_Precursors_Neurog2")

# Define base path for DEG files
basepath <- "/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/DEG objects/"

# Loop over clusters
for (cluster in clusters) {
  
  # Construct file path
  filepath <- paste0(basepath, cluster, "_CON_DEGs_Male_vs_Female_(no_p-value_cutoff).rds")
  
  # Read in DEGs file
  DEGs <- readRDS(filepath)
  
  # Add gene names
  DEGs$gene <- rownames(DEGs)
  
  # Categorize differential expression
  DEGs$diffexpressed <- "NO"
  DEGs$diffexpressed[DEGs$avg_log2FC > 0 & DEGs$p_val_adj < 0.05] <- "UP"
  DEGs$diffexpressed[DEGs$avg_log2FC < 0 & DEGs$p_val_adj < 0.05] <- "DOWN"
  
  # Create volcano plot
  p <- ggplot(DEGs, aes(x = avg_log2FC, y = -log10(p_val_adj), color = diffexpressed)) +
    geom_point(size = 2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
    scale_color_manual(values = c("DOWN" = "#00AFBB", "NO" = "gray", "UP" = "darkred"),
                       labels = c("DOWN" = "Downregulated", "NO" = "Not significant", "UP" = "Upregulated")) +
    theme_minimal() +
    labs(title = paste("Volcano Plot -", cluster, "(CON M vs CON F)"),
         x = "Log2 Fold Change",
         y = "-log10(Adjusted p-value)",
         color = "Expression") +
    geom_text_repel(data = subset(DEGs, diffexpressed != "NO"),
                    aes(label = gene), size = 3, max.overlaps = 20)
  
  # Print plot to viewer
  print(p)
  
  # Optional message
  message("Displayed volcano plot for: ", cluster)}


# CON M vs. COLD M (other clusters besides Cycling NSPCs) ------------------------------------------

# Define cluster names
clusters <- c("Neuronal_Precursors_Ascl1", 
              "Glial_Precursors",
              "NPCs-Glia", 
              "NPCs-Neuronal", 
              "Neuronal_Precursors_Neurog2")

# Define base path for DEG files
basepath <- "/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/DEG objects/"

# Loop over clusters
for (cluster in clusters) {
  
  # Construct file path
  filepath <- paste0(basepath, cluster, "_Male_DEGs_CON_vs_COLD_(no_p-value_cutoff).rds")
  
  # Read in DEGs file
  DEGs <- readRDS(filepath)
  
  # Add gene names
  DEGs$gene <- rownames(DEGs)
  
  # Categorize differential expression
  DEGs$diffexpressed <- "NO"
  DEGs$diffexpressed[DEGs$avg_log2FC > 0 & DEGs$p_val_adj < 0.05] <- "UP"
  DEGs$diffexpressed[DEGs$avg_log2FC < 0 & DEGs$p_val_adj < 0.05] <- "DOWN"
  
  # Create volcano plot
  p <- ggplot(DEGs, aes(x = avg_log2FC, y = -log10(p_val_adj), color = diffexpressed)) +
    geom_point(size = 2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
    scale_color_manual(values = c("DOWN" = "#00AFBB", "NO" = "gray", "UP" = "darkred"),
                       labels = c("DOWN" = "Downregulated", "NO" = "Not significant", "UP" = "Upregulated")) +
    theme_minimal() +
    labs(title = paste("Volcano Plot -", cluster, "(CON M vs COLD M)"),
         x = "Log2 Fold Change",
         y = "-log10(Adjusted p-value)",
         color = "Expression") +
    geom_text_repel(data = subset(DEGs, diffexpressed != "NO"),
                    aes(label = gene), size = 3, max.overlaps = 20)
  
  # Print plot to viewer
  print(p)
  
  # Optional message
  message("Displayed volcano plot for: ", cluster)}


# CON F vs. COLD F (other clusters besides Cycling NSPCs) ------------------------------------------

# Define cluster names
clusters <- c("Neuronal_Precursors_Ascl1", 
              "Glial_Precursors",
              "NPCs-Glia", 
              "NPCs-Neuronal", 
              "Neuronal_Precursors_Neurog2")

# Define base path for DEG files
basepath <- "/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/DEG objects/"

# Loop over clusters
for (cluster in clusters) {
  
  # Construct file path
  filepath <- paste0(basepath, cluster, "_Female_DEGs_CON_vs_COLD_(no_p-value_cutoff).rds")
  
  # Read in DEGs file
  DEGs <- readRDS(filepath)
  
  # Add gene names
  DEGs$gene <- rownames(DEGs)
  
  # Categorize differential expression
  DEGs$diffexpressed <- "NO"
  DEGs$diffexpressed[DEGs$avg_log2FC > 0 & DEGs$p_val_adj < 0.05] <- "UP"
  DEGs$diffexpressed[DEGs$avg_log2FC < 0 & DEGs$p_val_adj < 0.05] <- "DOWN"
  
  # Create volcano plot
  p <- ggplot(DEGs, aes(x = avg_log2FC, y = -log10(p_val_adj), color = diffexpressed)) +
    geom_point(size = 2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
    scale_color_manual(values = c("DOWN" = "#00AFBB", "NO" = "gray", "UP" = "darkred"),
                       labels = c("DOWN" = "Downregulated", "NO" = "Not significant", "UP" = "Upregulated")) +
    theme_minimal() +
    labs(title = paste("Volcano Plot -", cluster, "(CON F vs COLD F)"),
         x = "Log2 Fold Change",
         y = "-log10(Adjusted p-value)",
         color = "Expression") +
    geom_text_repel(data = subset(DEGs, diffexpressed != "NO"),
                    aes(label = gene), size = 3, max.overlaps = 20)
  
  # Print plot to viewer
  print(p)
  
  # Optional message
  message("Displayed volcano plot for: ", cluster)}
