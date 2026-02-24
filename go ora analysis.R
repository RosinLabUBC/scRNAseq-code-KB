# Using All Upregualted and Downregulated DEGs (and DEGs have parameters of log2FC > 0.5, < -0.5, adj p-value < 0.05) 

library(clusterProfiler)
library(org.Mm.eg.db)
library(Seurat)
library(dplyr)
library(tibble)
library(enrichplot)
library(ggplot2)

# GO ORA for all CON M vs CON F --------
# Load upregulated and downregulated genes
all_CON_DEGs_top <- read.csv("/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/DEGs for Analysis/all_upregulated_CON_M_vs_F.csv")

all_CON_DEGs_down <- read.csv("/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/DEGs for Analysis/all_downregulated_CON_M_vs_F.csv")

# load in object for gene universe (use subsetted object for universe, SF said makes sense when trying to explore)

clean_combined_object <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/clean_combined_object_CON&COLD_strat2.rds')

# Set default assay to RNA
DefaultAssay(clean_combined_object) <- "RNA"

# Normalize and scale RNA assay
clean_combined_object <- NormalizeData(clean_combined_object)
clean_combined_object <- FindVariableFeatures(clean_combined_object)
clean_combined_object <- ScaleData(clean_combined_object, features = row.names(clean_combined_object))

# Join Layers
clean_combined_object <- JoinLayers(clean_combined_object, layers = "data")

# Subset to just the CON cells
con_obj <- subset(clean_combined_object, subset = treatment == "CON")

# Extract gene symbols
# checked gene symbols are indeed in the X column
up_genes <- all_CON_DEGs_top$X
down_genes <- all_CON_DEGs_down$X

# Use RNA assay with normalized expression (log-normalized)
rna_data <- GetAssayData(con_obj, assay = "RNA", layer = "data")

# Calculate percent of cells expressing each gene
pct_expr <- rowSums(rna_data > 0) / ncol(rna_data)

# Filter genes expressed in >10% of cells (min.pct = 0.1)
gene_universe <- names(pct_expr[pct_expr > 0.1])

# Combine all genes for conversion
all_genes_to_convert <- unique(c(up_genes, down_genes, gene_universe))

# Convert to ENSEMBL IDs
gene_conversion <- bitr(all_genes_to_convert,
                        fromType = "SYMBOL",
                        toType = "ENSEMBL",
                        OrgDb = org.Mm.eg.db)

# Match converted gene lists
up_ensembl <- gene_conversion %>%
  filter(SYMBOL %in% up_genes) %>%
  pull(ENSEMBL)

down_ensembl <- gene_conversion %>%
  filter(SYMBOL %in% down_genes) %>%
  pull(ENSEMBL)

universe_ensembl <- gene_conversion %>%
  filter(SYMBOL %in% gene_universe) %>%
  pull(ENSEMBL)

# GO Enrichment for Upregulated Genes
ego_up <- enrichGO(gene          = up_ensembl,
                   universe      = universe_ensembl,
                   OrgDb         = org.Mm.eg.db,
                   keyType       = "ENSEMBL",
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   qvalueCutoff  = 0.05, #SF said to remove p-value, q-value is enough
                   readable      = TRUE)

head(ego_up)

# GO Enrichment for Downregulated Genes
ego_down <- enrichGO(gene          = down_ensembl,
                     universe      = universe_ensembl,
                     OrgDb         = org.Mm.eg.db,
                     keyType       = "ENSEMBL",
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     qvalueCutoff  = 0.05, #SF said to remove p-value, q-value is enough
                     readable      = TRUE)

head(ego_down)

# Save results to your specified folder
write.csv(as.data.frame(ego_up), "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/GO ORA/CON M vs F/ORA_BP_Upregulated_CON M_vs_CON F_minpct_01.csv")
write.csv(as.data.frame(ego_down), "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/GO ORA/CON M vs F/ORA_BP_Downregulated_CON M_vs_CON F_minpct_01.csv")

# Barplot
barplot(ego_up, showCategory = 20, title = "Top GO BP Terms (CON M vs. CON F) - Upregulated")
barplot(ego_down, showCategory = 20, title = "Top GO BP Terms (CON M vs. CON F) - Downregulated")


# GO ORA for for individual clusters CON M vs CON F --------
# Set default assay
DefaultAssay(con_obj) <- "RNA"

# Create a list of the clusters you want to analyze
CON_DEGs_clusters_to_analyze <- c("Neuronal Precursors (Ascl1+)", "Glial Precursors", "Cycling NSPCs", "NPCs-Glia", "NPCs-Neuronal", "Neuronal Precursors (Neurog2+)")

# Set paths
base_path <- "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/DEGs for Analysis/"
out_path <- "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/GO ORA/CON M vs F/"

for (cluster in CON_DEGs_clusters_to_analyze) {
  cat("Processing cluster:", cluster, "\n")
  
  # Clean cluster name for filenames
  clean_cluster <- gsub(" ", "_", gsub("[\\(\\)\\+]", "", cluster))
  
  # Load DEGs
  up_file <- paste0(base_path, "upregulated_", clean_cluster, "_CON_M_vs_F.csv")
  down_file <- paste0(base_path, "downregulated_", clean_cluster, "_CON_M_vs_F.csv")
  
  up_genes <- read.csv(up_file)$X
  down_genes <- read.csv(down_file)$X
  
  # Subset Seurat object to current cluster
  Idents(con_obj) <- "cluster_names"
  seurat_subset <- subset(con_obj, idents = cluster)
  
  # Calculate gene universe
  rna_data <- GetAssayData(seurat_subset, assay = "RNA", layer = "data")
  pct_expr <- rowSums(rna_data > 0) / ncol(rna_data)
  gene_universe <- names(pct_expr[pct_expr > 0.1])
  
  # Convert gene symbols to ENSEMBL
  all_genes_to_convert <- unique(c(up_genes, down_genes, gene_universe))
  gene_conversion <- bitr(all_genes_to_convert,
                          fromType = "SYMBOL",
                          toType = "ENSEMBL",
                          OrgDb = org.Mm.eg.db)
  
  # Match converted gene lists
  up_ensembl <- gene_conversion %>% filter(SYMBOL %in% up_genes) %>% pull(ENSEMBL)
  down_ensembl <- gene_conversion %>% filter(SYMBOL %in% down_genes) %>% pull(ENSEMBL)
  universe_ensembl <- gene_conversion %>% filter(SYMBOL %in% gene_universe) %>% pull(ENSEMBL)
  
  # Define a helper function to run and save enrichment results
  run_and_save_enrichment <- function(gene_set, universe, ont, direction) {
    ego <- enrichGO(gene = gene_set,
                    universe = universe,
                    OrgDb = org.Mm.eg.db,
                    keyType = "ENSEMBL",
                    ont = ont,
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.05,
                    readable = TRUE)
    
    # Check if valid enrichResult object and has results
    if (inherits(ego, "enrichResult") && nrow(as.data.frame(ego)) > 0) {
      
      # Save CSV
      csv_file <- paste0(out_path, "GO_", ont, "_", direction, "_", clean_cluster, "_CON_M_vs_F.csv")
      write.csv(as.data.frame(ego), file = csv_file, row.names = FALSE)
      
      # Save barplot
      plot_file <- paste0(out_path, "GO_", ont, "_", direction, "_", clean_cluster, "_barplot.pdf")
      pdf(plot_file, width = 10, height = 15)
      print(barplot(ego, showCategory = 20, title = paste("Top GO", ont, "-", direction, "-", cluster)))
      dev.off()
      } else {
      message("No GO terms found for ", direction, " genes in ", cluster, " (", ont, ")")}
    
    return(ego)}
  
  # Run GO BP
  ego_up_BP <- run_and_save_enrichment(up_ensembl, universe_ensembl, "BP", "Upregulated")
  ego_down_BP <- run_and_save_enrichment(down_ensembl, universe_ensembl, "BP", "Downregulated")}

# NOTES
# No GO terms found for Upregulated genes in Neuronal Precursors (Ascl1+) (BP)
# No GO terms found for Upregulated genes in Glial Precursors (BP)
# No GO terms found for Upregulated genes in NPCs-Glia (BP)
# No GO terms found for Upregulated genes in NPCs-Neuronal (BP)
# No GO terms found for Downregulated genes in NPCs-Neuronal (BP)

# GO ORA for all CON M vs COLD M --------
# Load upregulated and downregulated genes
all_Male_DEGs_top <- read.csv("/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/DEGs for Analysis/male only/all_Male_upregulated_CON_vs_COLD.csv")

all_Male_DEGs_down <- read.csv("/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/DEGs for Analysis/male only/all_Male_downregulated_CON_vs_COLD.csv")

# load in object for gene universe (use subsetted object for universe, SF said makese sense when trying to explore)

clean_combined_object <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/clean_combined_object_CON&COLD_strat2.rds')

# Set default assay to RNA
DefaultAssay(clean_combined_object) <- "RNA"

# Normalize and scale RNA assay
clean_combined_object <- NormalizeData(clean_combined_object)
clean_combined_object <- FindVariableFeatures(clean_combined_object)
clean_combined_object <- ScaleData(clean_combined_object, features = row.names(clean_combined_object))

# Join Layers
clean_combined_object <- JoinLayers(clean_combined_object, layers = "data")

# Subset to just the Male cells
male_obj <- subset(clean_combined_object, subset = sex == "Male")

# Extract gene symbols
up_genes <- all_Male_DEGs_top$X
down_genes <- all_Male_DEGs_down$X

# Use RNA assay with normalized expression (log-normalized)
rna_data <- GetAssayData(male_obj, assay = "RNA", layer = "data")

# Calculate percent of cells expressing each gene
pct_expr <- rowSums(rna_data > 0) / ncol(rna_data)

# Filter genes expressed in >10% of cells (min.pct = 0.1)
gene_universe <- names(pct_expr[pct_expr > 0.1])

# Combine all genes for conversion
all_genes_to_convert <- unique(c(up_genes, down_genes, gene_universe))

# Convert to ENSEMBL IDs
gene_conversion <- bitr(all_genes_to_convert,
                        fromType = "SYMBOL",
                        toType = "ENSEMBL",
                        OrgDb = org.Mm.eg.db)

# Match converted gene lists
up_ensembl <- gene_conversion %>%
  filter(SYMBOL %in% up_genes) %>%
  pull(ENSEMBL)

down_ensembl <- gene_conversion %>%
  filter(SYMBOL %in% down_genes) %>%
  pull(ENSEMBL)

universe_ensembl <- gene_conversion %>%
  filter(SYMBOL %in% gene_universe) %>%
  pull(ENSEMBL)

# GO Enrichment for Upregulated Genes
ego_up <- enrichGO(gene          = up_ensembl,
                   universe      = universe_ensembl,
                   OrgDb         = org.Mm.eg.db,
                   keyType       = "ENSEMBL",
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   qvalueCutoff  = 0.05, #SF said to remove p-value, q-value is enough
                   readable      = TRUE)

head(ego_up)

# GO Enrichment for Downregulated Genes
ego_down <- enrichGO(gene          = down_ensembl,
                     universe      = universe_ensembl,
                     OrgDb         = org.Mm.eg.db,
                     keyType       = "ENSEMBL",
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     qvalueCutoff  = 0.05, #SF said to remove p-value, q-value is enough
                     readable      = TRUE)

head(ego_down)

# Save results to your specified folder
write.csv(as.data.frame(ego_up), "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/GO ORA/CON M vs COLD M/ORA_BP_Upregulated_CON M_vs_COLD M_minpct_01.csv")
write.csv(as.data.frame(ego_down), "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/GO ORA/CON M vs COLD M/ORA_BP_Downregulated_CON M_vs_COLD M_minpct_01.csv")

# Barplot
barplot(ego_up, showCategory = 20, title = "Top GO BP Terms (CON M vs. COLD M) - Upregulated")
barplot(ego_down, showCategory = 20, title = "Top GO BP Terms (CON M vs. COLD M) - Downregulated")

# NOTES
# No Downregulated BP for All CON M vs COLD M

# GO ORA for for individual clusters CON M vs COLD M --------
# Set default assay
DefaultAssay(male_obj) <- "RNA"

# Create a list of the clusters you want to analyze
Male_DEGs_clusters_to_analyze <- c("Neuronal Precursors (Ascl1+)", "Glial Precursors", "Cycling NSPCs", "NPCs-Glia", "NPCs-Neuronal", "Neuronal Precursors (Neurog2+)")

# Set paths
base_path <- "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/DEGs for Analysis/male only/"
out_path <- "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/GO ORA/CON M vs COLD M/"

for (cluster in Male_DEGs_clusters_to_analyze) {
  cat("Processing cluster:", cluster, "\n")
  
  # Clean cluster name for filenames
  clean_cluster <- gsub(" ", "_", gsub("[\\(\\)\\+]", "", cluster))
  
  # Load DEGs
  up_file <- paste0(base_path, "upregulated_", clean_cluster, "_COLD_vs_CON_M.csv")
  down_file <- paste0(base_path, "downregulated_", clean_cluster, "_COLD_vs_CON_M.csv")
  
  up_genes <- read.csv(up_file)$X
  down_genes <- read.csv(down_file)$X
  
  # Subset Seurat object to current cluster
  Idents(male_obj) <- "cluster_names"
  seurat_subset <- subset(male_obj, idents = cluster)
  
  # Calculate gene universe
  rna_data <- GetAssayData(seurat_subset, assay = "RNA", layer = "data")
  pct_expr <- rowSums(rna_data > 0) / ncol(rna_data)
  gene_universe <- names(pct_expr[pct_expr > 0.1])
  
  # Convert gene symbols to ENSEMBL
  all_genes_to_convert <- unique(c(up_genes, down_genes, gene_universe))
  gene_conversion <- bitr(all_genes_to_convert,
                          fromType = "SYMBOL",
                          toType = "ENSEMBL",
                          OrgDb = org.Mm.eg.db)
  
  # Match converted gene lists
  up_ensembl <- gene_conversion %>% filter(SYMBOL %in% up_genes) %>% pull(ENSEMBL)
  down_ensembl <- gene_conversion %>% filter(SYMBOL %in% down_genes) %>% pull(ENSEMBL)
  universe_ensembl <- gene_conversion %>% filter(SYMBOL %in% gene_universe) %>% pull(ENSEMBL)
  
  # Define a helper function to run and save enrichment results
  run_and_save_enrichment <- function(gene_set, universe, ont, direction) {
    ego <- enrichGO(gene = gene_set,
                    universe = universe,
                    OrgDb = org.Mm.eg.db,
                    keyType = "ENSEMBL",
                    ont = ont,
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.05,
                    readable = TRUE)
    
    # Check if valid enrichResult object and has results
    if (inherits(ego, "enrichResult") && nrow(as.data.frame(ego)) > 0) {
      
      # Save CSV
      csv_file <- paste0(out_path, "GO_", ont, "_", direction, "_", clean_cluster, "_CON_M_vs_COLD_M.csv")
      write.csv(as.data.frame(ego), file = csv_file, row.names = FALSE)
      
      # Save barplot
      plot_file <- paste0(out_path, "GO_", ont, "_", direction, "_", clean_cluster, "_barplot.pdf")
      pdf(plot_file, width = 10, height = 15)
      print(barplot(ego, showCategory = 20, title = paste("Top GO", ont, "-", direction, "-", cluster)))
      dev.off()
    } else {
      message("No GO terms found for ", direction, " genes in ", cluster, " (", ont, ")")}
    
    return(ego)}
  
  # Run GO BP
  ego_up_BP <- run_and_save_enrichment(up_ensembl, universe_ensembl, "BP", "Upregulated")
  ego_down_BP <- run_and_save_enrichment(down_ensembl, universe_ensembl, "BP", "Downregulated")}

# NOTES
# No GO terms found for Downregulated genes in Glial Precursors (BP)
# No GO terms found for Downregulated genes in Cycling NSPCs (BP)
# No GO terms found for Upregulated genes in Neuronal Precursors (Neurog2+) (BP)
# No GO terms found for Downregulated genes in Neuronal Precursors (Neurog2+) (BP)

# GO ORA for all CON F vs COLD F --------
all_Female_DEGs_top <- read.csv("/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/DEGs for Analysis/female only/all_Female_upregulated_CON_vs_COLD.csv")

all_Female_DEGs_down <- read.csv("/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/DEGs for Analysis/female only/all_Female_downregulated_CON_vs_COLD.csv")

# load in object for gene universe (use subsetted object for universe, SF said makese sense when trying to explore)

clean_combined_object <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/clean_combined_object_CON&COLD_strat2.rds')

# Set default assay to RNA
DefaultAssay(clean_combined_object) <- "RNA"

# Normalize and scale RNA assay
clean_combined_object <- NormalizeData(clean_combined_object)
clean_combined_object <- FindVariableFeatures(clean_combined_object)
clean_combined_object <- ScaleData(clean_combined_object, features = row.names(clean_combined_object))

# Join Layers
clean_combined_object <- JoinLayers(clean_combined_object, layers = "data")

# Subset to just the Female cells
female_obj <- subset(clean_combined_object, subset = sex == "Female")

# Extract gene symbols
up_genes <- all_Female_DEGs_top$X
down_genes <- all_Female_DEGs_down$X

# Use RNA assay with normalized expression (log-normalized)
rna_data <- GetAssayData(female_obj, assay = "RNA", layer = "data")

# Calculate percent of cells expressing each gene
pct_expr <- rowSums(rna_data > 0) / ncol(rna_data)

# Filter genes expressed in >10% of cells (min.pct = 0.1)
gene_universe <- names(pct_expr[pct_expr > 0.1])

# Combine all genes for conversion
all_genes_to_convert <- unique(c(up_genes, down_genes, gene_universe))

# Convert to ENSEMBL IDs
gene_conversion <- bitr(all_genes_to_convert,
                        fromType = "SYMBOL",
                        toType = "ENSEMBL",
                        OrgDb = org.Mm.eg.db)

# Match converted gene lists
up_ensembl <- gene_conversion %>%
  filter(SYMBOL %in% up_genes) %>%
  pull(ENSEMBL)

down_ensembl <- gene_conversion %>%
  filter(SYMBOL %in% down_genes) %>%
  pull(ENSEMBL)

universe_ensembl <- gene_conversion %>%
  filter(SYMBOL %in% gene_universe) %>%
  pull(ENSEMBL)

# GO Enrichment for Upregulated Genes
ego_up <- enrichGO(gene          = up_ensembl,
                   universe      = universe_ensembl,
                   OrgDb         = org.Mm.eg.db,
                   keyType       = "ENSEMBL",
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   qvalueCutoff  = 0.05, #SF said to remove p-value, q-value is enough
                   readable      = TRUE)

head(ego_up)

# GO Enrichment for Downregulated Genes
ego_down <- enrichGO(gene          = down_ensembl,
                     universe      = universe_ensembl,
                     OrgDb         = org.Mm.eg.db,
                     keyType       = "ENSEMBL",
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     qvalueCutoff  = 0.05, #SF said to remove p-value, q-value is enough
                     readable      = TRUE)

head(ego_down)

# Save results to your specified folder
write.csv(as.data.frame(ego_up), "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/GO ORA/CON F vs COLD F/ORA_BP_Upregulated_CON F_vs_COLD F_minpct_01.csv")
write.csv(as.data.frame(ego_down), "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/GO ORA/CON F vs COLD F/ORA_BP_Downregulated_CON F_vs_COLD F_minpct_01.csv")

# Barplot
barplot(ego_up, showCategory = 20, title = "Top GO BP Terms (CON F vs. COLD F) - Upregulated")
barplot(ego_down, showCategory = 20, title = "Top GO BP Terms (CON F vs. COLD F) - Downregulated")


# GO ORA for for individual clusters CON F vs COLD F --------
# Set default assay
DefaultAssay(female_obj) <- "RNA"

# Create a list of the clusters you want to analyze
Female_DEGs_clusters_to_analyze <- c("Neuronal Precursors (Ascl1+)", "Glial Precursors", "Cycling NSPCs", "NPCs-Glia", "NPCs-Neuronal", "Neuronal Precursors (Neurog2+)")

# Set paths
base_path <- "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/DEGs for Analysis/female only/"
out_path <- "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/GO ORA/CON F vs COLD F/"

for (cluster in Female_DEGs_clusters_to_analyze) {
  cat("Processing cluster:", cluster, "\n")
  
  # Clean cluster name for filenames
  clean_cluster <- gsub(" ", "_", gsub("[\\(\\)\\+]", "", cluster))
  
  # Load DEGs
  up_file <- paste0(base_path, "upregulated_", clean_cluster, "_COLD_vs_CON_F.csv")
  down_file <- paste0(base_path, "downregulated_", clean_cluster, "_COLD_vs_CON_F.csv")
  
  up_genes <- read.csv(up_file)$X
  down_genes <- read.csv(down_file)$X
  
  # Subset Seurat object to current cluster
  Idents(female_obj) <- "cluster_names"
  seurat_subset <- subset(female_obj, idents = cluster)
  
  # Calculate gene universe
  rna_data <- GetAssayData(seurat_subset, assay = "RNA", layer = "data")
  pct_expr <- rowSums(rna_data > 0) / ncol(rna_data)
  gene_universe <- names(pct_expr[pct_expr > 0.1])
  
  # Convert gene symbols to ENSEMBL
  all_genes_to_convert <- unique(c(up_genes, down_genes, gene_universe))
  gene_conversion <- bitr(all_genes_to_convert,
                          fromType = "SYMBOL",
                          toType = "ENSEMBL",
                          OrgDb = org.Mm.eg.db)
  
  # Match converted gene lists
  up_ensembl <- gene_conversion %>% filter(SYMBOL %in% up_genes) %>% pull(ENSEMBL)
  down_ensembl <- gene_conversion %>% filter(SYMBOL %in% down_genes) %>% pull(ENSEMBL)
  universe_ensembl <- gene_conversion %>% filter(SYMBOL %in% gene_universe) %>% pull(ENSEMBL)
  
  # Define a helper function to run and save enrichment results
  run_and_save_enrichment <- function(gene_set, universe, ont, direction) {
    ego <- enrichGO(gene = gene_set,
                    universe = universe,
                    OrgDb = org.Mm.eg.db,
                    keyType = "ENSEMBL",
                    ont = ont,
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.05,
                    readable = TRUE)
    
    # Check if valid enrichResult object and has results
    if (inherits(ego, "enrichResult") && nrow(as.data.frame(ego)) > 0) {
      
      # Save CSV
      csv_file <- paste0(out_path, "GO_", ont, "_", direction, "_", clean_cluster, "_CON_F_vs_COLD_F.csv")
      write.csv(as.data.frame(ego), file = csv_file, row.names = FALSE)
      
      # Save barplot
      plot_file <- paste0(out_path, "GO_", ont, "_", direction, "_", clean_cluster, "_barplot.pdf")
      pdf(plot_file, width = 10, height = 15)
      print(barplot(ego, showCategory = 20, title = paste("Top GO", ont, "-", direction, "-", cluster)))
      dev.off()
    } else {
      message("No GO terms found for ", direction, " genes in ", cluster, " (", ont, ")")}
    
    return(ego)}
  
  # Run GO BP
  ego_up_BP <- run_and_save_enrichment(up_ensembl, universe_ensembl, "BP", "Upregulated")
  ego_down_BP <- run_and_save_enrichment(down_ensembl, universe_ensembl, "BP", "Downregulated")}

# NOTES
# No GO terms found for Downregulated genes in Neuronal Precursors (Ascl1+) (BP)
# No GO terms found for Downregulated genes in NPCs-Neuronal (BP)
# No GO terms found for Upregulated genes in Neuronal Precursors (Neurog2+) (BP)
# No GO terms found for Downregulated genes in Neuronal Precursors (Neurog2+) (BP)
