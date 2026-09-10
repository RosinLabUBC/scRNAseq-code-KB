library(clusterProfiler)
library(org.Mm.eg.db)
library(Seurat)
library(dplyr)
library(tibble)
library(enrichplot)
library(ggplot2)

# 1. Load Upregulated and Downregulated DEGs for Cycling NSPCs
cycling_CON_DEGs_top <- read.csv("/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/code for suppl 4/blackshaw + kb (analysis 1)/analysis on Vim+ subset/DEGs/Cycling NSPCs_upregulated_CON_M_vs_F.csv")

cycling_CON_DEGs_down <- read.csv("/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/code for suppl 4/blackshaw + kb (analysis 1)/analysis on Vim+ subset/DEGs/Cycling NSPCs_downregulated_CON_M_vs_F.csv")

# Extract gene symbols (assuming the gene column is named 'X' or 'gene' — change if needed)
up_genes <- cycling_CON_DEGs_top$X
down_genes <- cycling_CON_DEGs_down$X

# 2. Load Seurat Object
clean_combined_object <- readRDS('/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/code for suppl 4/blackshaw + kb (analysis 1)/clean_combined_object_blackshaw_and_kb_classified_cycling_nspcs.rds')

DefaultAssay(clean_combined_object) <- "RNA"
clean_combined_object <- NormalizeData(clean_combined_object)
clean_combined_object <- FindVariableFeatures(clean_combined_object)
clean_combined_object <- ScaleData(clean_combined_object, features = row.names(clean_combined_object))
clean_combined_object <- JoinLayers(clean_combined_object, layers = "data")

# 3. Subset to Control cells ONLY (matching your multi-cluster script logic)
con_obj <- subset(clean_combined_object, subset = treatment == "CON")

# Target only Cycling NSPCs
cluster <- "Cycling NSPCs"
cat("Processing cluster:", cluster, "\n")

# Clean cluster name for filenames
clean_cluster <- gsub(" ", "_", gsub("[\\(\\)\\+]", "", cluster))

# Set output path (update if you want these saved somewhere specific)
out_path <- "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/code for suppl 4/suppl 4 images/GO ORA/"
dir.create(out_path, showWarnings = FALSE, recursive = TRUE)

# 3. Subset Seurat object to Cycling NSPCs for gene universe calculation
# (Make sure "cluster_names" is your correct metadata column name containing cluster identities)
Idents(con_obj) <- "cluster_names"
seurat_subset <- subset(con_obj, idents = cluster)

# Calculate gene universe (>10% expression in the subset)
rna_data <- GetAssayData(seurat_subset, assay = "RNA", layer = "data")
pct_expr <- rowSums(rna_data > 0) / ncol(rna_data)
gene_universe <- names(pct_expr[pct_expr > 0.1])

# 4. Convert gene symbols to ENSEMBL
all_genes_to_convert <- unique(c(up_genes, down_genes, gene_universe))
gene_conversion <- bitr(all_genes_to_convert,
                        fromType = "SYMBOL",
                        toType = "ENSEMBL",
                        OrgDb = org.Mm.eg.db)

# Match converted gene lists
up_ensembl <- gene_conversion %>% filter(SYMBOL %in% up_genes) %>% pull(ENSEMBL)
down_ensembl <- gene_conversion %>% filter(SYMBOL %in% down_genes) %>% pull(ENSEMBL)
universe_ensembl <- gene_conversion %>% filter(SYMBOL %in% gene_universe) %>% pull(ENSEMBL)

# 5. Define helper function to run and save enrichment results
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
    
    message("Successfully saved GO results for ", direction, " genes.")
  } else {
    message("No GO terms found for ", direction, " genes in ", cluster, " (", ont, ")")
}
  
  return(ego)}

# 6. Run GO BP for Cycling NSPCs
ego_up_BP <- run_and_save_enrichment(up_ensembl, universe_ensembl, "BP", "Upregulated")
ego_down_BP <- run_and_save_enrichment(down_ensembl, universe_ensembl, "BP", "Downregulated")