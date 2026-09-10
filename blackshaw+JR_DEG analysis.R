
library(Seurat)
library(dplyr)
library(patchwork)
library(clustree)
library(pheatmap)
library(EnhancedVolcano)

# load in clean combined object -------------------------------------------

clean_combined_object <- readRDS('/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/blackshaw + kb/clean_combined_object_blackshaw_and_kb_classified_cycling_nspcs.rds')

# Set default assay to RNA
DefaultAssay(clean_combined_object) <- "RNA"

# Normalize and scale RNA assay
clean_combined_object <- NormalizeData(clean_combined_object)
clean_combined_object <- FindVariableFeatures(clean_combined_object)
clean_combined_object <- ScaleData(clean_combined_object, features = row.names(clean_combined_object))

# Join Layers
clean_combined_object <- JoinLayers(clean_combined_object, layers = "data")

# DEGs on Cycling NSPCs in CON M vs CON F ----------------------------------------------
# Subset to just the CON cells
con_obj <- subset(clean_combined_object, subset = treatment == "CON")

# Subset to just the cycling NSPCs
cycling_con_obj <- subset(con_obj, subset = cluster_names == "Cycling NSPCs")

# Confirm identities
Idents(cycling_con_obj) <- "sex"

# Run DE analysis for Male vs. Female within this cluster
cycling_DEGs_CON_M_vs_CON_F <- FindMarkers(
  cycling_con_obj,
  ident.1 = "Male",
  ident.2 = "Female",
  assay = "RNA",
  logfc.threshold = 0.25,
  min.pct = 0.1)

# Save as CSV
outdir <- "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/blackshaw + kb/analysis on Vim+ subset/DEGs/"

write.csv(
  cycling_DEGs_CON_M_vs_CON_F,
  file = paste0(outdir, "CON_DEGs_Male_vs_Female_CyclingNPCs_no_pcutoff.csv"))
  
  saveRDS(
    cycling_DEGs_CON_M_vs_CON_F,
    file = paste0(outdir, "CON_DEGs_Male_vs_Female_CyclingNSPCs_no_pcutoff.rds"))

# Save all upregulated DEGs (positive log2FC)
all_CON_DEGs_up <- cycling_DEGs_CON_M_vs_CON_F %>%
  filter(avg_log2FC > 0.5) %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC))
write.csv(all_CON_DEGs_up, "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/blackshaw + kb/analysis on Vim+ subset/DEGs/all_upregulated_CON_M_vs_F.csv")

# Save all downregulated DEGs (negative log2FC)
all_CON_DEGs_down <- cycling_DEGs_CON_M_vs_CON_F %>%
  filter(avg_log2FC < -0.5) %>%
  filter(p_val_adj < 0.05) %>%
  arrange(avg_log2FC)
write.csv(all_CON_DEGs_down, "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/blackshaw + kb/analysis on Vim+ subset/DEGs/all_downregulated_CON_M_vs_F.csv")


# just blackshaw con DEGs from this dataset -------------------------------

# Subset for CON, Blackshaw, and Cycling NSPCs all at once
cycling_con_blackshaw_obj <- subset(
  clean_combined_object, 
  subset = treatment == "CON" & dataset == "Blackshaw" & cluster_names == "Cycling NSPCs")

# Confirm identities
Idents(cycling_con_blackshaw_obj) <- "sex"

# Run DE analysis for Male vs. Female within this cluster
cycling_DEGs_CON_M_vs_CON_F <- FindMarkers(
  cycling_con_blackshaw_obj,
  ident.1 = "Male",
  ident.2 = "Female",
  assay = "RNA",
  logfc.threshold = 0.25,
  min.pct = 0.1)

# Save as CSV
outdir <- "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/blackshaw + kb/analysis on Vim+ subset/DEGs/"

write.csv(
  cycling_DEGs_CON_M_vs_CON_F,
  file = paste0(outdir, "CON_DEGs_Male_vs_Female_CyclingNPCs_no_pcutoff_just_blackshaw_degs.csv"))

# Force Seurat to run the DE test explicitly for just the Dlx genes
dlx_stats <- FindMarkers(
  object = cycling_con_blackshaw_obj,
  ident.1 = "Male",
  ident.2 = "Female",
  features = c("Dlx1", "Dlx2", "Dlx5"),
  assay = "RNA",
  min.pct = 0,           
  logfc.threshold = 0)

# Export the Dlx stats dataframe to a CSV file
write.csv(
  x = dlx_stats, 
  file = paste0(outdir, "cycling_nspcs_con_M_vs_F_dlx_stats_just_blackshaw_degs.csv"), 
  row.names = TRUE)

# comparing degs  ---------------------------------------------------------
# manually added "Gene" column title on csv files

# Load both CSV files
blackshaw_data <- read.csv('/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/blackshaw + kb/analysis on Vim+ subset/DEGs/CON_DEGs_Male_vs_Female_CyclingNPCs_no_pcutoff_just_blackshaw_degs.csv')

kodee_data <- read.csv('/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/DEGs/CON_DEGs_Male_vs_Female_CyclingNSPCs_kodee.csv')

# Filter both datasets for raw p-value < 0.05
blackshaw_sig <- subset(blackshaw_data, p_val < 0.05) # 551 / 3381
kodee_sig <- subset(kodee_data, p_val < 0.05) # 2191 / 5891

# Extract the gene names using your column name
blackshaw_sig_genes <- blackshaw_sig$Gene
kodee_sig_genes <- kodee_sig$Gene

# Find the overlap between the two lists
overlapping_genes <- intersect(blackshaw_sig_genes, kodee_sig_genes)

# Output the results
print(paste("Number of overlapping significant genes based on p-value (not adjusted):", length(overlapping_genes))) # 131 genes

# Save as CSV
outdir <- "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/blackshaw + kb/analysis on Vim+ subset/DEGs/"

write.csv(
  overlapping_genes,
  file = paste0(outdir, "CON_DEGs_Male_vs_Female_CyclingNSPCs_overlapping genes_blackshaw only.csv"))
