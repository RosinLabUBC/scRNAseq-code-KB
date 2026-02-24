# following 'DEGs for ORA workflow.R' from orignal analysis

# want DEGs list and volcano plots for CON M vs CON F, CON F vs COLD F
# Run DE analysis for CON F vs COLD F and CON F vs CON M (from original dataset)

library(Seurat)
library(dplyr)
library(patchwork)
library(clustree)
library(pheatmap)
library(EnhancedVolcano)

# load in clean combined object -------------------------------------------

clean_combined_object <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/saved objects/clean_combined_object_CON&COLD_F_verification.rds')

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
outdir <- "/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/Analysis on Vim+ subsetted cells/Cycling DEGs/"

write.csv(
  cycling_DEGs_CON_M_vs_CON_F,
  file = paste0(outdir, "CON_DEGs_Male_vs_Female_CyclingNSPCs_no_pcutoff.csv"))

saveRDS(
  cycling_DEGs_CON_M_vs_CON_F,
  file = paste0(outdir, "CON_DEGs_Male_vs_Female_CyclingNSPCs_no_pcutoff.rds"))


# DEGs on Cycling NSPCs in CON F vs COLD F ----------------------------------------------
# Subset to just the Female cells
female_obj <- subset(clean_combined_object, subset = sex == "Female")

# Subset to just the cycling NSPCs
cycling_female_obj <- subset(female_obj, subset = cluster_names == "Cycling NSPCs")

# Confirm identities
Idents(cycling_female_obj) <- "treatment"

# Run DE analysis for COLD F vs CON F within this cluster
cycling_DEGs_CON_F_vs_COLD_F <- FindMarkers(
  cycling_female_obj,
  ident.1 = "COLD", # comparing how Cold genes change compare to CON genes
  ident.2 = "CON",
  assay = "RNA",
  logfc.threshold = 0.25,
  min.pct = 0.1)

# Save as CSV
outdir <- "/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/Analysis on Vim+ subsetted cells/Cycling DEGs/"

write.csv(
  cycling_DEGs_CON_F_vs_COLD_F,
  file = paste0(outdir, "Female_DEGs_CON_vs_COLD_CyclingNSPCs_no_pcutoff.csv"))

saveRDS(
  cycling_DEGs_CON_F_vs_COLD_F,
  file = paste0(outdir, "Female_DEGs_CON_vs_COLD_CyclingNSPCs_no_pcutoff.rds"))


# overlapping genes from this dataset -------------------------------------
# data taken from sheets from 'volcano plot' folder under CON F Verification 2 folder

# CON M and COLD F --------------------------------------------------------
# Load upregulated genes
con_m_vs_f <- read.csv("/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/Analysis on Vim+ subsetted cells/volcano plot/upregulated_CON_M_vs_F.csv")
con_f_vs_cold_f <- read.csv("/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/Analysis on Vim+ subsetted cells/volcano plot/upregulated_CON_F_vs_COLD_F.csv")

# Extract gene names
genes_con_m <- con_m_vs_f$gene
genes_cold_f <- con_f_vs_cold_f$gene

# Count total genes
length(genes_con_m)       # 226
length(genes_cold_f)      # 281

# Find intersection
shared_genes <- intersect(genes_con_m, genes_cold_f)
length(shared_genes)      # 147

# Convert to a data frame without sorting
shared_genes_df <- data.frame(gene = shared_genes)

# Export to CSV in your desired folder
write.csv(shared_genes_df, 
          "/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/Analysis on Vim+ subsetted cells/Cycling DEGs/shared_upregulated_genes_CON_M_and_COLD_F.csv",
          row.names = FALSE)

# Load downregulated genes
con_m_vs_f <- read.csv("/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/Analysis on Vim+ subsetted cells/volcano plot/downregulated_CON_M_vs_F.csv")
con_f_vs_cold_f <- read.csv("/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/Analysis on Vim+ subsetted cells/volcano plot/downregulated_CON_F_vs_COLD_F.csv")

# Extract gene names
genes_con_m <- con_m_vs_f$gene
genes_cold_f <- con_f_vs_cold_f$gene

# Count total genes
length(genes_con_m)       # 52
length(genes_cold_f)      # 1736

# Find intersection
shared_genes <- intersect(genes_con_m, genes_cold_f)
length(shared_genes)      # 44

# Convert to a data frame without sorting
shared_genes_df <- data.frame(gene = shared_genes)

# Export to CSV in your desired folder
write.csv(shared_genes_df, 
          "/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/Analysis on Vim+ subsetted cells/Cycling DEGs/shared_downregulated_genes_CON_M_vs_F_and_CON_F_vs_COLD_F.csv",
          row.names = FALSE)


# overlapping genes with original dataset DEG analysis --------------------
# Load upregulated genes
original_shared_genes <- read.csv("/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/volcano plots2/check/shared_upregulated_genes_CON_M_vs_F_and_CON_F_vs_COLD_F.csv")
f_verif_shared_genes <- read.csv("/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/Analysis on Vim+ subsetted cells/Cycling DEGs/shared_upregulated_genes_CON_M_and_COLD_F.csv")

# Extract gene names
genes_orig <- original_shared_genes$gene
genes_f_verif <- f_verif_shared_genes$gene

# Count total genes
length(genes_orig)       # 51
length(genes_f_verif)      # 147

# Find intersection
shared_genes <- intersect(genes_orig, genes_f_verif)
length(shared_genes)      # 37

# Convert to a data frame without sorting
shared_genes_df <- data.frame(gene = shared_genes)

# Export to CSV in your desired folder
write.csv(shared_genes_df, 
          "/Users/kodeebao/Documents/Strategy 2 Final Analysis/CON F Verification 2/Analysis on Vim+ subsetted cells/Cycling DEGs/shared_upregulated_genes_orig_and_f_verif.csv",
          row.names = FALSE)


