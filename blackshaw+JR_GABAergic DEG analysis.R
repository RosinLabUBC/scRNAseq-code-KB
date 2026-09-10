# Upreg and Downreg DEGs for all CON M vs CON F Analysis --------

# load in object
merged_object_gad_pos_clean <- readRDS('/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/blackshaw + kb/merged_object_gad_pos_clean.rds')

# Set default assay to RNA
DefaultAssay(merged_object_gad_pos_clean) <- "RNA"

# Normalize and scale RNA assay
merged_object_gad_pos_clean <- NormalizeData(merged_object_gad_pos_clean)
merged_object_gad_pos_clean <- FindVariableFeatures(merged_object_gad_pos_clean)
merged_object_gad_pos_clean <- ScaleData(merged_object_gad_pos_clean, features = row.names(merged_object_gad_pos_clean))

# Join Layers
merged_object_gad_pos_clean <- JoinLayers(merged_object_gad_pos_clean, layers = "data")


# all CON cells -----------------------------------------------------------
# Subset to just the CON cells
con_obj <- subset(merged_object_gad_pos_clean, subset = treatment == "CON")

# Confirm identities
Idents(con_obj) <- "sex"

# calculate DEGs
all_CON_DEGs <- FindMarkers(
  con_obj,
  ident.1 = "Male", # comparing how Male genes change compare to female genes
  ident.2 = "Female",
  assay = "RNA",  
  logfc.threshold = 0.25,
  min.pct = 0.1)

# Save as CSV
outdir <- "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/blackshaw + kb/analysis on GABAergic cells/DEGs/"

write.csv(
  all_CON_DEGs,
  file = paste0(outdir, "CON_DEGs_Male_vs_Female_GABAergic_cells_no_pcutoff.csv"))


# Save all upregulated DEGs (positive log2FC)
all_CON_DEGs_up <- all_CON_DEGs %>%
  filter(avg_log2FC > 0.5) %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC))
write.csv(all_CON_DEGs_up, "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/blackshaw + kb/analysis on GABAergic cells/DEGs/all_upregulated_CON_M_vs_F.csv")

# Save all downregulated DEGs (negative log2FC)
all_CON_DEGs_down <- all_CON_DEGs %>%
  filter(avg_log2FC < -0.5) %>%
  filter(p_val_adj < 0.05) %>%
  arrange(avg_log2FC)
write.csv(all_CON_DEGs_down, "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/blackshaw + kb/analysis on GABAergic cells/DEGs/all_downregulated_CON_M_vs_F.csv")


# just blackshaw con DEGs from this dataset -------------------------------

# Subset for CON, Blackshaw all at once
con_blackshaw_obj <- subset(
  merged_object_gad_pos_clean, 
  subset = treatment == "CON" & dataset == "Blackshaw")

# Confirm identities
Idents(con_blackshaw_obj) <- "sex"

# Run DE analysis for Male vs. Female within this cluster
DEGs_CON_M_vs_CON_F <- FindMarkers(
  con_blackshaw_obj,
  ident.1 = "Male",
  ident.2 = "Female",
  assay = "RNA",
  logfc.threshold = 0.25,
  min.pct = 0.1)

# Save as CSV
outdir <- "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/blackshaw + kb/analysis on GABAergic cells/DEGs/"

write.csv(
  DEGs_CON_M_vs_CON_F,
  file = paste0(outdir, "CON_DEGs_Male_vs_Female_GABAergic_cells_no_pcutoff_just_blackshaw_degs.csv"))

# comparing degs  ---------------------------------------------------------

# manually added "Gene" column title on csv files

# Load both CSV files
blackshaw_data <- read.csv('/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/blackshaw + kb/analysis on GABAergic cells/DEGs/CON_DEGs_Male_vs_Female_GABAergic_cells_no_pcutoff_just_blackshaw_degs.csv')

kodee_data <- read.csv('/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/GABAergic cells/DEGs for comparison/CON_DEGs_Male_vs_Female_GABAergic_cells_kodee.csv')

# Filter both datasets for raw p-value < 0.05
blackshaw_sig <- subset(blackshaw_data, p_val < 0.05) # 368 / 825
kodee_sig <- subset(kodee_data, p_val < 0.05) # 1293 / 1392

# Extract the gene names using your column name
blackshaw_sig_genes <- blackshaw_sig$Gene
kodee_sig_genes <- kodee_sig$Gene

# Find the overlap between the two lists
overlapping_genes <- intersect(blackshaw_sig_genes, kodee_sig_genes)

# Output the results
print(paste("Number of overlapping significant genes based on p-value (not adjusted) in GABAergic cells:", length(overlapping_genes))) # 69 genes

# Save as CSV
outdir <- "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/scRNAseq revisions/blackshaw + kb/analysis on GABAergic cells/DEGs/"

write.csv(
  overlapping_genes,
  file = paste0(outdir, "CON_DEGs_Male_vs_Female_GABAergic_overlapping genes_blackshaw only.csv"))
