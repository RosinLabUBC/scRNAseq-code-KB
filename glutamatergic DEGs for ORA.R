# like in 'DEGs for ORA workflow' ave log2FC > 0.5, < -0.5, adj p-value < 0.05, using all of the DEGs to run pathway analysis

# Upreg and Downreg DEGs for all CON M vs CON F Analysis --------

# load in object
merged_object_vglut_pos_clean <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis of Glutamatergic Cells/objects/merged_object_vglut_pos_clean.rds')

# Set default assay to RNA
DefaultAssay(merged_object_vglut_pos_clean) <- "RNA"

# Normalize and scale RNA assay
merged_object_vglut_pos_clean <- NormalizeData(merged_object_vglut_pos_clean)
merged_object_vglut_pos_clean <- FindVariableFeatures(merged_object_vglut_pos_clean)
merged_object_vglut_pos_clean <- ScaleData(merged_object_vglut_pos_clean, features = row.names(merged_object_vglut_pos_clean))

# Join Layers
merged_object_vglut_pos_clean <- JoinLayers(merged_object_vglut_pos_clean, layers = "data")

# Subset to just the CON cells
con_obj <- subset(merged_object_vglut_pos_clean, subset = treatment == "CON")

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

# Save all upregulated DEGs (positive log2FC)
all_CON_DEGs_top <- all_CON_DEGs %>%
  filter(avg_log2FC > 0.5) %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC))
write.csv(all_CON_DEGs_top, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis of Glutamatergic Cells/Analysis on VGlut+ clean cells/DEGs/DEGs for ORA/all_upregulated_CON_M_vs_F.csv")

# Save all downregulated DEGs (negative log2FC)
all_CON_DEGs_down <- all_CON_DEGs %>%
  filter(avg_log2FC < -0.5) %>%
  filter(p_val_adj < 0.05) %>%
  arrange(avg_log2FC)
write.csv(all_CON_DEGs_down, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis of Glutamatergic Cells/Analysis on VGlut+ clean cells/DEGs/DEGs for ORA/all_downregulated_CON_M_vs_F.csv")


# Upreg and Downreg DEGs for all CON F vs COLD F Analysis ---------------------
female_obj <- subset(merged_object_vglut_pos_clean, subset = sex == "Female")

# Confirm identities
Idents(female_obj) <- "treatment"

# calcualte DEGs
all_Female_DEGs <- FindMarkers(
  female_obj,
  ident.1 = "COLD", # comparing how Cold genes change compare to CON genes
  ident.2 = "CON",
  assay = "RNA",  
  logfc.threshold = 0.25,
  min.pct = 0.1)

# Save all upregulated DEGs (positive log2FC)
all_Female_DEGs_top <- all_Female_DEGs %>%
  filter(avg_log2FC > 0.5) %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC))
write.csv(all_Female_DEGs_top, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis of Glutamatergic Cells/Analysis on VGlut+ clean cells/DEGs/DEGs for ORA/all_Female_upregulated_CON_vs_COLD.csv")

# Save all downregulated DEGs (negative log2FC)
all_Female_DEGs_down <- all_Female_DEGs %>%
  filter(avg_log2FC < -0.5) %>%
  filter(p_val_adj < 0.05) %>%
  arrange(avg_log2FC)
write.csv(all_Female_DEGs_down, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis of Glutamatergic Cells/Analysis on VGlut+ clean cells/DEGs/DEGs for ORA/all_Female_downregulated_CON_vs_COLD.csv")


# Upreg and Downreg DEGs for all CON M vs COLD M Analysis ---------------------
male_obj <- subset(merged_object_vglut_pos_clean, subset = sex == "Male")

# Confirm identities
Idents(male_obj) <- "treatment"

# calcualte DEGs
all_Male_DEGs <- FindMarkers(
  male_obj,
  ident.1 = "COLD", # comparing how Cold genes change compare to CON genes
  ident.2 = "CON",
  assay = "RNA",  
  logfc.threshold = 0.25,
  min.pct = 0.1)

# Save all upregulated DEGs (positive log2FC)
all_Male_DEGs_top <- all_Male_DEGs %>%
  filter(avg_log2FC > 0.5) %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC))
write.csv(all_Male_DEGs_top, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis of Glutamatergic Cells/Analysis on VGlut+ clean cells/DEGs/DEGs for ORA/all_Male_upregulated_CON_vs_COLD.csv")

# Save all downregulated DEGs (negative log2FC)
all_Male_DEGs_down <- all_Male_DEGs %>%
  filter(avg_log2FC < -0.5) %>%
  filter(p_val_adj < 0.05) %>%
  arrange(avg_log2FC)
write.csv(all_Male_DEGs_down, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis of Glutamatergic Cells/Analysis on VGlut+ clean cells/DEGs/DEGs for ORA/all_Male_downregulated_CON_vs_COLD.csv")
