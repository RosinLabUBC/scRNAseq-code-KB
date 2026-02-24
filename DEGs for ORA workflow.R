# want DEGs list and volcano plots for CON M vs CON F, CON M vs COLD M, CON F vs COLD F

library(Seurat)
library(dplyr)
library(patchwork)
library(clustree)
library(pheatmap)
library(EnhancedVolcano)

# load in object
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

# Confirm identities
Idents(con_obj) <- "sex"

# DEGs on all CON M vs CON F ----------------------------------------------

# calcualte DEGs
all_CON_DEGs <- FindMarkers(
  con_obj,
  ident.1 = "Male", # comparing how Male genes change compare to female genes
  ident.2 = "Female",
  assay = "RNA",  
  logfc.threshold = 0.25,
  min.pct = 0.1)

# Save
write.csv(all_CON_DEGs, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/DEGs/all_CON_DEGs_Male_vs_Female_(no_p-value_cutoff).csv")

saveRDS(all_CON_DEGs, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/DEG objects/all_CON_DEGs_Male_vs_Female_(no_p-value_cutoff).rds")


# DEGs on individual clusters in CON M vs CON F ----------------------------------------------
# Set the identities to the renamed clusters
Idents(con_obj) <- "cluster_names"

# Create a list of the clusters you want to analyze
CON_DEGs_clusters_to_analyze <- c("Neuronal Precursors (Ascl1+)", "Glial Precursors", "Cycling NSPCs", "NPCs-Glia", "NPCs-Neuronal", "Neuronal Precursors (Neurog2+)")

# Initialize a list to store DE results for each cluster
all_cluster_CON_DEGs <- list()

# Loop through each cluster to run DE analysis
for (cluster in CON_DEGs_clusters_to_analyze) {
  # Subset the object to include only cells from the current cluster
  cluster_cells_CON <- subset(con_obj, cluster_names == cluster)
  
  # Set identity to sex (Male vs Female)
  Idents(cluster_cells_CON) <- "sex"
  
  # Run DE analysis for Male vs. Female within this cluster
  DEGs_cluster_CON <- FindMarkers(
    cluster_cells_CON,
    ident.1 = "Male",  
    ident.2 = "Female",   
    assay = "RNA",         
    logfc.threshold = 0.25,
    min.pct = 0.1)
  
  # Store the DEGs for this cluster
  all_cluster_CON_DEGs[[cluster]] <- DEGs_cluster_CON
  
  # Clean the cluster name for file naming
  clean_cluster_name <- gsub("[\\(\\)+]", "", cluster)  # remove () and +
  clean_cluster_name <- gsub(" ", "_", clean_cluster_name)  # replace spaces with _
  
  # Save as CSV
  write.csv(
    DEGs_cluster_CON,
    paste0("/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/DEGs/",
           clean_cluster_name, "_CON_DEGs_Male_vs_Female_(no_p-value_cutoff).csv"))
  
  # Save as RDS
  saveRDS(
    DEGs_cluster_CON,
    paste0("/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/DEG objects/",
           clean_cluster_name, "_CON_DEGs_Male_vs_Female_(no_p-value_cutoff).rds"))}



# DEGs on all CON M vs COLD M -----------------------------------------------------

# Subset out male cells
male_obj <- subset(clean_combined_object, subset = sex == "Male")

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

# Save
write.csv(all_Male_DEGs, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/DEGs/all_Male_DEGs_CON_vs_COLD_(no_p-value_cutoff).csv")

saveRDS(all_Male_DEGs, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/DEG objects/all_Male_DEGs_CON_vs_COLD_(no_p-value_cutoff).rds")


# DEGs on individual clusters CON M vs COLD M -----------------------------------------------------

# Set the identities to the renamed clusters
Idents(male_obj) <- "cluster_names"

# Create a list of the clusters you want to analyze
Male_DEGs_clusters_to_analyze <- c("Neuronal Precursors (Ascl1+)", "Glial Precursors", "Cycling NSPCs", "NPCs-Glia", "NPCs-Neuronal", "Neuronal Precursors (Neurog2+)")

# Initialize a list to store DE results for each cluster
all_cluster_Male_DEGs <- list()

# Loop through each cluster to run DE analysis
for (cluster in Male_DEGs_clusters_to_analyze) {
  # Subset the object to include only cells from the current cluster
  cluster_cells_Male <- subset(male_obj, cluster_names == cluster)
  
  # Set identity to treatment (CON vs COLD)
  Idents(cluster_cells_Male) <- "treatment"
  
  # Run DE analysis for COLD vs. CON within this cluster
  DEGs_cluster_Male <- FindMarkers(
    cluster_cells_Male,
    ident.1 = "COLD",  
    ident.2 = "CON",   
    assay = "RNA",         
    logfc.threshold = 0.25,
    min.pct = 0.1)
  
  # Store the DEGs for this cluster
  all_cluster_Male_DEGs[[cluster]] <- DEGs_cluster_Male
  
  # Clean the cluster name for file naming
  clean_cluster_name <- gsub("[\\(\\)+]", "", cluster)  # remove () and +
  clean_cluster_name <- gsub(" ", "_", clean_cluster_name)  # replace spaces with _
  
  # Save as CSV
  write.csv(
    DEGs_cluster_Male,
    paste0("/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/DEGs/",
           clean_cluster_name, "_Male_DEGs_CON_vs_COLD_(no_p-value_cutoff).csv"))
  
  # Save as RDS
  saveRDS(
    DEGs_cluster_Male,
    paste0("/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/DEG objects/",
           clean_cluster_name, "_Male_DEGs_CON_vs_COLD_(no_p-value_cutoff).rds"))}


# DEGs on all CON F vs COLD F -----------------------------------------------------

# Subset out female cells
female_obj <- subset(clean_combined_object, subset = sex == "Female")

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

# Save
write.csv(all_Female_DEGs, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/DEGs/all_Female_DEGs_CON_vs_COLD_(no_p-value_cutoff).csv")

saveRDS(all_Female_DEGs, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/DEG objects/all_Female_DEGs_CON_vs_COLD_(no_p-value_cutoff).rds")


# DEGs on individual clusters CON F vs COLD F -----------------------------------------------------

# Set the identities to the renamed clusters
Idents(female_obj) <- "cluster_names"

# Create a list of the clusters you want to analyze
Female_DEGs_clusters_to_analyze <- c("Neuronal Precursors (Ascl1+)", "Glial Precursors", "Cycling NSPCs", "NPCs-Glia", "NPCs-Neuronal", "Neuronal Precursors (Neurog2+)")

# Initialize a list to store DE results for each cluster
all_cluster_Female_DEGs <- list()

# Loop through each cluster to run DE analysis
for (cluster in Female_DEGs_clusters_to_analyze) {
  # Subset the object to include only cells from the current cluster
  cluster_cells_Female <- subset(female_obj, cluster_names == cluster)
  
  # Set identity to treatment (CON vs COLD)
  Idents(cluster_cells_Female) <- "treatment"
  
  # Run DE analysis for COLD vs. CON within this cluster
  DEGs_cluster_Female <- FindMarkers(
    cluster_cells_Female,
    ident.1 = "COLD",  
    ident.2 = "CON",   
    assay = "RNA",         
    logfc.threshold = 0.25,
    min.pct = 0.1)
  
  # Store the DEGs for this cluster
  all_cluster_Female_DEGs[[cluster]] <- DEGs_cluster_Female
  
  # Clean the cluster name for file naming
  clean_cluster_name <- gsub("[\\(\\)+]", "", cluster)  # remove () and +
  clean_cluster_name <- gsub(" ", "_", clean_cluster_name)  # replace spaces with _
  
  # Save as CSV
  write.csv(
    DEGs_cluster_Female,
    paste0("/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/DEGs/",
           clean_cluster_name, "_Female_DEGs_CON_vs_COLD_(no_p-value_cutoff).csv"))
  
  # Save as RDS
  saveRDS(
    DEGs_cluster_Female,
    paste0("/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/DEG objects/",
           clean_cluster_name, "_Female_DEGs_CON_vs_COLD_(no_p-value_cutoff).rds"))}


# Upreg and Downreg DEGs for Analysis ----------------------------------------------------

# before, for initial just female we did avg_log2FC > 0, < 0, p_value_adj < 0.05

# ave log2FC > 0.5, < -0.5, adj p-value < 0.05, using all of the DEGs to run pathway analysis (will note down how many are up and down regulated for each cluster), talked to JR

# USE CODE BELOW

# Upreg and Downreg DEGs for all CON M vs CON F Analysis --------

# load in object
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
write.csv(all_CON_DEGs_top, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/DEGs for Analysis/all_upregulated_CON_M_vs_F.csv")

# Save all downregulated DEGs (negative log2FC)
all_CON_DEGs_down <- all_CON_DEGs %>%
  filter(avg_log2FC < -0.5) %>%
  filter(p_val_adj < 0.05) %>%
  arrange(avg_log2FC)
write.csv(all_CON_DEGs_down, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/DEGs for Analysis/all_downregulated_CON_M_vs_F.csv")


# Upreg and Downreg DEGs for individual clusters CON M vs CON F Analysis ---------------------------------
# Set the identities to the renamed clusters
Idents(con_obj) <- "cluster_names"

# Create a list of the clusters you want to analyze
CON_DEGs_clusters_to_analyze <- c("Neuronal Precursors (Ascl1+)", "Glial Precursors", "Cycling NSPCs", "NPCs-Glia", "NPCs-Neuronal", "Neuronal Precursors (Neurog2+)")

# Loop over clusters
for (cluster in CON_DEGs_clusters_to_analyze) {
  
  # Subset to just cells in this cluster
  cluster_obj <- subset(con_obj, idents = cluster)
  
  # Set identity to sex for DEGs
  Idents(cluster_obj) <- "sex"
  
  # Run differential expression (Male vs Female)
  degs <- FindMarkers(
    cluster_obj,
    ident.1 = "Male",
    ident.2 = "Female",
    assay = "RNA",
    logfc.threshold = 0.25,
    min.pct = 0.1)
  
  # Clean cluster name for filenames
  clean_cluster <- gsub(" ", "_", gsub("[\\(\\)\\+]", "", cluster))
  
  # Filter and save upregulated DEGs
  upregulated <- degs %>%
    filter(avg_log2FC > 0.5) %>%
    filter(p_val_adj < 0.05) %>%
    arrange(desc(avg_log2FC))
  
  write.csv(
    upregulated,
    paste0("/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/DEGs for Analysis/upregulated_", clean_cluster, "_CON_M_vs_F.csv"))
  
  # Filter and save downregulated DEGs
  downregulated <- degs %>%
    filter(avg_log2FC < -0.5) %>%
    filter(p_val_adj < 0.05) %>%
    arrange(avg_log2FC)
  
  write.csv(
    downregulated,
    paste0("/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/DEGs for Analysis/downregulated_", clean_cluster, "_CON_M_vs_F.csv"))}

# Upreg and Downreg DEGs for all CON M vs COLD M Analysis ---------------------
# Subset out male cells
male_obj <- subset(clean_combined_object, subset = sex == "Male")

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
write.csv(all_Male_DEGs_top, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/DEGs for Analysis/male only/all_Male_upregulated_CON_vs_COLD.csv")

# Save all downregulated DEGs (negative log2FC)
all_Male_DEGs_down <- all_Male_DEGs %>%
  filter(avg_log2FC < -0.5) %>%
  filter(p_val_adj < 0.05) %>%
  arrange(avg_log2FC)
write.csv(all_Male_DEGs_down, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/DEGs for Analysis/male only/all_Male_downregulated_CON_vs_COLD.csv")


# Upreg and Downreg DEGs for individual clusters CON M vs COLD M Analysis ---------------------
Idents(male_obj) <- "cluster_names"

# Create a list of the clusters you want to analyze
Male_DEGs_clusters_to_analyze <- c("Neuronal Precursors (Ascl1+)", "Glial Precursors", "Cycling NSPCs", "NPCs-Glia", "NPCs-Neuronal", "Neuronal Precursors (Neurog2+)")

# Loop over clusters
for (cluster in Male_DEGs_clusters_to_analyze) {
  # Subset cluster
  cluster_cells_Male <- subset(male_obj, idents = cluster)
  
  # Set identity to treatment (COLD vs CON)
  Idents(cluster_cells_Male) <- "treatment"
  
  # Find DEGs
  degs <- FindMarkers(
    cluster_cells_Male,
    ident.1 = "COLD",
    ident.2 = "CON",
    assay = "RNA",
    logfc.threshold = 0.25,
    min.pct = 0.1)
  
  # Clean cluster name
  clean_cluster <- gsub(" ", "_", gsub("[\\(\\)\\+]", "", cluster))
  
  # Save upregulated
  up <- degs %>%
    filter(avg_log2FC > 0.5) %>%
    filter(p_val_adj < 0.05) %>%
    arrange(desc(avg_log2FC))
  write.csv(up, paste0("/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/DEGs for Analysis/male only/upregulated_", clean_cluster, "_COLD_vs_CON_M.csv"))
  
  # Save downregulated
  down <- degs %>%
    filter(avg_log2FC < -0.5) %>%
    filter(p_val_adj < 0.05) %>%
    arrange(avg_log2FC)
  write.csv(down, paste0("/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/DEGs for Analysis/male only/downregulated_", clean_cluster, "_COLD_vs_CON_M.csv"))}


# Upreg and Downreg DEGs for all CON F vs COLD F Analysis ---------------------
female_obj <- subset(clean_combined_object, subset = sex == "Female")

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
write.csv(all_Female_DEGs_top, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/DEGs for Analysis/female only/all_Female_upregulated_CON_vs_COLD.csv")

# Save all downregulated DEGs (negative log2FC)
all_Female_DEGs_down <- all_Female_DEGs %>%
  filter(avg_log2FC < -0.5) %>%
  filter(p_val_adj < 0.05) %>%
  arrange(avg_log2FC)
write.csv(all_Female_DEGs_down, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/DEGs for Analysis/female only/all_Female_downregulated_CON_vs_COLD.csv")


# Upreg and Downreg DEGs for individual clusters CON F vs COLD F Analysis ---------------------
Idents(female_obj) <- "cluster_names"

# Create a list of the clusters you want to analyze
Female_DEGs_clusters_to_analyze <- c("Neuronal Precursors (Ascl1+)", "Glial Precursors", "Cycling NSPCs", "NPCs-Glia", "NPCs-Neuronal", "Neuronal Precursors (Neurog2+)")

# Loop over clusters
for (cluster in Female_DEGs_clusters_to_analyze) {
  # Subset cluster
  cluster_cells_Female <- subset(female_obj, idents = cluster)
  
  # Set identity to treatment (COLD vs CON)
  Idents(cluster_cells_Female) <- "treatment"
  
  # Find DEGs
  degs <- FindMarkers(
    cluster_cells_Female,
    ident.1 = "COLD",
    ident.2 = "CON",
    assay = "RNA",
    logfc.threshold = 0.25,
    min.pct = 0.1)
  
  # Clean cluster name
  clean_cluster <- gsub(" ", "_", gsub("[\\(\\)\\+]", "", cluster))
  
  # Save upregulated
  up <- degs %>%
    filter(avg_log2FC > 0.5) %>%
    filter(p_val_adj < 0.05) %>%
    arrange(desc(avg_log2FC))
  write.csv(up, paste0("/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/DEGs for Analysis/female only/upregulated_", clean_cluster, "_COLD_vs_CON_F.csv"))
  
  # Save downregulated
  down <- degs %>%
    filter(avg_log2FC < -0.5) %>%
    filter(p_val_adj < 0.05) %>%
    arrange(avg_log2FC)
  write.csv(down, paste0("/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/DEGs for Analysis/female only/downregulated_", clean_cluster, "_COLD_vs_CON_F.csv"))}

