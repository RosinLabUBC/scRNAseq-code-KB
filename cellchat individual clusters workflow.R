library(CellChat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(Seurat)
library(future)


# load NSC objects but need to change to RNA assay from SCT and normalize ----------
# load in object
clean_combined_object <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/clean_combined_object_CON&COLD_strat2.rds')

# Set default assay to RNA
DefaultAssay(clean_combined_object) <- "RNA"

# Normalize and scale RNA assay
clean_combined_object <- NormalizeData(clean_combined_object)
clean_combined_object <- FindVariableFeatures(clean_combined_object)
clean_combined_object <- ScaleData(clean_combined_object, features = row.names(clean_combined_object))

# subsetting
# Female + CON
female_con_NSC_obj <- subset(clean_combined_object, subset = sex == "Female" & treatment == "CON")

# Female + COLD
female_cold_NSC_obj <- subset(clean_combined_object, subset = sex == "Female" & treatment == "COLD")

# Male + CON
male_con_NSC_obj <- subset(clean_combined_object, subset = sex == "Male" & treatment == "CON")

# Male + COLD
male_cold_NSC_obj <- subset(clean_combined_object, subset = sex == "Male" & treatment == "COLD")


# load microglia objects default is already RNA assay and has been normalized ----------
# load in microglia objects
female_con_mg_obj <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/microglia data/female_con_mg_obj.rds')
female_cold_mg_obj <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/microglia data/female_cold_mg_obj.rds')
male_con_mg_obj <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/microglia data/male_con_mg_obj.rds')
male_cold_mg_obj <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/microglia data/male_cold_mg_obj.rds')

# Keep cluster annotations as celltype for NSCs (instead of generic "NSC")
female_con_NSC_obj$celltype <- female_con_NSC_obj$cluster_names
female_cold_NSC_obj$celltype <- female_cold_NSC_obj$cluster_names
male_con_NSC_obj$celltype <- male_con_NSC_obj$cluster_names
male_cold_NSC_obj$celltype <- male_cold_NSC_obj$cluster_names

# Keep "Microglia" for microglia object
female_con_mg_obj$celltype <- "Microglia"
female_cold_mg_obj$celltype <- "Microglia"
male_con_mg_obj$celltype <- "Microglia"
male_cold_mg_obj$celltype <- "Microglia"

# Load the mouse database
CellChatDB <- CellChatDB.mouse


# Cell Chat for Female CON NSCs (with clusters) + Female CON Microglia ------------------------
# Expression matrices
female_con_NSC_mtx <- GetAssayData(female_con_NSC_obj, assay = "RNA", slot = "data")
female_con_mg_mtx  <- GetAssayData(female_con_mg_obj, assay = "RNA", slot = "data")

# Convert to data.frame for merging
df_fem_con_NSC <- as.data.frame(female_con_NSC_mtx)
df_fem_con_mg  <- as.data.frame(female_con_mg_mtx)

# Add gene column
df_fem_con_NSC$gene <- rownames(df_fem_con_NSC)
df_fem_con_mg$gene  <- rownames(df_fem_con_mg)

# Merge by gene
merged_fem_con <- merge(df_fem_con_NSC, df_fem_con_mg, by = "gene", all = TRUE)
merged_fem_con[is.na(merged_fem_con)] <- 0
rownames(merged_fem_con) <- merged_fem_con$gene
merged_fem_con$gene <- NULL
merged_matrix_fem_con <- as.matrix(merged_fem_con)

# Merge metadata (NSC + MG)
merged_meta_fem_con <- bind_rows(female_con_NSC_obj@meta.data, female_con_mg_obj@meta.data)

# Create CellChat object (cluster-level NSCs vs. Microglia)
female_con_NSCclusters_mg_cellchat <- createCellChat(
  object = merged_matrix_fem_con, 
  meta   = merged_meta_fem_con, 
  group.by = "celltype")

# Assign mouse LR database
female_con_NSCclusters_mg_cellchat@DB <- CellChatDB.mouse

# Standard CellChat pipeline
female_con_NSCclusters_mg_cellchat <- subsetData(female_con_NSCclusters_mg_cellchat)
female_con_NSCclusters_mg_cellchat <- identifyOverExpressedGenes(female_con_NSCclusters_mg_cellchat)
female_con_NSCclusters_mg_cellchat <- identifyOverExpressedInteractions(female_con_NSCclusters_mg_cellchat)
female_con_NSCclusters_mg_cellchat <- smoothData(female_con_NSCclusters_mg_cellchat, adj = PPI.mouse)

female_con_NSCclusters_mg_cellchat <- computeCommunProb(female_con_NSCclusters_mg_cellchat)
female_con_NSCclusters_mg_cellchat <- filterCommunication(female_con_NSCclusters_mg_cellchat, min.cells = 10)
female_con_NSCclusters_mg_cellchat <- computeCommunProbPathway(female_con_NSCclusters_mg_cellchat)
female_con_NSCclusters_mg_cellchat <- aggregateNet(female_con_NSCclusters_mg_cellchat)

# Extract results
lr_interactions <- subsetCommunication(female_con_NSCclusters_mg_cellchat, slot.name = "net")
pathway_interactions <- subsetCommunication(female_con_NSCclusters_mg_cellchat, slot.name = "netP")

# Save
write.csv(lr_interactions, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Ligand-Receptor Analysis/clusters/Female CON/female_con_NSCclusters_mg_LR_interactions.csv", row.names = FALSE)
#1743 results


saveRDS(female_con_NSCclusters_mg_cellchat, file = "/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/CellChat Objects/NSC clusters vs microglia/female_con_NSCclusters_mg_cellchat.rds")


# Cell Chat for Female COLD NSCs (with clusters) + Female COLD Microglia ------------------------
# Expression matrices
female_cold_NSC_mtx <- GetAssayData(female_cold_NSC_obj, assay = "RNA", slot = "data")
female_cold_mg_mtx  <- GetAssayData(female_cold_mg_obj, assay = "RNA", slot = "data")

# Convert to data.frame for merging
df_fem_cold_NSC <- as.data.frame(female_cold_NSC_mtx)
df_fem_cold_mg  <- as.data.frame(female_cold_mg_mtx)

# Add gene column
df_fem_cold_NSC$gene <- rownames(df_fem_cold_NSC)
df_fem_cold_mg$gene  <- rownames(df_fem_cold_mg)

# Merge by gene
merged_fem_cold <- merge(df_fem_cold_NSC, df_fem_cold_mg, by = "gene", all = TRUE)
merged_fem_cold[is.na(merged_fem_cold)] <- 0
rownames(merged_fem_cold) <- merged_fem_cold$gene
merged_fem_cold$gene <- NULL
merged_matrix_fem_cold <- as.matrix(merged_fem_cold)

# Merge metadata (NSC + MG)
merged_meta_fem_cold <- bind_rows(female_cold_NSC_obj@meta.data, female_cold_mg_obj@meta.data)

# Create CellChat object (cluster-level NSCs vs. Microglia)
female_cold_NSCclusters_mg_cellchat <- createCellChat(
  object = merged_matrix_fem_cold, 
  meta   = merged_meta_fem_cold, 
  group.by = "celltype")

# Assign mouse LR database
female_cold_NSCclusters_mg_cellchat@DB <- CellChatDB.mouse

# Standard CellChat pipeline
female_cold_NSCclusters_mg_cellchat <- subsetData(female_cold_NSCclusters_mg_cellchat)
female_cold_NSCclusters_mg_cellchat <- identifyOverExpressedGenes(female_cold_NSCclusters_mg_cellchat)
female_cold_NSCclusters_mg_cellchat <- identifyOverExpressedInteractions(female_cold_NSCclusters_mg_cellchat)
female_cold_NSCclusters_mg_cellchat <- smoothData(female_cold_NSCclusters_mg_cellchat, adj = PPI.mouse)

female_cold_NSCclusters_mg_cellchat <- computeCommunProb(female_cold_NSCclusters_mg_cellchat)
female_cold_NSCclusters_mg_cellchat <- filterCommunication(female_cold_NSCclusters_mg_cellchat, min.cells = 10)
female_cold_NSCclusters_mg_cellchat <- computeCommunProbPathway(female_cold_NSCclusters_mg_cellchat)
female_cold_NSCclusters_mg_cellchat <- aggregateNet(female_cold_NSCclusters_mg_cellchat)

# Extract results
lr_interactions <- subsetCommunication(female_cold_NSCclusters_mg_cellchat, slot.name = "net")
pathway_interactions <- subsetCommunication(female_cold_NSCclusters_mg_cellchat, slot.name = "netP")

# Save
write.csv(lr_interactions, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Ligand-Receptor Analysis/clusters/Female COLD/female_cold_NSCclusters_mg_LR_interactions.csv", row.names = FALSE)
#1111 results

saveRDS(female_cold_NSCclusters_mg_cellchat, file = "/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/CellChat Objects/NSC clusters vs microglia/female_cold_NSCclusters_mg_cellchat.rds")


# Cell Chat for Male CON NSCs (with clusters) + Male CON Microglia ------------------------
# Expression matrices
male_con_NSC_mtx <- GetAssayData(male_con_NSC_obj, assay = "RNA", slot = "data")
male_con_mg_mtx  <- GetAssayData(male_con_mg_obj, assay = "RNA", slot = "data")

# Convert to data.frame for merging
df_male_con_NSC <- as.data.frame(male_con_NSC_mtx)
df_male_con_mg  <- as.data.frame(male_con_mg_mtx)

# Add gene column
df_male_con_NSC$gene <- rownames(df_male_con_NSC)
df_male_con_mg$gene  <- rownames(df_male_con_mg)

# Merge by gene
merged_male_con <- merge(df_male_con_NSC, df_male_con_mg, by = "gene", all = TRUE)
merged_male_con[is.na(merged_male_con)] <- 0
rownames(merged_male_con) <- merged_male_con$gene
merged_male_con$gene <- NULL
merged_matrix_male_con <- as.matrix(merged_male_con)

# Merge metadata (NSC + MG)
merged_meta_male_con <- bind_rows(male_con_NSC_obj@meta.data, male_con_mg_obj@meta.data)

# Create CellChat object (cluster-level NSCs vs. Microglia)
male_con_NSCclusters_mg_cellchat <- createCellChat(
  object = merged_matrix_male_con, 
  meta   = merged_meta_male_con, 
  group.by = "celltype")

# Assign mouse LR database
male_con_NSCclusters_mg_cellchat@DB <- CellChatDB.mouse

# Standard CellChat pipeline
male_con_NSCclusters_mg_cellchat <- subsetData(male_con_NSCclusters_mg_cellchat)
male_con_NSCclusters_mg_cellchat <- identifyOverExpressedGenes(male_con_NSCclusters_mg_cellchat)
male_con_NSCclusters_mg_cellchat <- identifyOverExpressedInteractions(male_con_NSCclusters_mg_cellchat)
male_con_NSCclusters_mg_cellchat <- smoothData(male_con_NSCclusters_mg_cellchat, adj = PPI.mouse)

male_con_NSCclusters_mg_cellchat <- computeCommunProb(male_con_NSCclusters_mg_cellchat)
male_con_NSCclusters_mg_cellchat <- filterCommunication(male_con_NSCclusters_mg_cellchat, min.cells = 10)
male_con_NSCclusters_mg_cellchat <- computeCommunProbPathway(male_con_NSCclusters_mg_cellchat)
male_con_NSCclusters_mg_cellchat <- aggregateNet(male_con_NSCclusters_mg_cellchat)

# Extract results
lr_interactions <- subsetCommunication(male_con_NSCclusters_mg_cellchat, slot.name = "net")
pathway_interactions <- subsetCommunication(male_con_NSCclusters_mg_cellchat, slot.name = "netP")

# Save
write.csv(lr_interactions, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Ligand-Receptor Analysis/clusters/Male CON/male_con_NSCclusters_mg_LR_interactions.csv", row.names = FALSE)
#2134 results

saveRDS(male_con_NSCclusters_mg_cellchat, file = "/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/CellChat Objects/NSC clusters vs microglia/male_con_NSCclusters_mg_cellchat.rds")


# Cell Chat for Male COLD NSCs (with clusters) + Male COLD Microglia ------------------------
# Expression matrices
male_cold_NSC_mtx <- GetAssayData(male_cold_NSC_obj, assay = "RNA", slot = "data")
male_cold_mg_mtx  <- GetAssayData(male_cold_mg_obj, assay = "RNA", slot = "data")

# Convert to data.frame for merging
df_male_cold_NSC <- as.data.frame(male_cold_NSC_mtx)
df_male_cold_mg  <- as.data.frame(male_cold_mg_mtx)

# Add gene column
df_male_cold_NSC$gene <- rownames(df_male_cold_NSC)
df_male_cold_mg$gene  <- rownames(df_male_cold_mg)

# Merge by gene
merged_male_cold <- merge(df_male_cold_NSC, df_male_cold_mg, by = "gene", all = TRUE)
merged_male_cold[is.na(merged_male_cold)] <- 0
rownames(merged_male_cold) <- merged_male_cold$gene
merged_male_cold$gene <- NULL
merged_matrix_male_cold <- as.matrix(merged_male_cold)

# Merge metadata (NSC + MG)
merged_meta_male_cold <- bind_rows(male_cold_NSC_obj@meta.data, male_cold_mg_obj@meta.data)

# Create CellChat object (cluster-level NSCs vs. Microglia)
male_cold_NSCclusters_mg_cellchat <- createCellChat(
  object = merged_matrix_male_cold, 
  meta   = merged_meta_male_cold, 
  group.by = "celltype")

# Assign mouse LR database
male_cold_NSCclusters_mg_cellchat@DB <- CellChatDB.mouse

# Standard CellChat pipeline
male_cold_NSCclusters_mg_cellchat <- subsetData(male_cold_NSCclusters_mg_cellchat)
male_cold_NSCclusters_mg_cellchat <- identifyOverExpressedGenes(male_cold_NSCclusters_mg_cellchat)
male_cold_NSCclusters_mg_cellchat <- identifyOverExpressedInteractions(male_cold_NSCclusters_mg_cellchat)
male_cold_NSCclusters_mg_cellchat <- smoothData(male_cold_NSCclusters_mg_cellchat, adj = PPI.mouse)

male_cold_NSCclusters_mg_cellchat <- computeCommunProb(male_cold_NSCclusters_mg_cellchat)
male_cold_NSCclusters_mg_cellchat <- filterCommunication(male_cold_NSCclusters_mg_cellchat, min.cells = 10)
male_cold_NSCclusters_mg_cellchat <- computeCommunProbPathway(male_cold_NSCclusters_mg_cellchat)
male_cold_NSCclusters_mg_cellchat <- aggregateNet(male_cold_NSCclusters_mg_cellchat)

# Extract results
lr_interactions <- subsetCommunication(male_cold_NSCclusters_mg_cellchat, slot.name = "net")
pathway_interactions <- subsetCommunication(male_cold_NSCclusters_mg_cellchat, slot.name = "netP")

# Save
write.csv(lr_interactions, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Ligand-Receptor Analysis/clusters/Male COLD/male_cold_NSCclusters_mg_LR_interactions.csv", row.names = FALSE)
#1259 results

saveRDS(male_cold_NSCclusters_mg_cellchat, file = "/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/CellChat Objects/NSC clusters vs microglia/male_cold_NSCclusters_mg_cellchat.rds")
