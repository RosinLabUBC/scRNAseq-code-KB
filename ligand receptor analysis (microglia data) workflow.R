# subsetting microglia CON data and microglia COLD data (following workflow for NSCs) to prepare for ligand receptor analysis

library(dplyr)
library(Seurat)
library(patchwork)
library(clustree)

# quality control and merging objects -------------------------------------

# load data
microglia.CON <- Read10X(data.dir = "/Users/kodeebao/Documents/Microglia scRNAseq Data_CON")
microglia.COLD <- Read10X(data.dir = "/Users/kodeebao/Documents/Microglia scRNAseq Data_COLD")

# create Seurat objects and filter
Microglia_CON <- CreateSeuratObject(counts = microglia.CON, project = "Microglia_CON", min.cells = 3, min.features = 200)
Microglia_CON[["percent.mt"]] <- PercentageFeatureSet(Microglia_CON, pattern = "^mt-")
Microglia_CON <- subset(Microglia_CON, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)					
VlnPlot(Microglia_CON, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)					

Microglia_COLD <- CreateSeuratObject(counts = microglia.COLD, project = "Microglia_COLD", min.cells = 3, min.features = 200)
Microglia_COLD[["percent.mt"]] <- PercentageFeatureSet(Microglia_COLD, pattern = "^mt-")
Microglia_COLD <- subset(Microglia_COLD, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)					
VlnPlot(Microglia_COLD, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)					

# Number of cells after filtering
ncol(Microglia_CON)   # 1143
ncol(Microglia_COLD)  # 1424

#Tagging each object with a treatment label 
Microglia_CON$treatment <- "CON"
Microglia_COLD$treatment <- "COLD"

# process each object (Normalize and scale RNA assay)
Microglia_CON <- NormalizeData(Microglia_CON)
Microglia_CON <- FindVariableFeatures(Microglia_CON)
Microglia_CON <- ScaleData(Microglia_CON, features = row.names(Microglia_CON))

Microglia_COLD <- NormalizeData(Microglia_COLD)
Microglia_COLD <- FindVariableFeatures(Microglia_COLD)
Microglia_COLD <- ScaleData(Microglia_COLD, features = row.names(Microglia_COLD))

# Merge COLD_NSC and CON_NSC
merged_object_mg <- merge(
  x = Microglia_CON,                          # reference object
  y = Microglia_COLD,                         # object to merge
  add.cell.ids = c("CON", "COLD"),     # prefixes to distinguish cells
  project = "Merged_MG")

# check number of cells total
ncol(merged_object_mg)   # 2567

saveRDS(merged_object_mg, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/microglia data/CON_COLD_merged_microglia.rds")

# analysis without integration --------------------------------------------------

# load in merged_object
merged_object_mg <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/microglia data/CON_COLD_merged_microglia.rds')

# Normalize merged object
merged_object_mg <- NormalizeData(merged_object_mg)
merged_object_mg <- FindVariableFeatures(merged_object_mg)
merged_object_mg <- ScaleData(merged_object_mg, features = row.names(merged_object_mg))

DefaultAssay(merged_object_mg) <- "RNA"
merged_object_mg <- RunPCA(merged_object_mg)
ElbowPlot(merged_object_mg, ndims = 40, reduction = "pca")		

saveRDS(merged_object_mg, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/microglia data/CON_COLD_merged_microglia_normalized.rds")


# clustering merged and normalized microglia data -------------------------
# load in merged_object
merged_object_mg <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/microglia data/CON_COLD_merged_microglia_normalized.rds')

# Confirm number of cells
ncol(merged_object_mg)
# 2567 cells

# From JR CELL Paper
# just know that they used 0.5 res in methods
# supplemental showed they had 7 clusters total

# AL clustered 11 with merged object

# going with 30 dims, 0.5 res (checked with JR)
# Clusters to keep: 0,1,2,3,5,7
# Macrophage = 9
# Proliferative / NSC = 6
# Pericyte / RBC = 8
# Unknown = 4


# 30 dims, 0.5 res to cluster ----------------------------------------------
# keeping clusters 0,1,2,3,5,7

merged_object_mg <- FindNeighbors(merged_object_mg, dims = 1:30)	
merged_object_mg <- FindClusters(merged_object_mg, resolution = 0.5)		
merged_object_mg <- RunUMAP(merged_object_mg, dims = 1:30)		

DimPlot(merged_object_mg, reduction = "umap")
DimPlot(merged_object_mg, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"), split.by = "orig.ident")

# Specify the clusters to keep
keep_clusters <- c(0, 1, 2, 3, 5, 7)

# Subset the object to keep only those clusters
merged_object_mg_clean <- subset(merged_object_mg, idents = keep_clusters)

# Save original cluster identities if needed
merged_object_mg_clean$original_cluster <- Idents(merged_object_mg_clean)

# Assign all cells the identity "microglia"
Idents(merged_object_mg_clean) <- "microglia"

# Check the original sample identities (should still be CON and COLD)
table(merged_object_mg_clean$treatment)
#    COLD        CON 
#    1200        824      (Total 2024)

DefaultAssay(merged_object_mg_clean) <- "RNA"
merged_object_mg_clean <- RunPCA(merged_object_mg_clean)
ElbowPlot(merged_object_mg_clean, ndims = 40, reduction = "pca")		

saveRDS(merged_object_mg_clean, "/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/microglia data/merged_object_mg_clean.rds")

# removing male and female double positive and non-expressing cells -------
# following same workflow as NSCs
# going to just use "Xist" for females and other cells not expressing are male cells (confirmed with JR)

# load in object
merged_object_mg_clean <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/microglia data/merged_object_mg_clean.rds')

# define markers and expression threshold
female_genes <- c("Xist")
expr_threshold <- 0.1

# Use normalized RNA assay for expression values
DefaultAssay(merged_object_mg_clean) <- "RNA"

# get data from each matrix (CON and COLD)
expr_con <- GetAssayData(merged_object_mg_clean, assay = "RNA", layer = "data.Microglia_CON")
expr_cold <- GetAssayData(merged_object_mg_clean, assay = "RNA", layer = "data.Microglia_COLD")

# Identify CON female and male cells
female_mg_con <- colnames(expr_con)[colSums(expr_con[female_genes, , drop = FALSE] > expr_threshold) > 0]
male_mg_con <- setdiff(colnames(expr_con), female_mg_con)

# Identify COLD female and male cells
female_mg_cold <- colnames(expr_cold)[colSums(expr_cold[female_genes, , drop = FALSE] > expr_threshold) > 0]
male_mg_cold <- setdiff(colnames(expr_cold), female_mg_cold)

#confirm number of cells
length(female_mg_con) # 269
length(male_mg_con)   # 555

length(female_mg_cold) # 1110
length(male_mg_cold)   # 90

# checked and there are no cells overlapping 

# Store in metadata
# Create a named vector to store sex labels
sex_labels <- rep(NA, ncol(merged_object_mg_clean))
names(sex_labels) <- colnames(merged_object_mg_clean)

# Assign labels
sex_labels[female_mg_con] <- "Female"
sex_labels[male_mg_con]   <- "Male"
sex_labels[female_mg_cold] <- "Female"
sex_labels[male_mg_cold]   <- "Male"

# Check for NAs (e.g., if some cells are not part of either group)
table(is.na(sex_labels))  # FALSE (2024), no cells were left unclassified

# Add to Seurat metadata
merged_object_mg_clean$sex <- sex_labels

saveRDS(merged_object_mg_clean, file = "/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/microglia data/merged_object_mg_clean_labelled.rds") 


# creating separate suerat objects ----------------------------------------

# load in clean, merged, labelled microglia object
merged_object_mg_clean_labelled <- readRDS("/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/microglia data/merged_object_mg_clean_labelled.rds")

# splitting "merged_object_mg_clean_labelled" into by sex and CON/COLD

# Female + CON
female_con_mg_obj <- subset(merged_object_mg_clean_labelled, subset = sex == "Female" & treatment == "CON")

# Female + COLD
female_cold_mg_obj <- subset(merged_object_mg_clean_labelled, subset = sex == "Female" & treatment == "COLD")

# Male + CON
male_con_mg_obj <- subset(merged_object_mg_clean_labelled, subset = sex == "Male" & treatment == "CON")

# Male + COLD
male_cold_mg_obj <- subset(merged_object_mg_clean_labelled, subset = sex == "Male" & treatment == "COLD")

# check
sapply(list(female_con_mg_obj, female_cold_mg_obj, male_con_mg_obj, male_cold_mg_obj), ncol) # 269, 1110, 555, 90

# save
saveRDS(female_con_mg_obj, file = "/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/microglia data/female_con_mg_obj.rds")
saveRDS(female_cold_mg_obj, file = "/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/microglia data/female_cold_mg_obj.rds")
saveRDS(male_con_mg_obj, file = "/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/microglia data/male_con_mg_obj.rds")
saveRDS(male_cold_mg_obj, file = "/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/microglia data/male_cold_mg_obj.rds")

# load
female_con_mg_obj <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/microglia data/female_con_mg_obj.rds')
female_cold_mg_obj <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/microglia data/female_cold_mg_obj.rds')
male_con_mg_obj <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/microglia data/male_con_mg_obj.rds')
male_cold_mg_obj <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/microglia data/male_cold_mg_obj.rds')
