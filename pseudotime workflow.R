# SF ran RNA velocity on CON vs COLD (M and F), CON M vs COLD M, CON F vs COLD F 

library(dplyr)
library(clustree)
library(Signac)
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)


# this is for clean combined object (CON and COLD, Male and Female) --------

# load in objects
clean_combined_object <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/clean_combined_object_CON&COLD_strat2.rds')

# Pseudotime Workflow
# Step 1: Convert Seurat object to Monocle 3 CDS
clean_combined_object_cds <- as.cell_data_set(clean_combined_object)

# Step 2: Transfer PCA and UMAP embeddings from Seurat to Monocle
reducedDims(clean_combined_object_cds) <- SimpleList(PCA = Embeddings(clean_combined_object, "pca"))
reducedDims(clean_combined_object_cds)$UMAP <- Embeddings(clean_combined_object, "umap")

# Step 3: Run clustering (Monocle needs cluster information) 
# tried clustering within Monocle (overriding Seurat clusters, but did not work, only made 1 partition/cluster)
colData(clean_combined_object_cds)$seurat_clusters <- Idents(clean_combined_object)
clean_combined_object_cds <- cluster_cells(clean_combined_object_cds, reduction_method = "UMAP", cluster_cells = FALSE)

# Step 4: Learn the trajectory (pseudotime)
# prune_graph = FALSE tells Monocle not to trim "weak" connections when learning the trajectory. We have diverse or distant clusters. Want more complete graph coverage. Default graph doesn’t reach all clusters. We do know the cycling RGC (cluster 2) is biologically relevant so should be included 

clean_combined_object_cds <- learn_graph(
  clean_combined_object_cds,
  use_partition = FALSE,  # ← this forces Monocle to build ONE connected graph
  learn_graph_control = list(prune_graph = FALSE))  # disables Monocle's pruning

# Step 4b: choosing root nodes manually (will discuss with JR)
clean_combined_object_cds <- order_cells(clean_combined_object_cds)

# Step 5: Visualize cells along pseudotime on UMAP
plot_cells(clean_combined_object_cds, color_cells_by = "pseudotime", 
           label_groups_by_cluster = TRUE, show_trajectory_graph = F)

# JR wants these plotted without the lines

