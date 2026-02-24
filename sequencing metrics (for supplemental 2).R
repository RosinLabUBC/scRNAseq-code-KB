# workflow from "workflow to cluster and subset.R"

library(Seurat)
library(dplyr)
library(patchwork)
library(clustree)

# clustering merged object that is SCTransformed --------------------------
# merged object contains data sets COLD M & F and CON M & F, they have been dieted then SCTransformed again

# load in merged_object
merged_object <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/merged_object_SCTransformed.rds')

# Confirm number of cells
ncol(merged_object)
# 22643 cells

# Number of cells by condition
table(merged_object$orig.ident)
# COLD_NSC  CON_NSC 
# 12564    10079 

# Mean UMI counts per cell by condition
merged_object@meta.data %>%
  group_by(orig.ident) %>%
  summarise(mean_UMI_per_cell = mean(nCount_RNA))
# COLD_NSC                 8204
# CON_NSC                 10042

# Median genes per cell by condition
merged_object@meta.data %>%
  group_by(orig.ident) %>%
  summarise(median_genes_per_cell = median(nFeature_RNA))
# COLD_NSC                    2987
# CON_NSC                     3577


# umap showing cells
DefaultAssay(merged_object) <- "SCT"
merged_object <- RunPCA(merged_object)
ElbowPlot(merged_object, ndims = 50, reduction = "pca")

merged_object <- FindNeighbors(merged_object, dims = 1:35)	
merged_object <- FindClusters(merged_object, resolution = 0.15)		
merged_object <- RunUMAP(merged_object, dims = 1:35)
DimPlot(merged_object, reduction = "umap")
DimPlot(merged_object, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"), split.by = "orig.ident")
DimPlot(merged_object, reduction = "umap", group.by = "orig.ident", label = F)

FeaturePlot(merged_object, features = "Vim")

# load in Vim+ subsetted cells --------------------------------------------

merged_object_vim_pos <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/merged_object_vim_pos.rds')

# Confirm number of cells
ncol(merged_object_vim_pos)
# 13624 cells

# using 35 dims, 0.1 res to cluster (found on clustree on separate workflow "Testing Parameters to cluster COLD and CON datasets v2, Vim+ only.R" and discussed with JR)
# 0.1 res looked to have more cells in cluster 1 (fewer red cells in cluster 1)
merged_object_vim_pos <- FindNeighbors(merged_object_vim_pos, dims = 1:35)	
merged_object_vim_pos <- FindClusters(merged_object_vim_pos, resolution = 0.1)		
merged_object_vim_pos <- RunUMAP(merged_object_vim_pos, dims = 1:35)
DimPlot(merged_object_vim_pos, reduction = "umap")
DimPlot(merged_object_vim_pos, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"), split.by = "orig.ident")

FeaturePlot(merged_object_vim_pos, features = "Ascl1")

# load in Vim+ cells that have been subsetted (cluster 1 or Ascl1+) -----------------------------------------

merged_object_vim_pos_subset <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/merged_object_vim_pos_subset.rds')

# testing parameters on cluster 1 + Ascl1+ subset on separate workflow 

# using 30 dims, 0.2 res -- confirmed with JR -----------------------------
merged_object_vim_pos_subset <- FindNeighbors(merged_object_vim_pos_subset, dims = 1:30)	
merged_object_vim_pos_subset <- FindClusters(merged_object_vim_pos_subset, resolution = 0.2)		
merged_object_vim_pos_subset <- RunUMAP(merged_object_vim_pos_subset, dims = 1:30)
DimPlot(merged_object_vim_pos_subset, reduction = "umap")
DimPlot(merged_object_vim_pos_subset, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"), split.by = "orig.ident")


# clean combined object analysis ------------------------------------------

clean_combined_object <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/clean_combined_object_CON&COLD_strat2.rds')

DimPlot(clean_combined_object, reduction = "umap")
DimPlot(clean_combined_object, reduction = "umap", group.by = "orig.ident", split.by = "orig.ident")
DimPlot(clean_combined_object, reduction = "umap", group.by = "cluster_names", cols = cluster_colors, split.by = "orig.ident")

table(clean_combined_object$sex)
# Female   Male 
# 2364    2121 

table(clean_combined_object$treatment)
# COLD  CON 
# 2392 2093 

table(clean_combined_object$cluster_names)
# Cycling NSPCs               Glial Precursors   Neuronal Precursors (Ascl1+) 
#  756                           1137                           1839 
# Neuronal Precursors (Neurog2+)                      NPCs-Glia                  NPCs-Neuronal 
#             77                                      408                            268 

# Median & Mean UMI counts per cell by condition
clean_combined_object@meta.data %>%
  group_by(orig.ident) %>%
  summarise(median_UMI_per_cell = median(nCount_RNA))

# COLD_NSC               7310
# CON_NSC               9623

clean_combined_object@meta.data %>%
  group_by(orig.ident) %>%
  summarise(mean_UMI_per_cell = mean(nCount_RNA))

# COLD_NSC               8749
# CON_NSC               1093

# Median & Mean genes per cell by condition
clean_combined_object@meta.data %>%
  group_by(orig.ident) %>%
  summarise(median_genes_per_cell = median(nFeature_RNA))
# COLD_NSC                    3079
# CON_NSC                     3683

clean_combined_object@meta.data %>%
  group_by(orig.ident) %>%
  summarise(mean_genes_per_cell = mean(nFeature_RNA))
# COLD_NSC                    3183
# CON_NSC                     3662

# vln plot for UMI counts split by orig.ident

VlnPlot( clean_combined_object,
         features = "nCount_RNA",
         group.by = "orig.ident") +
         labs(
           x = "Treatment",
           y = "UMI Count")


# vln plot for gene counts split by orig.ident

VlnPlot( clean_combined_object,
         features = "nFeature_RNA",
         group.by = "orig.ident") +
         labs(
          x = "Treatment",
          y = "Gene Count")

# vln plot for UMI counts split by cluster_names

VlnPlot( clean_combined_object,
         features = "nCount_RNA",
         group.by = "cluster_names",
         cols = cluster_colors) +
         labs(
          x = "Treatment",
          y = "UMI Count")

# vln plot for gene counts split by cluster_names

VlnPlot( clean_combined_object,
         features = "nFeature_RNA",
         group.by = "cluster_names",
         cols = cluster_colors) +
         labs(
          x = "Treatment",
          y = "Gene Count")

# mito genes vs UMI graph

df <- clean_combined_object@meta.data

ggplot(df, aes(x = nCount_RNA, y = percent.mt, color = orig.ident)) +
  geom_point(alpha = 0.6, size = 0.8) +
  theme_classic() +
  scale_x_log10() +  # log scale for UMI counts (helps spread out the data)
  labs(
    x = "UMI Counts per Cell",
    y = "% Mitochondrial Genes",
    color = "orig.ident")


# gene count vs UMI graph

ggplot(df, aes(x = nCount_RNA, y = nFeature_RNA, color = orig.ident)) +
  geom_point(alpha = 0.6, size = 0.8) +
  theme_classic() +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = "UMI Counts per Cell",
    y = "Genes Counts per Cell",
    color = "orig.ident")

# mitochondrial genes graph
VlnPlot(
  clean_combined_object,
  features = "percent.mt",
  group.by = "cluster_names",
  split.by = "orig.ident",
  split.plot = T) +
  theme_classic() +
  labs(
    x = "Cluster",
    y = "% Mitochondrial Genes") + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1))

# apoptosis graph
apoptosis_genes <- c(
  "Aifm1", "Akt1", "Akt2", "Akt3", "Apaf1", "Atm", "Bad", "Bax", "Bcl2", 
  "Bcl2l1", "Bid", "Birc2", "Birc3", "Capn1", "Capn2", "Casp10", "Casp3", 
  "Casp6", "Casp7", "Casp8", "Casp9", "Cflar", "Chp1", "Chp2", "Chuk", 
  "Csf2rb", "Cycs", "Dffa", "Dffb", "Endod1", "Endog", "Exog", "Fadd", 
  "Fas", "Faslg", "Ikbkb", "Ikbkg", "Il1a", "Il1b", "Il1r1", "Il1rap", 
  "Il3", "Il3ra", "Irak1", "Irak2", "Irak3", "Irak4", "Map3k14", "Myd88", 
  "Nfkb1", "Nfkbia", "Ngf", "Ntrk1", "Pik3ca", "Pik3cb", "Pik3cd", "Pik3cg", 
  "Pik3r1", "Pik3r2", "Pik3r3", "Pik3r5", "Ppp3ca", "Ppp3cb", "Ppp3cc", 
  "Ppp3r1", "Ppp3r2", "Prkaca", "Prkacb", "Prkacg", "Prkar1a", "Prkar1b", 
  "Prkar2a", "Prkar2b", "Prkx", "Rela", "Ripk1", "Tnf", "Tnfrsf10a", 
  "Tnfrsf10b", "Tnfrsf10c", "Tnfrsf10d", "Tnfrsf1a", "Tnfsf10", "Tp53", 
  "Tradd", "Traf2", "Xiap")

clean_combined_object <- AddModuleScore(
  object = clean_combined_object,
  features = list(apoptosis_genes),
  name = "Apoptosis")

VlnPlot(
  clean_combined_object,
  features = "Apoptosis1",
  group.by = "cluster_names",
  split.by = "orig.ident",
  split.plot = T) +
  theme_classic() +
  labs(
    x = "Cluster",
    y = "Apoptosis (KEGG M8492) Module Score") + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1))


# % mito genes per cluster and apoptosis score with t-test ----------------

# pull metadata
meta <- clean_combined_object@meta.data

# make sure treatment is a factor (CON vs COLD)
meta$treatment <- recode(
  meta$orig.ident,
  "CON_NSC"  = "CON",
  "COLD_NSC" = "COLD") %>% factor(levels = c("CON", "COLD"))

# check it
table(meta$treatment)

# Function: run t-test per cluster
run_cluster_ttest <- function(df, feature) {
  df %>%
    group_by(cluster_names) %>%
    summarise(
      p_value = tryCatch(
        pairwise.t.test(
          x = get(feature),
          g = treatment,
          p.adjust.method = "none"   # TRUE pairwise t-test
        )$p.value[1, 1],
        error = function(e) NA),
      mean_CON  = mean(get(feature)[treatment == "CON"], na.rm = TRUE),
      mean_COLD = mean(get(feature)[treatment == "COLD"], na.rm = TRUE),
      n_CON     = sum(treatment == "CON"),
      n_COLD    = sum(treatment == "COLD")
    ) %>%
    mutate(
      diff_COLD_vs_CON = mean_COLD - mean_CON,
      direction = case_when(
        diff_COLD_vs_CON > 0 ~ "Higher in COLD",
        diff_COLD_vs_CON < 0 ~ "Lower in COLD",
        TRUE ~ "No difference"),
      padj = p.adjust(p_value, method = "fdr"),
      interpretation = ifelse(padj < 0.05, direction, "Not significant"))}  

# Run t-tests

# 1. Mitochondrial %
mito_results <- run_cluster_ttest(meta, "percent.mt")
print("Mitochondrial % pair-wise t-test results")
print(mito_results)

# 2. Apoptosis module score
apop_results <- run_cluster_ttest(meta, "Apoptosis1")
print("Apoptosis module score pair-wise t-test results")
print(apop_results)

# Save results (optional)
# mito t-test results
write.csv(
  mito_results,
  "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/supplemental 3 figures/pairwise_t_test_mito_results.csv",
  row.names = FALSE)

# apoptosis t-test results
write.csv(
  apop_results,
  "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/supplemental 3 figures/pairwise_t_test_apoptosis_results.csv",
  row.names = FALSE)
