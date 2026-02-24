# Using All Upregualted and Downregulated DEGs (and DEGs have parameters of log2FC > 0.5, < -0.5, adj p-value < 0.05) -- same workflow as 'GO ORA Analysis.R'

library(clusterProfiler)
library(org.Mm.eg.db)
library(Seurat)
library(dplyr)
library(tibble)
library(enrichplot)
library(ggplot2)

# GO ORA for all CON M vs CON F --------
# Load upregulated and downregulated genes
all_CON_DEGs_top <- read.csv("/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Gad1+ : Gad2+ cells/Analysis on Gad+ clean cells/DEGs/DEGs for ORA/all_upregulated_CON_M_vs_F.csv")

all_CON_DEGs_down <- read.csv("/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Gad1+ : Gad2+ cells/Analysis on Gad+ clean cells/DEGs/DEGs for ORA/all_downregulated_CON_M_vs_F.csv")

# load in object for gene universe (use subsetted object for universe, SF said makes sense when trying to explore)

# load in object
merged_object_gad_pos_clean <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Gad1+ : Gad2+ cells/objects/merged_object_gad_pos_clean.rds')

# Set default assay to RNA
DefaultAssay(merged_object_gad_pos_clean) <- "RNA"

# Normalize and scale RNA assay
merged_object_gad_pos_clean <- NormalizeData(merged_object_gad_pos_clean)
merged_object_gad_pos_clean <- FindVariableFeatures(merged_object_gad_pos_clean)
merged_object_gad_pos_clean <- ScaleData(merged_object_gad_pos_clean, features = row.names(merged_object_gad_pos_clean))

# Join Layers
merged_object_gad_pos_clean <- JoinLayers(merged_object_gad_pos_clean, layers = "data")

# Subset to just the CON cells
con_obj <- subset(merged_object_gad_pos_clean, subset = treatment == "CON")

# Extract gene symbols
# checked gene symbols are indeed in the X column
up_genes <- all_CON_DEGs_top$X
down_genes <- all_CON_DEGs_down$X

# Use RNA assay with normalized expression (log-normalized)
rna_data <- GetAssayData(con_obj, assay = "RNA", layer = "data")

# Calculate percent of cells expressing each gene
pct_expr <- rowSums(rna_data > 0) / ncol(rna_data)

# Filter genes expressed in >10% of cells (min.pct = 0.1)
gene_universe <- names(pct_expr[pct_expr > 0.1])

# Combine all genes for conversion
all_genes_to_convert <- unique(c(up_genes, down_genes, gene_universe))

# Convert to ENSEMBL IDs
gene_conversion <- bitr(all_genes_to_convert,
                        fromType = "SYMBOL",
                        toType = "ENSEMBL",
                        OrgDb = org.Mm.eg.db)

# Match converted gene lists
up_ensembl <- gene_conversion %>%
  filter(SYMBOL %in% up_genes) %>%
  pull(ENSEMBL)

down_ensembl <- gene_conversion %>%
  filter(SYMBOL %in% down_genes) %>%
  pull(ENSEMBL)

universe_ensembl <- gene_conversion %>%
  filter(SYMBOL %in% gene_universe) %>%
  pull(ENSEMBL)

# GO Enrichment for Upregulated Genes
ego_up <- enrichGO(gene          = up_ensembl,
                   universe      = universe_ensembl,
                   OrgDb         = org.Mm.eg.db,
                   keyType       = "ENSEMBL",
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   qvalueCutoff  = 0.05, #SF said to remove p-value, q-value is enough
                   readable      = TRUE)

head(ego_up)

# GO Enrichment for Downregulated Genes
ego_down <- enrichGO(gene          = down_ensembl,
                     universe      = universe_ensembl,
                     OrgDb         = org.Mm.eg.db,
                     keyType       = "ENSEMBL",
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     qvalueCutoff  = 0.05, #SF said to remove p-value, q-value is enough
                     readable      = TRUE)

head(ego_down)

# Save results to your specified folder
write.csv(as.data.frame(ego_up), "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Gad1+ : Gad2+ cells/Analysis on Gad+ clean cells/GO ORA/CON M vs CON F/ORA_BP_Upregulated_CON M_vs_CON F_minpct_01.csv")
write.csv(as.data.frame(ego_down), "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Gad1+ : Gad2+ cells/Analysis on Gad+ clean cells/GO ORA/CON M vs CON F/ORA_BP_Downregulated_CON M_vs_CON F_minpct_01.csv")


# GO ORA for all CON F vs COLD F --------
# Load upregulated and downregulated genes
all_Female_DEGs_top <- read.csv("/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Gad1+ : Gad2+ cells/Analysis on Gad+ clean cells/DEGs/DEGs for ORA/all_Female_upregulated_CON_vs_COLD.csv")

all_Female_DEGs_down <- read.csv("/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Gad1+ : Gad2+ cells/Analysis on Gad+ clean cells/DEGs/DEGs for ORA/all_Female_downregulated_CON_vs_COLD.csv")

# load in object for gene universe (use subsetted object for universe, SF said makes sense when trying to explore)

# load in object
merged_object_gad_pos_clean <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Gad1+ : Gad2+ cells/objects/merged_object_gad_pos_clean.rds')

# Set default assay to RNA
DefaultAssay(merged_object_gad_pos_clean) <- "RNA"

# Normalize and scale RNA assay
merged_object_gad_pos_clean <- NormalizeData(merged_object_gad_pos_clean)
merged_object_gad_pos_clean <- FindVariableFeatures(merged_object_gad_pos_clean)
merged_object_gad_pos_clean <- ScaleData(merged_object_gad_pos_clean, features = row.names(merged_object_gad_pos_clean))

# Join Layers
merged_object_gad_pos_clean <- JoinLayers(merged_object_gad_pos_clean, layers = "data")

# Subset to just the Female cells
female_obj <- subset(merged_object_gad_pos_clean, subset = sex == "Female")

# Extract gene symbols
# checked gene symbols are indeed in the X column
up_genes <- all_Female_DEGs_top$X
down_genes <- all_Female_DEGs_down$X

# Use RNA assay with normalized expression (log-normalized)
rna_data <- GetAssayData(female_obj, assay = "RNA", layer = "data")

# Calculate percent of cells expressing each gene
pct_expr <- rowSums(rna_data > 0) / ncol(rna_data)

# Filter genes expressed in >10% of cells (min.pct = 0.1)
gene_universe <- names(pct_expr[pct_expr > 0.1])

# Combine all genes for conversion
all_genes_to_convert <- unique(c(up_genes, down_genes, gene_universe))

# Convert to ENSEMBL IDs
gene_conversion <- bitr(all_genes_to_convert,
                        fromType = "SYMBOL",
                        toType = "ENSEMBL",
                        OrgDb = org.Mm.eg.db)

# Match converted gene lists
up_ensembl <- gene_conversion %>%
  filter(SYMBOL %in% up_genes) %>%
  pull(ENSEMBL)

down_ensembl <- gene_conversion %>%
  filter(SYMBOL %in% down_genes) %>%
  pull(ENSEMBL)

universe_ensembl <- gene_conversion %>%
  filter(SYMBOL %in% gene_universe) %>%
  pull(ENSEMBL)

# GO Enrichment for Upregulated Genes
ego_up <- enrichGO(gene          = up_ensembl,
                   universe      = universe_ensembl,
                   OrgDb         = org.Mm.eg.db,
                   keyType       = "ENSEMBL",
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   qvalueCutoff  = 0.05, #SF said to remove p-value, q-value is enough
                   readable      = TRUE)

head(ego_up)

# GO Enrichment for Downregulated Genes
ego_down <- enrichGO(gene          = down_ensembl,
                     universe      = universe_ensembl,
                     OrgDb         = org.Mm.eg.db,
                     keyType       = "ENSEMBL",
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     qvalueCutoff  = 0.05, #SF said to remove p-value, q-value is enough
                     readable      = TRUE)

head(ego_down)

# Save results to your specified folder
write.csv(as.data.frame(ego_up), "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Gad1+ : Gad2+ cells/Analysis on Gad+ clean cells/GO ORA/CON F vs COLD F/ORA_BP_Upregulated_CON F_vs_COLD F_minpct_01.csv")
write.csv(as.data.frame(ego_down), "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Gad1+ : Gad2+ cells/Analysis on Gad+ clean cells/GO ORA/CON F vs COLD F/ORA_BP_Downregulated_CON F_vs_COLD F_minpct_01.csv")

# no downregulated GO terms for COLD F vs CON F 


# GO ORA for all CON M vs COLD M --------
# Load upregulated and downregulated genes
all_Male_DEGs_top <- read.csv("/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Gad1+ : Gad2+ cells/Analysis on Gad+ clean cells/DEGs/DEGs for ORA/all_Male_upregulated_CON_vs_COLD.csv")

all_Male_DEGs_down <- read.csv("/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Gad1+ : Gad2+ cells/Analysis on Gad+ clean cells/DEGs/DEGs for ORA/all_Male_downregulated_CON_vs_COLD.csv")

# load in object for gene universe (use subsetted object for universe, SF said makes sense when trying to explore)

# load in object
merged_object_gad_pos_clean <- readRDS('/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Gad1+ : Gad2+ cells/objects/merged_object_gad_pos_clean.rds')

# Set default assay to RNA
DefaultAssay(merged_object_gad_pos_clean) <- "RNA"

# Normalize and scale RNA assay
merged_object_gad_pos_clean <- NormalizeData(merged_object_gad_pos_clean)
merged_object_gad_pos_clean <- FindVariableFeatures(merged_object_gad_pos_clean)
merged_object_gad_pos_clean <- ScaleData(merged_object_gad_pos_clean, features = row.names(merged_object_gad_pos_clean))

# Join Layers
merged_object_gad_pos_clean <- JoinLayers(merged_object_gad_pos_clean, layers = "data")

# Subset to just the Male cells
male_obj <- subset(merged_object_gad_pos_clean, subset = sex == "Male")

# Extract gene symbols
# checked gene symbols are indeed in the X column
up_genes <- all_Male_DEGs_top$X
down_genes <- all_Male_DEGs_down$X

# Use RNA assay with normalized expression (log-normalized)
rna_data <- GetAssayData(male_obj, assay = "RNA", layer = "data")

# Calculate percent of cells expressing each gene
pct_expr <- rowSums(rna_data > 0) / ncol(rna_data)

# Filter genes expressed in >10% of cells (min.pct = 0.1)
gene_universe <- names(pct_expr[pct_expr > 0.1])

# Combine all genes for conversion
all_genes_to_convert <- unique(c(up_genes, down_genes, gene_universe))

# Convert to ENSEMBL IDs
gene_conversion <- bitr(all_genes_to_convert,
                        fromType = "SYMBOL",
                        toType = "ENSEMBL",
                        OrgDb = org.Mm.eg.db)

# Match converted gene lists
up_ensembl <- gene_conversion %>%
  filter(SYMBOL %in% up_genes) %>%
  pull(ENSEMBL)

down_ensembl <- gene_conversion %>%
  filter(SYMBOL %in% down_genes) %>%
  pull(ENSEMBL)

universe_ensembl <- gene_conversion %>%
  filter(SYMBOL %in% gene_universe) %>%
  pull(ENSEMBL)

# GO Enrichment for Upregulated Genes
ego_up <- enrichGO(gene          = up_ensembl,
                   universe      = universe_ensembl,
                   OrgDb         = org.Mm.eg.db,
                   keyType       = "ENSEMBL",
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   qvalueCutoff  = 0.05, #SF said to remove p-value, q-value is enough
                   readable      = TRUE)

head(ego_up)

# GO Enrichment for Downregulated Genes
ego_down <- enrichGO(gene          = down_ensembl,
                     universe      = universe_ensembl,
                     OrgDb         = org.Mm.eg.db,
                     keyType       = "ENSEMBL",
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     qvalueCutoff  = 0.05, #SF said to remove p-value, q-value is enough
                     readable      = TRUE)

head(ego_down)

# Save results to your specified folder
write.csv(as.data.frame(ego_up), "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Gad1+ : Gad2+ cells/Analysis on Gad+ clean cells/GO ORA/CON M vs COLD M/ORA_BP_Upregulated_CON M_vs_COLD M_minpct_01.csv")
write.csv(as.data.frame(ego_down), "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Gad1+ : Gad2+ cells/Analysis on Gad+ clean cells/GO ORA/CON M vs COLD M/ORA_BP_Downregulated_CON M_vs_COLD M_minpct_01.csv")
