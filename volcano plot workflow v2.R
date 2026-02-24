# want DEGs list and volcano plots for CON M vs CON F, CON M vs COLD M, CON F vs COLD F
# for cycling NSPCs (drawing file DEGs from "DEGs for ORA workflow")

library(Seurat)
library(dplyr)
library(patchwork)
library(clustree)
library(pheatmap)
library(EnhancedVolcano)
library(ggplot2)
library(ggrepel)


# CON M vs. CON F (Cycling NSPCs) ------------------------------------------

# Load each DEG .rds file into the list 
cluster <- "Cycling_NSPCs"

filepath <- paste0("/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/DEG objects/", cluster, "_CON_DEGs_Male_vs_Female_(no_p-value_cutoff).rds")

DEGs <- readRDS(filepath)

# Add gene names
DEGs$gene <- rownames(DEGs)

# Categorize differential expression
DEGs$diffexpressed <- "NO"
DEGs$diffexpressed[DEGs$avg_log2FC > 0 & DEGs$p_val_adj < 0.05] <- "UP"
DEGs$diffexpressed[DEGs$avg_log2FC < 0 & DEGs$p_val_adj < 0.05] <- "DOWN"


# Plot
ggplot(DEGs, aes(x = avg_log2FC, y = -log10(p_val_adj), color = diffexpressed)) +
  geom_point(size = 2) +
  geom_vline(xintercept = c(0, 0), linetype = "dashed", color = "gray") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
  scale_color_manual(values = c("DOWN" = "#00AFBB", "NO" = "gray", "UP" = "darkred"),
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  theme_minimal() +
  labs(title = paste("Volcano Plot -", cluster),
       x = "Log2 Fold Change",
       y = "-log10(Adjusted p-value)",
       color = "Expression") +
  geom_text_repel(data = subset(DEGs, diffexpressed != "NO"),
                  aes(label = gene), size = 3, max.overlaps = 20)

# CHECK -------------------------------------------------------------------

# All upregulated
upregulated <- DEGs %>%
  filter(avg_log2FC > 0, p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC))

write.csv(upregulated,
          paste0("/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/volcano plots2/check/all_upregulated_", cluster, "_CON_M_vs_F.csv"),
          row.names = FALSE)

# All downregulated
downregulated <- DEGs %>%
  filter(avg_log2FC < 0, p_val_adj < 0.05) %>%
  arrange(avg_log2FC) 

write.csv(downregulated,
          paste0("/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/volcano plots2/check/all_downregulated_", cluster, "_CON_M_vs_F.csv"),
          row.names = FALSE)


# CON M vs. COLD M (Cycling NSPCs) ------------------------------------------

# Load each DEG .rds file into the list 
cluster <- "Cycling_NSPCs"

filepath <- paste0("/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/DEG objects/", cluster, "_Male_DEGs_CON_vs_COLD_(no_p-value_cutoff).rds")

DEGs <- readRDS(filepath)

# Add gene names
DEGs$gene <- rownames(DEGs)

# Categorize differential expression
DEGs$diffexpressed <- "NO"
DEGs$diffexpressed[DEGs$avg_log2FC > 0 & DEGs$p_val_adj < 0.05] <- "UP"
DEGs$diffexpressed[DEGs$avg_log2FC < 0 & DEGs$p_val_adj < 0.05] <- "DOWN"


# Plot
ggplot(DEGs, aes(x = avg_log2FC, y = -log10(p_val_adj), color = diffexpressed)) +
  geom_point(size = 2) +
  geom_vline(xintercept = c(0, 0), linetype = "dashed", color = "gray") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
  scale_color_manual(values = c("DOWN" = "#00AFBB", "NO" = "gray", "UP" = "darkred"),
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  theme_minimal() +
  labs(title = paste("Volcano Plot -", cluster),
       x = "Log2 Fold Change",
       y = "-log10(Adjusted p-value)",
       color = "Expression") +
  geom_text_repel(data = subset(DEGs, diffexpressed != "NO"),
                  aes(label = gene), size = 3, max.overlaps = 20)

# CHECK -------------------------------------------------------------------

# All upregulated
upregulated <- DEGs %>%
  filter(avg_log2FC > 0, p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) 

write.csv(upregulated,
          paste0("/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/volcano plots2/check/all_upregulated_", cluster, "_CON_M_vs_COLD_M.csv"),
          row.names = FALSE)

# All downregulated
downregulated <- DEGs %>%
  filter(avg_log2FC < 0, p_val_adj < 0.05) %>%
  arrange(avg_log2FC) 

write.csv(downregulated,
          paste0("/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/volcano plots2/check/all_downregulated_", cluster, "_CON_M_vs_COLD_M.csv"),
          row.names = FALSE)


# CON F vs. COLD F (Cycling NSPCs) ------------------------------------------

# Load each DEG .rds file into the list 
cluster <- "Cycling_NSPCs"

filepath <- paste0("/Users/kodeebao/Documents/Strategy 2 Final Analysis/saved objects/DEG objects/", cluster, "_Female_DEGs_CON_vs_COLD_(no_p-value_cutoff).rds")

DEGs <- readRDS(filepath)

# Add gene names
DEGs$gene <- rownames(DEGs)

# Categorize differential expression
DEGs$diffexpressed <- "NO"
DEGs$diffexpressed[DEGs$avg_log2FC > 0 & DEGs$p_val_adj < 0.05] <- "UP"
DEGs$diffexpressed[DEGs$avg_log2FC < 0 & DEGs$p_val_adj < 0.05] <- "DOWN"


# Plot
ggplot(DEGs, aes(x = avg_log2FC, y = -log10(p_val_adj), color = diffexpressed)) +
  geom_point(size = 2) +
  geom_vline(xintercept = c(0, 0), linetype = "dashed", color = "gray") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
  scale_color_manual(values = c("DOWN" = "#00AFBB", "NO" = "gray", "UP" = "darkred"),
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  theme_minimal() +
  labs(title = paste("Volcano Plot -", cluster),
       x = "Log2 Fold Change",
       y = "-log10(Adjusted p-value)",
       color = "Expression") +
  geom_text_repel(data = subset(DEGs, diffexpressed != "NO"),
                  aes(label = gene), size = 3, max.overlaps = 20)

# CHECK -------------------------------------------------------------------

# All upregulated
upregulated <- DEGs %>%
  filter(avg_log2FC > 0, p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) 

write.csv(upregulated,
          paste0("/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/volcano plots2/check/all_upregulated_", cluster, "_CON_F_vs_COLD_F.csv"),
          row.names = FALSE)

# All downregulated
downregulated <- DEGs %>%
  filter(avg_log2FC < 0, p_val_adj < 0.05) %>%
  arrange(avg_log2FC) 

write.csv(downregulated,
          paste0("/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/volcano plots2/check/all_downregulated_", cluster, "_CON_F_vs_COLD_F.csv"),
          row.names = FALSE)


# for presentation --------------------------------------------------------

# load in the data you want to analyze, CON M vs CON F, CON M vs COLD M, CON F vs COLD F

# CON M vs CON F
# Define genes to label
# got genes from "DEG analysis for thesis"

genes_to_label <- c("Cbfa2t3","Sp9","Meis2","Abracl","Arl4d","Rbp1","Dlx5","Dlx1","Epha4","Gm17750",
                     "Cited2","Cdca7","Pou3f4","Bcl11b","Zeb2","Plxna2","Dlx2","Lmo1","Nrxn3","Foxg1",
                     "Gm4258","H4c4","Ttc9b","Sfta3-ps","Tmem131l","Sox11","Pou3f3","Arx","Lima1","Dclk2",  
                     "Zc3h13","Egr1","Nkx2-1","Ldhb","Scg5","Pantr1","Dlx6","Mex3b","Tcf4","Elavl3",
                     "Pou3f2","Auts2","Ccdc88a","Rrm2","Gdpd1","Tulp4","Myo1b","Ank2","Map4k4","Enc1","Chd3"  )

ggplot(DEGs, aes(x = avg_log2FC, y = -log10(p_val_adj), color = diffexpressed)) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
  scale_color_manual(values = c("DOWN" = "#00AFBB", "NO" = "gray", "UP" = "darkred"),
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  theme_minimal() +
  labs(title = paste("Volcano Plot -", cluster),
       x = "Log2 Fold Change",
       y = "-log10(Adjusted p-value)",
       color = "Expression") +
  geom_text_repel(
    data = subset(DEGs, gene %in% genes_to_label),
    aes(label = gene),
    size = 3,
    max.overlaps = 20,
    color = "darkred",
    segment.color = "darkred")


# CON F vs COLD F
# Define genes to label
genes_to_label <- c("Cbfa2t3","Sp9","Meis2","Abracl","Arl4d","Rbp1","Dlx5","Dlx1","Epha4","Gm17750",
                    "Cited2","Cdca7","Pou3f4","Bcl11b","Zeb2","Plxna2","Dlx2","Lmo1","Nrxn3","Foxg1",
                    "Gm4258","H4c4","Ttc9b","Sfta3-ps","Tmem131l","Sox11","Pou3f3","Arx","Lima1","Dclk2",  
                    "Zc3h13","Egr1","Nkx2-1","Ldhb","Scg5","Pantr1","Dlx6","Mex3b","Tcf4","Elavl3",
                    "Pou3f2","Auts2","Ccdc88a","Rrm2","Gdpd1","Tulp4","Myo1b","Ank2","Map4k4","Enc1","Chd3"  )

ggplot(DEGs, aes(x = avg_log2FC, y = -log10(p_val_adj), color = diffexpressed)) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
  scale_color_manual(values = c("DOWN" = "#00AFBB", "NO" = "gray", "UP" = "darkred"),
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  theme_minimal() +
  labs(title = paste("Volcano Plot -", cluster),
       x = "Log2 Fold Change",
       y = "-log10(Adjusted p-value)",
       color = "Expression") +
  geom_text_repel(
    data = subset(DEGs, gene %in% genes_to_label),
    aes(label = gene),
    size = 3,
    max.overlaps = 20,
    color = "darkred",
    segment.color = "darkred")
