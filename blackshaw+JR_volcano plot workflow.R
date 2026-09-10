# drawing files from blackshaw & kb analysis 1 folder

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
filepath <- paste0("/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/code for suppl 4/blackshaw + kb (analysis 1)/analysis on Vim+ subset/CON_DEGs_Male_vs_Female_CyclingNSPCs_no_pcutoff.rds")

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
  labs(title = paste("Volcano Plot - CON M vs CON F (Blackshaw Dataset)"),
       x = "Log2 Fold Change",
       y = "-log10(Adjusted p-value)",
       color = "Expression") +
  geom_text_repel(data = subset(DEGs, diffexpressed != "NO"),
                  aes(label = gene), size = 3, max.overlaps = 50)

# CHECK -------------------------------------------------------------------

# All upregulated
upregulated <- DEGs %>%
  filter(avg_log2FC > 0, p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC))

write.csv(upregulated,
          paste0("/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/code for suppl 4/suppl 4 images/DEGs check/upregulated_CON_M_vs_F.csv"),
          row.names = FALSE)

# All downregulated
downregulated <- DEGs %>%
  filter(avg_log2FC < 0, p_val_adj < 0.05) %>%
  arrange(avg_log2FC) 

write.csv(downregulated,
          paste0("/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/code for suppl 4/suppl 4 images/DEGs check/downregulated_CON_M_vs_F.csv"),
          row.names = FALSE)

