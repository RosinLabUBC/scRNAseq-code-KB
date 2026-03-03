# looking into all the genes overlapping in CON M and COLD F for Glutamatergic cells

# Load necessary library
library(dplyr)


# CON M and COLD F --------------------------------------------------------

# Load upregulated genes
con_m_vs_f <- read.csv("/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis of Glutamatergic Cells/Analysis on VGlut+ clean cells/volcano plots/check/all_upregulated_CON_M_vs_F.csv")
con_f_vs_cold_f <- read.csv("/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis of Glutamatergic Cells/Analysis on VGlut+ clean cells/volcano plots/check/all_upregulated_CON_F_vs_COLD_F.csv")

# Extract gene names
genes_con_m <- con_m_vs_f$gene
genes_cold_f <- con_f_vs_cold_f$gene

# Count total genes
length(genes_con_m)       # 67
length(genes_cold_f)      # 157

# Find intersection
shared_genes <- intersect(genes_con_m, genes_cold_f)
length(shared_genes)      # 13

# View shared genes
shared_genes

# Convert to a data frame without sorting
shared_genes_df <- data.frame(gene = shared_genes)

# Export to CSV in your desired folder
write.csv(shared_genes_df, 
          "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis of Glutamatergic Cells/Analysis on VGlut+ clean cells/volcano plots/check/shared_upregulated_genes_CON_M_vs_F_and_CON_F_vs_COLD_F.csv",
          row.names = FALSE)


