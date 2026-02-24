# looking into all the genes overlapping in CON M and COLD F

# Load necessary library
library(dplyr)


# CON M and COLD F --------------------------------------------------------

# Load upregulated genes
con_m_vs_f <- read.csv("/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/volcano plots2/check/all_upregulated_Cycling_NSPCs_CON_M_vs_F.csv")
con_f_vs_cold_f <- read.csv("/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/volcano plots2/check/all_upregulated_Cycling_NSPCs_CON_F_vs_COLD_F.csv")

# Extract gene names
genes_con_m <- con_m_vs_f$gene
genes_cold_f <- con_f_vs_cold_f$gene

# Count total genes
length(genes_con_m)       # 123
length(genes_cold_f)      # 109

# Find intersection
shared_genes <- intersect(genes_con_m, genes_cold_f)
length(shared_genes)      # 51

# View shared genes
shared_genes

# Convert to a data frame without sorting
shared_genes_df <- data.frame(gene = shared_genes)

# Export to CSV in your desired folder
write.csv(shared_genes_df, 
          "/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/volcano plots2/check/shared_upregulated_genes_CON_M_vs_F_and_CON_F_vs_COLD_F.csv",
          row.names = FALSE)


# CON M and COLD M --------------------------------------------------------

# Load upregulated genes
con_m_vs_f <- read.csv("/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/volcano plots2/check/all_upregulated_Cycling_NSPCs_CON_M_vs_F.csv")
con_m_vs_cold_m <- read.csv("/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/volcano plots2/check/all_upregulated_Cycling_NSPCs_CON_M_vs_COLD_M.csv")

# Extract gene names
genes_con_m <- con_m_vs_f$gene
genes_cold_m <- con_m_vs_cold_m$gene

# Count total genes
length(genes_con_m)       # 123
length(genes_cold_m)      # 20

# Find intersection
shared_genes <- intersect(genes_con_m, genes_cold_m)
length(shared_genes)      # 0

# View shared genes
shared_genes


# COLD M and COLD F -------------------------------------------------------

# Load upregulated genes
con_m_vs_cold_m <- read.csv("/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/volcano plots2/check/all_upregulated_Cycling_NSPCs_CON_M_vs_COLD_M.csv")
con_f_vs_cold_f <- read.csv("/Users/kodeebao/Documents/Strategy 2 Final Analysis/Analysis on Vim+ and subsetted cells/volcano plots2/check/all_upregulated_Cycling_NSPCs_CON_F_vs_COLD_F.csv")

# Extract gene names
genes_cold_m <- con_m_vs_cold_m$gene
genes_cold_f <- con_f_vs_cold_f$gene

# Count total genes
length(genes_cold_m)       # 20
length(genes_cold_f)      # 109

# Find intersection
shared_genes <- intersect(genes_cold_m, genes_cold_f)
length(shared_genes)      # 8

# View shared genes
shared_genes
# ribosomal proteins and uncharacterized genes
