library(Seurat)
library(velocyto.R)
library(SeuratWrappers)


# for thesis --------------------------------------------------------------

HY <- readRDS("clean_combined_CON_COLD_strat2.rds")
HY <- RunVelocity(object=HY, spliced="SCTspliced", unspliced="SCTunspliced", ambiguous="SCTambiguous", deltaT=1, kCells=25, fit.quantile=0.02)
ident.colors <- (scales::hue_pal())(n = length(x = levels(x = HY)))
names(x = ident.colors) <- levels(x = HY)
cell.colors <- ident.colors[Idents(object = HY)]
names(x = cell.colors) <- colnames(x = HY)
pdf("clean_combined_CON_COLD_strat2.pdf")
show.velocity.on.embedding.cor(emb = Embeddings(object=HY, reduction="umap"), vel = Tool(object=HY, slot = "RunVelocity"), n=200, scale="sqrt", cell.colors = ac(x=cell.colors, alpha=0.5), cex=0.8, arrow.scale=3, show.grid.flow=TRUE, min.grid.cell.mass=0.5, grid.n=40, arrow.lwd=1, do.par=FALSE, cell.border.alpha=0.1)
dev.off()
