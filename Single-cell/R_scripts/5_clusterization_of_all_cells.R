# Libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

options(ggrepel.max.overlaps = Inf)
set.seed(1)

# set working directory
wd <- "~/B_Memory_master/Master_Diploma_Private/Output/Single-cell/"
setwd(wd)

load("4_Filetered_QC.RData")

# Utility function to save plots to PDF
save_plot <- function(plot_obj, filename) {
  pdf(paste0("Graphs/5_", filename, ".pdf"), height=7, width=14)
  print(plot_obj)
  dev.off()
}

# Analysis
LC <- FindVariableFeatures(LC, selection.method = "vst", nfeatures = 2000)

# Show and save most variable features
a <- LabelPoints(plot = VariableFeaturePlot(LC), points = head(VariableFeatures(LC), 30), repel = FALSE)
save_plot(a, "Dotplot_variable_features")

# Filter genes using dplyr
data(cc.genes)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

var.genes_pattern <- "TR[ABGD][CV]|RP[SL]|IG[HKL][VDJ]|IGHM|IGHD|IGHE|IGHA[1-2]|IGHG[1-4]|IGKC|IGLC[1-7]|AC233755.1|IGLL|JCHAIN|MT-|MALAT1|MTRNR|XIST|IGHGP|HSP|DNAJ|JUN|DUSP|IER|FOS|BTG|ATF3|NR4A1"
var.genes <- grep(pattern = var.genes_pattern, x = rownames(LC), value = TRUE)

VariableFeatures(LC) %>% 
  setdiff(var.genes) %>% 
  setdiff(s.genes) %>% 
  setdiff(g2m.genes) -> VariableFeatures(LC)

# Show and save most variable features after filtration
b <- LabelPoints(plot = VariableFeaturePlot(LC), points = head(VariableFeatures(LC), 30), repel = FALSE)
save_plot(b, "Dotplot_variable_features_after_filtration")

LC <- ScaleData(LC, features = rownames(LC))
LC <- RunPCA(LC, features = VariableFeatures(object = LC))
LC <- FindNeighbors(LC, dims = 1:40)
LC <- FindClusters(LC, resolution = 2.5)
LC <- RunUMAP(LC, dims = 1:40)

# Save DimPlot
getPalette_clust = colorRampPalette(brewer.pal(9, "Set1"))(length(unique(LC$seurat_clusters)))
c <- DimPlot(LC, reduction = "umap", label = TRUE, cols = getPalette_clust) + NoLegend()
save_plot(c, "Dim_plot")

# rm all accept LC
rm(list = ls()[-match("LC", ls())])

# Save workspace
save.image("5_Clustered_all_cells.RData")
