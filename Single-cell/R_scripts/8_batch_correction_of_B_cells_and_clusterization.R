library(Seurat)
library(tidyverse)
library(patchwork)
library(RColorBrewer)
library(ggplot2)
library(reshape2)

options(ggrepel.max.overlaps = Inf)
set.seed(1)

# set working directory
wd <- "~/B_Memory_master/Master_Diploma_Private/Output/Single-cell/"
setwd(wd)

# Create directory for output if it doesn't exist
if (!dir.exists("Graphs_png")) {
  dir.create("Graphs_png")
}

# Helper function to get color palette
getPalette <- function(data_col, n) {
  len <- length(unique(data_col))
  return(colorRampPalette(brewer.pal(min(n, 9), "Set1"))(len))
}

# Helper function for feature plotting with custom colors
featurePlotCustom <- function(data, feature) {
  FeaturePlot(data, features = feature, reduction = "umap", pt.size = 0.05, raster=FALSE) + 
    scale_colour_gradientn(colours = c("grey80", "#3366CC", "#0000FF", "#0000CC", "#000099"))
}

# visualize the percent of patients in clusters
plot_integrated_clusters = function (srat) {
  
  count_table <- table(srat@meta.data$seurat_clusters, srat@meta.data$patient) # patient, lib, cancer_type
  count_mtx   <- as.data.frame.matrix(count_table)
  count_mtx$cluster <- rownames(count_mtx)
  melt_mtx    <- melt(count_mtx)
  melt_mtx$cluster <- as.factor(melt_mtx$cluster)
  
  cluster_size   <- aggregate(value ~ cluster, data = melt_mtx, FUN = sum)
  
  sorted_labels <- paste(sort(as.integer(levels(cluster_size$cluster)),decreasing = T))
  cluster_size$cluster <- factor(cluster_size$cluster,levels = sorted_labels)
  melt_mtx$cluster <- factor(melt_mtx$cluster,levels = sorted_labels)
  
  colnames(melt_mtx)[2] <- "dataset"
  
  # vis
  p1 <- ggplot(cluster_size, aes(y= cluster,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
    theme_bw() + scale_x_log10() + xlab("Cells per cluster, log10 scale") + ylab("")
  
  set.seed(5)
  getPalette_perc = sample(c(colorRampPalette(brewer.pal(12, 'Paired'))(length(unique(melt_mtx$dataset)) - 2), "grey80", "grey20"))
  p2 <- ggplot(melt_mtx,aes(x=cluster,y=value,fill=dataset)) + 
    geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
    ylab("Fraction of cells in each dataset") + xlab("Cluster number") +
    scale_fill_manual(values = getPalette_perc) + theme(legend.position="top")
  
  p2 + p1 + plot_layout(widths = c(4,1))
  
}

# Function to save plot to png
saveToPng <- function(plot, path_suffix) {
  png(paste0("Graphs_png/8_", path_suffix, ".png"), 
      width = 10, 
      height = 7, 
      res = 450, 
      units = "in")
  print(plot)
  dev.off()
}

# Function to save plot to png at lower resolution
saveToPng_small <- function(plot, path_suffix) {
  png(paste0("Graphs_png/8_", path_suffix, "_small.png"), 
      width = 6, 
      height = 4, 
      res = 450, 
      units = "in")
  print(plot)
  dev.off()
}

# Load data and modify object
load("7_Clustered_B_cells.RData") 
LC_B <- LC_B %>% 
  subset(nFeature_RNA < 3000 & percent_mito < 15 & percent_ribo > 10) %>% 
  subset(tissue_type == "tum") %>% 
  subset(doublet_finder == "Singlet")

LC_B@meta.data <- LC_B@meta.data[c(1:4, 6:7, 9:13)]

# select var feat
data(cc.genes)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
var.genes_pattern <- "TR[ABGD][CV]|RP[SL]|IG[HKL][VDJ]|IGHM|IGHD|IGHE|IGHA[1-2]|IGHG[1-4]|IGKC|IGLC_B[1-7]|AC233755.1|IGLL|JCHAIN|MT-|MALAT1|MTRNR|XIST|IGHGP|HSP|DNAJ|JUN|DUSP|IER|FOS|BTG|ATF3|NR4A1"
var.genes <- grep(pattern = var.genes_pattern, x = rownames(LC_B), value = TRUE)
excl_genes <- c(s.genes, g2m.genes, var.genes)

# Preprocess data for integration
Idents(LC_B) <- LC_B$patient
# LC_B.list <- SplitObject(LC_B, split.by = "patient")
# LC_B.list
LC_B <- subset(LC_B, idents = c(406, 458, 464, 460, 569, 630, 593), invert = TRUE) # https://bioinformatics.stackexchange.com/questions/15720/error-in-findintegrationanchors-seurat-package
Idents(LC_B) <- LC_B$orig.ident

LC_B.list <- SplitObject(LC_B, split.by = "patient")
rm(LC_B)

# normalize and identify variable features for each dataset independently
LC_B.list <- lapply(X = LC_B.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 1000)
  VariableFeatures(x) <- VariableFeatures(x)[-which(VariableFeatures(x) %in% excl_genes)]
  x
})

# Select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = LC_B.list)
immune.anchors <- FindIntegrationAnchors(object.list = LC_B.list, anchor.features = features, normalization.method = "LogNormalize")
rm(LC_B.list)

# Integrate Data
immune.combined <- IntegrateData(anchorset = immune.anchors, normalization.method = "LogNormalize", k.weight = 15)

# Without integration
DefaultAssay(immune.combined) <- "RNA"

immune.combined <- NormalizeData(immune.combined, verbose = F)
immune.combined <- FindVariableFeatures(immune.combined, selection.method = "vst", nfeatures = 500, verbose = F)

VariableFeatures(immune.combined) <- VariableFeatures(immune.combined)[-which(VariableFeatures(immune.combined) %in% excl_genes)]

immune.combined <- ScaleData(immune.combined, verbose = F)
immune.combined <- RunPCA(immune.combined, npcs = 20, verbose = F)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20, verbose = F)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20) 
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

saveToPng(DimPlot(immune.combined, reduction = "umap", label = TRUE, group.by = "patient", cols = getPalette(immune.combined$patient, 9), raster=FALSE) + plot_annotation(title = "without integration"), "DimPlot_patient_without_integration")
saveToPng_small(DimPlot(immune.combined, reduction = "umap", label = TRUE, group.by = "patient", cols = getPalette(immune.combined$patient, 9), raster=FALSE) + plot_annotation(title = "without integration"), "DimPlot_patient_without_integration")

# Saving BarPlot to png
saveToPng(plot_integrated_clusters(immune.combined) , "Barplot_patients_without_integration")
saveToPng_small(plot_integrated_clusters(immune.combined) , "Barplot_patients_without_integration")

# Switch to integrated assay. The variable features of this assay are automatically calculated

DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = F) 
immune.combined <- RunPCA(immune.combined, verbose = F)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.6)

# Saving DimPlots to png
dim_plots_integrated <- list(
  DimPlot(immune.combined, reduction = "umap", group.by = "patient", cols = getPalette(immune.combined$patient, 9), raster=FALSE),
  DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE, cols = getPalette(immune.combined$seurat_clusters, 9), raster=FALSE) + plot_annotation(title = "with integration")
)
dim_paths_integrated <- c("patient", "seurat_clusters")

for(i in seq_along(dim_paths_integrated)) {
  saveToPng(dim_plots_integrated[[i]], paste0("DimPlot_", dim_paths_integrated[i], "_integrated"))
  saveToPng_small(dim_plots_integrated[[i]], paste0("DimPlot_", dim_paths_integrated[i], "_integrated"))
}

# Saving FeaturePlot to png
saveToPng(featurePlotCustom(immune.combined, "percent_mito"), "FeaturePlot_percent_mito_integrated")
saveToPng_small(featurePlotCustom(immune.combined, "percent_mito"), "FeaturePlot_percent_mito_integrated")

# For performing differential expression after integration, we switch back 
# to the original data

DefaultAssay(immune.combined) <- "RNA"

# Saving VlnPlots to png
vln_paths_integrated <- c("percent_mito", "nFeature_RNA", "percent_ribo")
lapply(vln_paths_integrated, function(feat) {
  plot <- VlnPlot(immune.combined, feat, pt.size = 0.01) +
    stat_summary(fun = mean, geom='point', size = 10, colour = "black", shape = 95)
  saveToPng(plot, paste0("VlnPlot_", feat, "_integrated"))
  saveToPng_small(plot, paste0("VlnPlot_", feat, "_integrated"))
})

# Saving FeaturePlots to png
featureNames_integrated <- c("FCRL4", "CD3E", "CD38", "TCL1A", "IGHG1", "IGHD", "CD27", "PDCD1")
lapply(featureNames_integrated, function(feat) {
  plot <- featurePlotCustom(immune.combined, feat)
  saveToPng(plot, paste0("FeaturePlot_", feat, "_integrated"))
  saveToPng_small(plot, paste0("FeaturePlot_", feat, "_integrated"))
})

# Saving VlnPlots to png
vln_paths_integrated_counts <- c("FCRL4", "CD3E", "CD38", "TCL1A", "IGHG1", "IGHD", "CD27", "PDCD1")
lapply(vln_paths_integrated_counts, function(feat) {
  plot <- VlnPlot(immune.combined, feat) +
    stat_summary(fun = mean, geom='point', size = 5, colour = "black", shape = 95)
  saveToPng(plot, paste0("VlnPlot_", feat, "_integrated"))
  saveToPng_small(plot, paste0("VlnPlot_", feat, "_integrated"))
})

# Saving BarPlot to png
saveToPng(plot_integrated_clusters(immune.combined) , "Barplot_patients_integrated")
saveToPng_small(plot_integrated_clusters(immune.combined) , "Barplot_patients_integrated")

# search for DE genes among all clusters
all_markers <- FindAllMarkers(immune.combined) %>% filter(p_val_adj < 0.001)
write.table(all_markers[, c(ncol(all_markers), 1:(ncol(all_markers)-1))], "Tables/8_all_DE_genes.tsv", row.names = T, col.names = T, sep = "\t")

# rm all accept immune.combined
rm(list = ls()[-match("immune.combined", ls())])

# Save workspace
save.image("8_Corrected_B_cells.RData")
