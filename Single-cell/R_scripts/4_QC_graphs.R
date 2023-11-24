# Libraries
library(Seurat)
library(tidyverse)

options(ggrepel.max.overlaps = Inf)
set.seed(1)

# set working directory
wd <- "~/B_Memory_master/Master_Diploma_Private/Output/Single-cell/"
setwd(wd)

# Create directory for output if it doesn't exist
if (!dir.exists("Graphs_png")) {
  dir.create("Graphs_png")
}

load("3_QC.RData")

# Function to save plot to png
plot_scatter <- function(name, plot) {
  png(paste0("Graphs_png/4_", name, ".png"), 
      width = 10, 
      height = 7, 
      res = 450, 
      units = "in")
  print(plot)
  dev.off()
}

# Function to save plot to png at lower resolution
plot_scatter_small <- function(name, plot) {
  png(paste0("Graphs_png/4_", name, "_small.png"), 
      width = 6, 
      height = 4, 
      res = 450, 
      units = "in")
  print(plot)
  dev.off()
}

# Part 1: Visualizing QC Metrics with Violin Plots
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb", "percent_largest_gene")
for (i in feats) {
  graph <- VlnPlot(LC, group.by = "patient", features = i, pt.size = 0.01, raster=FALSE) + NoLegend()
  plot_scatter(paste0("QC_vln_", i), graph)
  plot_scatter_small(paste0("QC_vln_", i), graph)
}

# Part 2: Scatterplots for QC Metrics
qc.metrics <- as_tibble(LC[[]], rownames="Cell.Barcode")

plot_scatter("QC_dotplot_count_feature_doublets", ggplot(qc.metrics, aes(nCount_RNA,nFeature_RNA,colour=doublet_finder)) + geom_point(size = 0.001) + scale_color_manual(values=c("blue","red")) + ggtitle("Example of plotting QC metrics") + geom_hline(yintercept = c(200, 4000)) + theme_bw())
plot_scatter_small("QC_dotplot_count_feature_doublets", ggplot(qc.metrics, aes(nCount_RNA,nFeature_RNA,colour=doublet_finder)) + geom_point(size = 0.001) + scale_color_manual(values=c("blue","red")) + ggtitle("Example of plotting QC metrics") + geom_hline(yintercept = c(200, 4000)) + theme_bw())


# Part 3: Complexity Calculations and Visualizations
qc.metrics <- qc.metrics %>%
  mutate(complexity=log10(nFeature_RNA) / log10(nCount_RNA))

complexity.lm <- lm(log10(qc.metrics$nFeature_RNA)~log10(qc.metrics$nCount_RNA)) #  The lm() function is used to fit linear models to data frames

qc.metrics <- qc.metrics %>%
  mutate(complexity_diff = log10(nFeature_RNA) - ((log10(qc.metrics$nCount_RNA)*complexity.lm$coefficients[2])+complexity.lm$coefficients[1]))

plot_scatter("QC_area_complexity", ggplot(qc.metrics, aes(x=complexity_diff)) + geom_density(fill="grey") + theme_bw())
plot_scatter_small("QC_area_complexity", ggplot(qc.metrics, aes(x=complexity_diff)) + geom_density(fill="grey") + theme_bw())

plot_scatter("QC_dotplot_count_feature_complexity_log_scale", ggplot(qc.metrics %>% mutate(complexity_diff = replace(complexity_diff, complexity_diff < -0.1, -0.1)), aes(x=log10(nCount_RNA), y=log10(nFeature_RNA), colour=complexity_diff)) + geom_point(size=0.5) + geom_abline(slope=complexity.lm$coefficients[2], intercept = complexity.lm$coefficients[1]) + scale_colour_gradient2(low="blue2",mid="grey",high="red2") + theme_bw())
plot_scatter_small("QC_dotplot_count_feature_complexity_log_scale", ggplot(qc.metrics %>% mutate(complexity_diff = replace(complexity_diff, complexity_diff < -0.1, -0.1)), aes(x=log10(nCount_RNA), y=log10(nFeature_RNA), colour=complexity_diff)) + geom_point(size=0.5) + geom_abline(slope=complexity.lm$coefficients[2], intercept = complexity.lm$coefficients[1]) + scale_colour_gradient2(low="blue2",mid="grey",high="red2") + theme_bw())

# Part 4: Analyzing and Visualizing Largest Genes
largest_gene_list <- qc.metrics %>% group_by(largest_gene) %>% summarise(count = n()) %>% arrange(desc(count))
largest_genes_to_plot <- largest_gene_list %>% filter(count > 5000) %>% pull(largest_gene)

plot_scatter("QC_dotplot_count_feature_largest_gene_log_scale", ggplot(qc.metrics %>% filter(largest_gene %in% largest_genes_to_plot) %>% mutate(largest_gene = factor(largest_gene, levels = largest_genes_to_plot)) %>% arrange(largest_gene), aes(x=log10(nCount_RNA), y=log10(nFeature_RNA), colour=largest_gene)) + geom_point(size=0.5) + scale_colour_manual(values=RColorBrewer::brewer.pal(9, "Set1")) + theme_bw())
plot_scatter_small("QC_dotplot_count_feature_largest_gene_log_scale", ggplot(qc.metrics %>% filter(largest_gene %in% largest_genes_to_plot) %>% mutate(largest_gene = factor(largest_gene, levels = largest_genes_to_plot)) %>% arrange(largest_gene), aes(x=log10(nCount_RNA), y=log10(nFeature_RNA), colour=largest_gene)) + geom_point(size=0.5) + scale_colour_manual(values=RColorBrewer::brewer.pal(9, "Set1")) + theme_bw())

# Part 5: Further QC Metrics Visualizations
plot_scatter("QC_dotplot_mito_ribo_feature", ggplot(qc.metrics, aes(percent_ribo,percent_mito,colour=nFeature_RNA)) + geom_point(size = 0.001) + scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) + ggtitle("Example of plotting QC metrics") + theme_bw())
plot_scatter_small("QC_dotplot_mito_ribo_feature", ggplot(qc.metrics, aes(percent_ribo,percent_mito,colour=nFeature_RNA)) + geom_point(size = 0.001) + scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) + ggtitle("Example of plotting QC metrics") + theme_bw())

# Part 6: Additional QC and Saving the Workspace
LC <- subset(LC, nFeature_RNA < 4000 & nFeature_RNA > 200 & percent_mito < 20 & percent_ribo > 5)

# rm all accept LC
rm(list = ls()[-match("LC", ls())])

# save workspace
save.image("4_Filetered_QC.RData")
