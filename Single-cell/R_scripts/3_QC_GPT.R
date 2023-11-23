# Libraries
library(Seurat)

options(ggrepel.max.overlaps = Inf)
set.seed(1)

# set working directory
wd <- "~/B_Memory_master/Master_Diploma_Private/Output/Single-cell/"
setwd(wd)

load("2_Labeled_doublets.RData")

LC <- PercentageFeatureSet(LC, pattern = "^MT-", col.name = "percent_mito")
LC <- PercentageFeatureSet(LC, pattern = "^RP[SL]", col.name = "percent_ribo")
LC <- PercentageFeatureSet(LC, pattern = "^HB[^(P)]", col.name = "percent_hb")

# calculate % of the largest gene
data.temp <- LC[rownames(LC) != "MALAT1",]
data.temp$largest_count <- apply(data.temp@assays$RNA@counts, 2, max)
data.temp$largest_index <- apply(data.temp@assays$RNA@counts, 2, which.max)
data.temp$largest_gene <- rownames(data.temp)[data.temp$largest_index]
data.temp$percent_largest_gene <- 100 * data.temp$largest_count / data.temp$nCount_RNA
LC$largest_gene <- data.temp$largest_gene
LC$percent_largest_gene <- data.temp$percent_largest_gene

# rm all accept LC
rm(list = ls()[-match("LC", ls())])

# save workspace
save.image("3_QC.RData")
