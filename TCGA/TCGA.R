# Load libraries
library(TCGAbiolinks)
library(survminer)
library(survival)
library(SummarizedExperiment)
library(tidyverse)
library(DESeq2)

# set working directory
wd <- "~/B_Memory_master/Master_Diploma_Private/Output/Bulk/"
setwd(wd)

# Create directory for output if it doesn't exist
if (!dir.exists("Graphs_png")) {
  dir.create("Graphs_png")
}

# Download clinical data for lung adenocarcinoma (LUAD)
clinical_luad <- GDCquery_clinic("TCGA-LUAD")

# Display relevant survival-related variables
clinical_luad[,c("vital_status", "days_to_last_follow_up", "days_to_death")]

# Check for NA values in vital_status
is.na(clinical_luad[,"vital_status"])

# Remove rows with NA values in vital_status
clinical_luad <- clinical_luad[complete.cases(clinical_luad$vital_status), ]

# Verify that NA rows have been removed
print(is.na(clinical_luad[, "vital_status"]))

# Check distribution of vital status
table(clinical_luad$vital_status)

# Encode deceased status
clinical_luad$deceased <- ifelse(clinical_luad$vital_status == "Alive", FALSE, TRUE)

# Check the distribution of deceased status
table(clinical_luad$deceased)

# Create an "overall survival" variable
clinical_luad$overall_survival <- ifelse(clinical_luad$deceased,
                                         clinical_luad$days_to_death,
                                         clinical_luad$days_to_last_follow_up)


# Gene expression data retrieval ------------------------------------------


# Build a query to get gene expression data for the entire cohort
query_luad_all = GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  data.type = "Gene Expression Quantification",
  sample.type = "Primary Tumor",
  access = "open")

output_luad <- getResults(query_luad_all)

# # Get primary tissue sample barcodes
tumor <- output_luad$cases
tumor

# Build a query to get gene expression data from 50 primary tumors
query_luad <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  data.type = "Gene Expression Quantification",
  sample.type = c("Primary Tumor", "Solid Tissue Normal"),
  access = "open",
  barcode = tumor)

# Download data
GDCdownload(query_luad)

# Prepare gene expression data
tcga_luad_data <- GDCprepare(query_luad, summarizedExperiment = TRUE)
luad_matrix <- assay(tcga_luad_data, "unstranded")

# Extract gene and sample metadata
gene_metadata <- as.data.frame(rowData(tcga_luad_data))
coldata <- as.data.frame(colData(tcga_luad_data))


# VST transform counts for survival analysis ------------------------------


# Setting up countData object   
dds <- DESeqDataSetFromMatrix(countData = luad_matrix,
                              colData = coldata,
                              design = ~ 1)

# Remove genes with a sum total of less than 10 reads across all samples
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# VST transformation
vsd <- vst(dds, blind=FALSE)
luad_matrix_vst <- assay(vsd)


# Extract data for FCRL4 gene ---------------------------------------------


luad_genes_FCRL4 <- luad_matrix_vst %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'case_id', value = 'counts', -gene_id) %>% 
  left_join(., gene_metadata, by = "gene_id") %>% 
  filter(gene_name == "FCRL4")

# Check for zero counts
any(luad_genes_FCRL4$counts == 0)

# Get median value
median_value_FCRL4 <- median(luad_genes_FCRL4$counts)

# Categorize cases based on gene expression
luad_genes_FCRL4$strata <- ifelse(luad_genes_FCRL4$counts >= median_value_FCRL4, "HIGH", "LOW")

# Add clinical information to luad_genes_FCRL4
luad_genes_FCRL4$case_id <- gsub('-01.*', '', luad_genes_FCRL4$case_id)
luad_genes_FCRL4 <- merge(luad_genes_FCRL4, clinical_luad, by.x = 'case_id', by.y = 'submitter_id')


# Fit survival curve for FCRL4 --------------------------------------------


fit_FCRL4 <- survfit(Surv(overall_survival, deceased) ~ strata, data = luad_genes_FCRL4)

# Display survival plot without risk table
surv_plot <- ggsurvplot(fit_FCRL4,
                        data = luad_genes_FCRL4,
                        pval = T,
                        risk.table = F) 

ggsave(filename = "Graphs_png/surv_plot.png", plot = surv_plot$plot, device = "png", height=10, width=10)
ggsave(filename = "Graphs_png/surv_plot_small.png", plot = surv_plot$plot, device = "png", height=5, width=5)

# Perform log-rank test
fit2_FCRL4 <- survdiff(Surv(overall_survival, deceased) ~ strata, 
                       data = luad_genes_FCRL4)
fit2_FCRL4


# Extract data for FCRL4 gene normalized to the MS4A1 ---------------------


luad_genes_FCRL4_norm <- luad_matrix_vst %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'case_id', value = 'counts', -gene_id) %>% 
  left_join(., gene_metadata, by = "gene_id") %>% 
  filter(gene_name == "FCRL4" | gene_name == "MS4A1")

# Check for zero counts
any(luad_genes_FCRL4_norm$counts == 0)

luad_genes_FCRL4_norm_1 <- luad_genes_FCRL4_norm[, c("counts", "gene_name", "case_id")]


# Normalization step

# Check whether every value goes twice
if (nrow(luad_genes_FCRL4_norm_1) %% 2 == 0) {
  if (all(table(luad_genes_FCRL4_norm_1$case_id) == 2)) {
    print("The lines come in pairs.")
  }
}

# Filter the dataframe for FCRL4 and MS4A1 rows
fcrl4_data <- subset(luad_genes_FCRL4_norm_1, gene_name == "FCRL4")
ms4a1_data <- subset(luad_genes_FCRL4_norm_1, gene_name == "MS4A1")

# Merge the data based on case_id
merged_data <- merge(fcrl4_data, ms4a1_data, by = "case_id", suffixes = c("_fcrl4", "_ms4a1"))

# Create a new column 'counts' by dividing FCRL4 counts by MS4A1 counts
merged_data$counts <- merged_data$counts_fcrl4 / merged_data$counts_ms4a1

# Select only the relevant columns 'case_id' and 'counts'
result_table_FCRL4_norm <- merged_data[, c("case_id", "counts")]

# Get median value
median_value_FCRL4_norm <- median(result_table_FCRL4_norm$counts)

# Categorize cases based on gene expression
result_table_FCRL4_norm$strata <- ifelse(result_table_FCRL4_norm$counts >= median_value_FCRL4_norm, "HIGH", "LOW")

# Add clinical information to result_table_FCRL4_norm
result_table_FCRL4_norm$case_id <- gsub('-01.*', '', result_table_FCRL4_norm$case_id)
result_table_FCRL4_norm <- merge(result_table_FCRL4_norm, clinical_luad, by.x = 'case_id', by.y = 'submitter_id')


# Fit survival curve for normalizaed FCRL4 --------------------------------


fit_FCRL4_norm <- survfit(Surv(overall_survival, deceased) ~ strata, data = result_table_FCRL4_norm)

# Display survival plot without risk table
surv_plot_norm <- ggsurvplot(fit_FCRL4_norm,
                        data = result_table_FCRL4_norm,
                        pval = T,
                        risk.table = F) 

ggsave(filename = "Graphs_png/surv_plot_norm.png", plot = surv_plot_norm$plot, device = "png", height=10, width=10)
ggsave(filename = "Graphs_png/surv_plot_norm_small.png", plot = surv_plot_norm$plot, device = "png", height=5, width=5)


# Perform log-rank test
fit2_FCRL4_norm <- survdiff(Surv(overall_survival, deceased) ~ strata, 
                            data = result_table_FCRL4_norm)
fit2_FCRL4_norm

# Save workspace
save.image("TCGA.RData")