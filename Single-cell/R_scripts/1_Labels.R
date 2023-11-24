# Libraries
library(Seurat)
library(DoubletFinder)
library(readr)  # for read_csv
library(tidyverse) # contains stringr

options(ggrepel.max.overlaps = Inf)
set.seed(1)

# set working directory
wd <- "~/B_Memory_master/Master_Diploma_Private/Output/Single-cell/"
setwd(wd)

# Function to download a file if it doesn't exist
download_if_not_exists <- function(url, destfile) {
  if (!file.exists(destfile)) {
    message(paste("Downloading", destfile))
    mode_option <- if (Sys.info()["sysname"] == "Windows") "wb" else ""
    download.file(url = url, destfile = destfile, mode = mode_option)
  }
}

data_dir <- file.path(wd, "Data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# Download and load the main data
data_file <- file.path(data_dir, "lung_ldm.rd")
data_url <- "https://www.dropbox.com/s/vjbide8ro5iwrfh/lung_ldm.rd?dl=1"
download_if_not_exists(data_url, data_file)

load(data_file)

counts <- lung_ldm[["dataset"]][["umitab"]]
rm(lung_ldm)

LC <- CreateSeuratObject(counts = counts, project = "LC", min.cells = 3)
rm(counts)

# Download and read metadata
metadata_file <- file.path(data_dir, "table_s1_sample_table.csv")
metadata_url <- "https://raw.githubusercontent.com/effiken/Leader_et_al/4a884161d50ed768963603c8a0aea38ea4c9299b/input_tables/table_s1_sample_table.csv"
download_if_not_exists(metadata_url, metadata_file)

metadata <- read_csv(metadata_file)

# Filter metadata according to your needs
metadata <- metadata %>%
  filter(!sample_ID %in% c(115, 116, 343, 344, 480, 481),
         !patient_ID %in% str_extract(patient_ID, "(Lambrechts_\\d|zilionis_\\d)"))

# Prepare look-up tables for fast access
lookup_patient <- setNames(metadata$patient_ID, metadata$sample_ID)
lookup_tissue <- setNames(metadata$tissue, metadata$sample_ID)
lookup_disease <- setNames(metadata$disease, metadata$patient_ID)
lookup_prep <- setNames(metadata$prep, metadata$patient_ID)

# 1. patient ID
LC$patient <- lookup_patient[as.character(LC$orig.ident)]

# 2. cancer/not
LC$tissue_type <- lookup_tissue[as.character(LC$orig.ident)]

# 3. LUAD/LUSC
LC$cancer_type <- lookup_disease[as.character(LC$patient)]

# 4. Library preparation method
LC$lib <- sapply(LC$patient, function(p) {
  if (lookup_prep[p] == "digest/dead cell") "dead"
  else if (lookup_prep[p] == "sort") "sort"
  else "beads"
})

# rm all accept LC
rm(list = ls()[-match("LC", ls())])

# save workspace
save.image("1_Labeled.RData")
