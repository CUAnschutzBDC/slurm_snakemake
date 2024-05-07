library(Seurat)
library(tidyverse)
library(cowplot)
library(here)
library(scAnalysisR)
library(muscat)
library(scran)
library(DESeq2)
library(UpSetR)

source(here("src/scripts/muscat_plotting_functions.R"))

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

#args <- commandArgs(trailingOnly = TRUE)

#sample <- args[[1]]
#sample <- gsub("__.*", "", sample)
sample <- "merged"

#sample_info <- args[[4]]
sample_info <- here("files/sample_info.tsv")

#results_dir <- args[[2]]
results_dir <- here("results")

#sample_metadata <- args[[3]]
sample_metadata <- here("files/Deidentified_donor_metadata.xlsx")
sample_metadata <- openxlsx::readWorkbook(sample_metadata, detectDates = TRUE,
                                          sheet = "New_metadata") %>%
  dplyr::select(Sample.Name, millions.of.cells.frozen, dplyr::contains("num_"))

sample_metadata$Sample.Name <- gsub("JH313-15_JBM",
                                    "JH313-15_JB",
                                    sample_metadata$Sample.Name)

sample_info <- read.table(sample_info, fill = TRUE, header = TRUE)

samples_use <- sample_info[sample_info$sample != sample, ]$sample

sample_info <- sample_info[sample_info$sample == sample,]

HTO <- sample_info$HTO
ADT <- sample_info$ADT
hash_ident <- sample_info$hash_ident


ADT_pca <- sample_info$adt_pca
run_adt_umap <- sample_info$adt_umap

RNA_pcs <- sample_info$PCs
resolution <- sample_info$resolution
batch_correction <- sample_info$batch_correction



if(normalization_method == "SCT"){
  SCT <- TRUE
  seurat_assay <- "SCT"
} else {
  SCT <- FALSE
  seurat_assay <- "RNA"
}

# Set directories
save_dir <- file.path(results_dir, "R_analysis", sample)

# Read in data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed_no_doublet.rds"))

seurat_data$age_bin <- ifelse(seurat_data$Age.at.Collection..years. <= 10,
                              "less_than_10",
                              ifelse(seurat_data$Age.at.Collection..years. <= 20,
                                     "10_to_20", "more_than_20"))

sample_metadata <- sample_metadata %>%
  dplyr::filter(Sample.Name %in% seurat_data$Sample.Name)


full_seurat_meta <- seurat_data[[]] %>%
  dplyr::select(Sample.Name) %>%
  tibble::rownames_to_column("barcode")

sample_metadata <- merge(sample_metadata, full_seurat_meta, all.x = FALSE,
                         all.y = TRUE)

sample_metadata <- sample_metadata %>%
  tibble::column_to_rownames("barcode") %>%
  dplyr::select(!Sample.Name)

seurat_data <- AddMetaData(seurat_data, metadata = sample_metadata)

# Colors -----------------------------------------------------------------------
all_samples <- unique(seurat_data$sample)
sample_colors <- MetBrewer::met.brewer(name = "Archambault", 
                                       n = length(all_samples),
                                       type = "continuous")

names(sample_colors) <- all_samples

antigen_colors <- MetBrewer::met.brewer(name = "Demuth", n = 4,
                                        type = "discrete")

names(antigen_colors) <- c("Negative", "Tet_antigen", "other",
                           "diabetes_antigen")

status_colors <- MetBrewer::met.brewer(name = "Egypt", n = 3,
                                       type = "discrete")

names(status_colors) <- c("nd", "no", "aab")


status_colors2 <- MetBrewer::met.brewer(name = "Egypt", n = 4,
                                        type = "discrete")

names(status_colors2) <- c("nd", "no", "aab_stage_1", "aab_stage_2")

celltype_colors <- readRDS(file.path(save_dir, "color_palette.rds"))

# Update seurat object ---------------------------------------------------------

# First remove all non-b cells
celltypes_keep <- c("Naive_1", "Naive_2", "Naive_3",
                    "BND_cluster", "B.intermediate",
                    "Early_memory", "Memory_IgE_IgG", "Resting_memory",
                    "DN2", "Early_memory", "BND2",
                    "Memory_IgA", "Activated_memory", "Activated_naive")

seurat_data <- subset(seurat_data, subset = RNA_celltype %in% celltypes_keep)

# Any diabetes antigen and no binding
seurat_data$test_id <- ifelse(seurat_data$hash.ID %in% c("INS-tet",
                                                         "GAD-tet",
                                                         "IA2-tet"),
                              "diabetes_antigen", ifelse(seurat_data$hash.ID 
                                                         == "Negative",
                                                         "Negative", "other"))

seurat_data$frozen_cell_bin <- ifelse(seurat_data$millions.of.cells.frozen < 11,
                                      "10_or_less",
                                      ifelse(is.na(seurat_data$millions.of.cells.frozen),
                                             "unknown", "10_or_more"))

Idents(seurat_data) <- "test_id"

name_mapping <- c("INS-tet" = "diabetes_antigen",
                  "GAD-tet" = "diabetes_antigen",
                  "IA2-tet" = "diabetes_antigen",
                  "TET-tet" = "Tet_antigen",
                  "Negative" = "Negative",
                  "Doublet" = "other",
                  "DNA-tet" = "other")

seurat_data$test_id <- name_mapping[as.character(seurat_data$hash.ID)]

Idents(seurat_data) <- "test_id"

seurat_data$Status <- gsub(" ", "_", seurat_data$Status)

seurat_data$test_full <- paste(seurat_data$Status, seurat_data$test_id,
                               sep = "_")

seurat_data$new_status <- gsub("aab stage [1|2]", "aab", seurat_data$Status)

seurat_data$test_full_new_status <- paste(seurat_data$new_status, seurat_data$test_id,
                                          sep = "_")

# All B cells ------------------------------------------------------------------
se <- as.SingleCellExperiment(seurat_data)

all_res <- scoreMarkers(se, groups = se$test_full_new_status, 
                        block = se$date.processed.for.scSeq)


# sig_res <- all_res$aab_stage_1_diabetes_antigen[order(all_res$aab_stage_1_diabetes_antigen$mean.AUC,
#                                               decreasing=TRUE),]

all_res2 <- scoreMarkers(se, groups = se$Status, 
                         block = se$date.processed.for.scSeq)




sig_res <- all_res2$no[order(all_res2$no$mean.AUC, decreasing=TRUE),]

featDistPlot(seurat_data, geneset = "HLA-C", sep_by = "Status",
             combine = FALSE)

# Consider taking the top 100 genes from each comparison?

all_res3 <- scoreMarkers(se, groups = se$test_id,
                         block = se$date.processed.for.scSeq)


# Make some plots --------------------------------------------------------------


sample_order <- c("JH310-12_AP",
                  "JH310-12_TS",
                  "JH313-15_OF",
                  "JH310-12_MV",
                  "JH310-12_NG",
                  "JH313-15_DB",
                  "JH313-15_MH",
                  "JH310-12_BH",
                  "JH313-15_AG",
                  "JH313-15_NF",
                  "JH313-15_PH",
                  "JH310-12_FG",
                  "JH310-12_EB",
                  "JH310-12_CP",
                  "JH313-15_FD",
                  "JH313-15_JB")

antigen_order <- c("diabetes_antigen", "other", "Negative", "Tet_antigen")

seurat_data$meta_col <- paste(seurat_data$sample, seurat_data$test_id,
                                sep = "_")
meta_df <- seurat_data[[c("sample", "test_id", "Status", "meta_col")]]
meta_ave <- meta_df
rownames(meta_ave) <- NULL
meta_ave <- distinct(meta_ave)
rownames(meta_ave) <- meta_ave$meta_col

meta_ave$sample <- factor(meta_ave$sample, levels = sample_order)
meta_ave$test_id <- factor(meta_ave$test_id, levels = antigen_order)

meta_ave <- meta_ave %>%
  dplyr::arrange(sample, test_id)

color_list <- list("sample" = sample_colors[names(sample_colors) %in% 
                                              unique(meta_df$sample)],
                   "test_id" = antigen_colors,
                   "Status" = status_colors2)

new_meta <- meta_df %>%
  dplyr::select(sample, Status) %>%
  dplyr::distinct()

rownames(new_meta) <- new_meta$sample

new_meta$sample <- factor(new_meta$sample, levels = sample_order)

new_meta <- new_meta %>%
  arrange(sample)

gene_list <- lapply(names(all_res2), function(x){
  sig_res <- all_res2[[x]] %>%
    data.frame() %>%
    dplyr::filter(median.logFC.cohen > 0.1)
  if(nrow(sig_res) > 100){
    sig_res <- sig_res[order(sig_res$median.AUC, decreasing=TRUE),][1:100, ]
  }
  return(rownames(sig_res))
})

gene_list <- unique(unlist(gene_list))

plot_heatmap(seurat_object = seurat_data, gene_list = gene_list,
             meta_df = new_meta, average_expression = TRUE, 
             meta_col = "sample", cluster_rows = TRUE, plot_rownames = FALSE,
             color_list = color_list)

dev.off()

plot_heatmap(seurat_object = seurat_data, gene_list = gene_list,
             meta_df = new_meta, average_expression = TRUE, 
             meta_col = "sample", cluster_rows = TRUE, cluster_cols = TRUE,
             plot_rownames = FALSE,
             color_list = color_list)


featDistPlot(seurat_data, geneset = "MT-CO1", sep_by = "sample", combine = FALSE)


muscat_genes <- read.csv(file.path(save_dir, "files",
                                             "pseudobulk_all_bcell_de.csv"))

keep_genes <- muscat_genes %>%
  dplyr::filter(contrast == "NO - ND")

length(gene_list)

length(keep_genes$gene)

length(intersect(gene_list, keep_genes$gene))
