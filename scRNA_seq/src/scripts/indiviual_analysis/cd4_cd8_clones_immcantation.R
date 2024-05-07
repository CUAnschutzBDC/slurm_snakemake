library(Seurat)
library(tidyverse)
library(cowplot)
library(openxlsx)
library(here)
library(scAnalysisR)
library(clustifyr)
library(splitstackshape)
library(data.table)


# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log


cell_types <- "combined_celltype"
clusters <- "combined_cluster"
pval <- 0.05
logfc <- 0.5

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

args <- commandArgs(trailingOnly = TRUE)

sample <- args[[1]]
sample <- gsub("__.*", "", sample)
#sample <- "Npod6456_PLN"
#sample <- "Npod6553_PLN"

sample_info <- args[[3]]
#sample_info <- here("files/sample_info.tsv")

results_dir <- args[[2]]
#results_dir <- here("results")

sample_info <- read.table(sample_info, fill = TRUE, header = TRUE)

sample_info <- sample_info[sample_info$sample == sample,]

HTO <- sample_info$HTO
ADT <- sample_info$ADT
hash_ident <- sample_info$hash_ident

ADT_pca <- sample_info$adt_pca
run_adt_umap <- sample_info$adt_umap

RNA_pcs <- sample_info$PCs
resolution <- sample_info$resolution

if(normalization_method == "SCT"){
  SCT <- TRUE
  seurat_assay <- "SCT"
} else {
  SCT <- FALSE
  seurat_assay <- "RNA"
}

# Set directories
save_dir <- file.path(results_dir, "R_analysis")

# Read in data
cd4_save_dir <- file.path(save_dir, paste0(sample, "_cd4"))
cd8_save_dir <- file.path(save_dir, paste0(sample, "_cd8"))
seurat_cd4 <- readRDS(file.path(cd4_save_dir, "rda_obj", "seurat_processed.rds"))
seurat_cd8 <- readRDS(file.path(cd8_save_dir, "rda_obj", "seurat_processed.rds"))
immcantation_clones <- read.table(file.path(results_dir, sample, "outs",
                                          "immcantation", "tcr",
                                          "filtered_contig_igblast_db-pass.tsv"),
                                row.names = 1, header = TRUE, sep = "\t")

immcantation_cd4 <- immcantation_clones[immcantation_clones$cell_id %in%
                                          colnames(seurat_cd4),]

immcantation_cd8 <- immcantation_clones[immcantation_clones$cell_id %in%
                                          colnames(seurat_cd8),]

if(sample == "Npod6553_PLN"){
  # Remove any clones with low UMI count as they could be background
  immcantation_cd4 <- immcantation_cd4[immcantation_cd4$umi_count > 1,]
  immcantation_cd8 <- immcantation_cd8[immcantation_cd8$umi_count > 1,]
  

}

order_cells <- function(immcantation_df){
  keep_cols <- c("productive", "v_call", "d_call", "j_call",
                 "junction", "junction_aa", "locus", "c_call",
                 "fwr1", "fwr2", "fwr3", "fwr4", "cdr1", "cdr2",
                 "cdr3", "cell_id")
  
  
  `%notin%` <- Negate(`%in%`)
  pivot_cols <- keep_cols[keep_cols %notin% c("locus", "cell_id")]
  
  
  
  immcantation_df <- immcantation_df[ , keep_cols]
  
  return_df <- immcantation_df %>%
    group_by(cell_id, locus) %>%
    arrange(locus) %>%
    mutate(locus = paste(locus, row_number(), sep = "_")) %>%
    tidyr::pivot_wider(names_from = locus,
                       names_sep = "-",
                       values_from = dplyr::all_of(pivot_cols),
                       id_cols = cell_id)
  
  grouping_columns <- c(colnames(return_df)[grepl("v_call",
                                            colnames(return_df))],
                        colnames(return_df)[grepl("j_call",
                                            colnames(return_df))], 
                        colnames(return_df)[grepl("junction",
                                            colnames(return_df))])
  
  return_df <- return_df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(grouping_columns))) %>%
    dplyr::arrange() %>%
    dplyr::mutate(clone_id = cur_group_id())
  
  col_order <- c("cell_id", "clone_id",
                 colnames(return_df)[grepl("TRA_1", colnames(return_df))],
                 colnames(return_df)[grepl("TRA_2", colnames(return_df))],
                 colnames(return_df)[grepl("TRB_1", colnames(return_df))],
                 colnames(return_df)[grepl("TRB_2", colnames(return_df))])
  
  
  return_df <- return_df[, col_order]
}

cd4_df <- order_cells(immcantation_cd4)
cd8_df <- order_cells(immcantation_cd8)

save_wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb = save_wb, sheetName = "cd8")
openxlsx::writeData(wb = save_wb, sheet = "cd8", x = cd8_df)

openxlsx::addWorksheet(wb = save_wb, sheetName = "cd4")
openxlsx::writeData(wb = save_wb, sheet = "cd4", x = cd4_df)

openxlsx::saveWorkbook(wb = save_wb, file = file.path(save_dir, sample, "files",
                                                      paste0(sample, 
                                                             "_cd4_cd8_cells.xlsx")),
                       overwrite = TRUE)

