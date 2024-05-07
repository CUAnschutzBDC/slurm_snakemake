library(Seurat)
library(tidyverse)
library(cowplot)
library(openxlsx)
library(here)
library(scAnalysisR)
library(clustifyr)
library(splitstackshape)
library(viridis)
library(circlize)
library(pheatmap)



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

#sample <- args[[1]]
#sample <- gsub("__.*", "", sample)
sample <- "Npod6456_PLN"

#sample_info <- args[[3]]
sample_info <- here("files/sample_info.tsv")

#results_dir <- args[[2]]
results_dir <- here("results")

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
save_dir <- file.path(results_dir, "R_analysis", sample)

seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))

# Subset to B cells
seurat_data <- subset(seurat_data, subset = updated_celltype %in%
                        c("B_Memory", "B_Naive", "Plasmablast"))

vdj_dir <- file.path(save_dir, "images", "bcell_vdj_plots")

ifelse(!dir.exists(vdj_dir), dir.create(vdj_dir), FALSE)

# Colors -----------------------------------------------------------------------
final_colors <- c("B_Memory" = "#924bdb", # Resting memory
                  "B_Naive" = "#69ba3d", # Naive 1
                  "Plasmablast" = "#ffac14") 

status_colors <- MetBrewer::met.brewer(name = "Egypt", n = 4,
                                       type = "discrete")

seurat_data$Status <- "no"

names(status_colors) <- c("nd", "no", "aab_stage_1", "aab_stage_2")

# Analysis ---------------------------------------------------------------------

sep_columns <- c("bcr_chains", "bcr_cdr3", "bcr_cdr3_length",
                 "bcr_cdr3_nt_length", "bcr_v_gene", "bcr_d_gene",
                 "bcr_j_gene", "bcr_c_gene", "bcr_reads", "bcr_umis",
                 "bcr_productive", "bcr_full_length",
                 "bcr_v_ins", "bcr_v_del", "bcr_v_mis", 
                 "bcr_d_ins", "bcr_d_del", "bcr_d_mis", "bcr_j_ins",
                 "bcr_j_del", "bcr_j_mis", "bcr_c_ins", "bcr_c_del",
                 "bcr_c_mis", "bcr_all_ins", "bcr_all_del",
                 "bcr_all_mis", "bcr_vd_ins", "bcr_vd_del", "bcr_dj_ins",
                 "bcr_dj_del", "bcr_v_mis_freq", "bcr_d_mis_freq",
                 "bcr_j_mis_freq", "bcr_c_mis_freq", "bcr_all_mis_freq")
keep_columns <- c("bcr_isotype", "updated_celltype", "bcr_paired",
                  "Status", "all_chains")

all_info <- seurat_data[[]] %>%
  dplyr::mutate(all_chains = bcr_chains) %>%
  dplyr::select(dplyr::all_of(c(sep_columns, keep_columns)), bcr_n_chains) %>%
  tibble::rownames_to_column("barcode") %>%
  dplyr::filter(!is.na(bcr_n_chains))


all_info_split <- cSplit(all_info, sep_columns, sep = ";", direction = "long") %>%
  dplyr::filter(!is.na(bcr_chains))

colnames(all_info_split) <- gsub("bcr_", "", colnames(all_info_split))

# Circos plots to make
# 1. Heavy V+J pairs
# 2. Light V+J pairs
# 3. Heavy + Light V+J pairs
# 4. Heavy VJ pairs between samples
# 5. Light VJ pairs between samples
# 6. Heavy + Light VJ pairs between samples
# For each, also separate by chain...

# Functions

make_circos_plot <- function(save_name, circos_df, color = NULL,
                             grid_color = NULL){
  if(nrow(circos_df) > 0){
    pdf(save_name, height = 17, width = 17)
    
    chordDiagram(circos_df, annotationTrack = "grid", 
                 preAllocateTracks = 1, col = color,
                 grid.col = grid_color)
    
    # we go back to the first track and customize sector labels
    circos.track(track.index = 1, panel.fun = function(x, y) {
      circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                  facing = "clockwise", niceFacing = TRUE, adj = c(-0.1, 0.5))
    }, bg.border = NA) # here set bg.border to NA is important
    
    dev.off()
    
    graphics.off()
    
    circos.clear()    
  }
  
}

count_genes <- function(starting_df, group_by = "all",
                        chain = "IGH",
                        keep_chains = c("IGH", "IGH;IGK",
                                        "IGH;IGK;IGK", "IGH;IGK;IGL",
                                        "IGH;IGL", "IGH;IGL;IGL"),
                        color_list = NULL, subset_counts = 0){
  if(group_by == "all"){
    select_cols <- c("v_gene", "j_gene")
  } else if (group_by == "sample") {
    select_cols <- c("v_gene", "j_gene", "sample")
  } else if (group_by == "Status") {
    select_cols <- c("v_gene", "j_gene", "Status")
  } else {
    stop("group_by must be 'all', 'sample', or 'Status'")
  }
  
  return_df <- starting_df %>%
    dplyr::filter(all_chains %in% keep_chains) %>%
    dplyr::filter(chains == chain) %>%
    dplyr::select(dplyr::all_of(select_cols))
  
  return_df <- return_df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(select_cols))) %>% 
    dplyr::add_count(name = "value") %>%
    dplyr::distinct() %>%
    dplyr::ungroup() %>%
    dplyr::filter(value >= subset_counts)
  
  if(group_by != "all"){
    return_df$color <- color_list[return_df[[group_by]]]
    color_values <- return_df$color
    return_df <- return_df %>%
      dplyr::select(-dplyr::all_of(group_by), -color)
  } else{
    color_values <- NULL
  }
  
  return(list("df" = return_df, "color_list" = color_values))
}

all_v <- unique(all_info_split$v_gene)
all_j <- unique(all_info_split$j_gene)

palette1 <- colorRampPalette(colors = 
                               RColorBrewer::brewer.pal(name = "Set1", n = 9))

palette2 <- colorRampPalette(colors = 
                               RColorBrewer::brewer.pal(name = "Set2", n = 8))

v_colors <- palette1(length(all_v))
names(v_colors) <- all_v

j_colors <- palette2(length(all_j))
names(j_colors) <- all_j

all_colors <- c(v_colors, j_colors)

graphics.off()
# Circos -----------------------------------------------------------------------

test_circos <- list("all" = list(directory = file.path(vdj_dir, "all_pairs"),
                                 subset_counts = 0),
                    "10_plus" = list(directory = file.path(vdj_dir, 
                                                           "vj_counts_10_plus"),
                                     subset_counts = 10))

for(i in names(test_circos)){
  directory <- test_circos[[i]]$directory
  subset_counts <- test_circos[[i]]$subset_counts
  
  ifelse(!dir.exists(directory), dir.create(directory), FALSE)
  
  ## Heavy V+J pairs -------------------------------------------------------------
  # I'll start with 10x
  all_data <- count_genes(starting_df = all_info_split,
                          group_by = "all",
                          subset_counts = subset_counts)
  
  make_circos_plot(save_name = file.path(directory,
                                         "heavy_v_j_circos.pdf"),
                   circos_df = all_data$df,
                   color = all_data$color_list, 
                   grid_color = all_colors)
  
  ## Isotype switched status ---------------------------------------------------
  isotypes_use <- c("IGHA", "IGHD", "IGHG", "IGHM")
  
  invisible(lapply(isotypes_use, function(x){
    test_data <- all_info_split %>%
      dplyr::filter(isotype == x)
    
    ## All samples
    all_data <- count_genes(starting_df = test_data,
                            group_by = "all",
                            subset_counts = subset_counts)
    
    make_circos_plot(save_name = file.path(directory,
                                           paste0(x, "_heavy_v_j_circos.pdf")),
                     circos_df = all_data$df,
                     color = all_data$color_list, 
                     grid_color = all_colors)
  }))
  
}
graphics.off()


# Heatmap ----------------------------------------------------------------------
make_heatmap_info <- function(all_info_split, select_cols,
                              group_cols = NULL, type = NULL,
                              subset_counts = 0,
                              keep_chains = c("IGH", "IGH;IGK", "IGH;IGK;IGK",
                                              "IGH;IGK;IGL", "IGH;IGL", 
                                              "IGH;IGL;IGL"),
                              chains_use = c("IGH")){
  
  if(is.null(group_cols)){
    group_cols <- select_cols
  }
  all_info_split_heatmap <- all_info_split %>%
    dplyr::filter(all_chains %in% keep_chains) %>%
    dplyr::filter(chains %in% chains_use) %>%
    dplyr::select(dplyr::all_of(select_cols))
  
  all_info_split_heatmap <- all_info_split_heatmap %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>% 
    dplyr::add_count(name = "value") %>%
    dplyr::distinct()
  
  if(is.null(type)){
    all_info_split_heatmap <- all_info_split_heatmap %>%
      dplyr::mutate(j_group = j_gene) %>%
      ungroup()
  } else if(type == "sample"){
    all_info_split_heatmap <- all_info_split_heatmap %>%
      dplyr::mutate(j_group = paste(j_gene, sample, sep = "_")) %>%
      ungroup()
  } else if (type == "status") {
    all_info_split_heatmap <- all_info_split_heatmap %>%
      dplyr::mutate(j_group = paste(j_gene, Status, sep = "_")) %>%
      ungroup()
  } else {
    stop("type arugment only excepts NULL, status or sample")
  }
  
  new_info <- all_info_split_heatmap %>%
    dplyr::select(v_gene, j_group, value) %>%
    tidyr::pivot_wider(names_from = j_group, values_from = value) %>%
    tibble::column_to_rownames("v_gene")
  
  new_info[is.na(new_info)] <- 0
  
  # Find max for each row
  row_max <- apply(new_info, 1, max)
  
  new_info <- new_info[row_max > subset_counts, ]
  
  if (is.null(type)){
    sample_info <- NULL
  } else if(type %in% c("sample", "status")){
    select_cols_sample_info <- c("j_group", select_cols)
    select_cols_sample_info <- 
      select_cols_sample_info[select_cols_sample_info != "v_gene"]
    sample_info <- all_info_split_heatmap %>%
      dplyr::select(dplyr::all_of(select_cols_sample_info)) %>%
      dplyr::distinct() %>%
      tibble::column_to_rownames("j_group")
    
    if(!identical(rownames(sample_info), colnames(new_info))){
      new_info <- new_info[ , rownames(sample_info)]
    }  
  }
  
  return(list(heatmap_plot = new_info, sample_info = sample_info))
  
}

plot_heatmap <- function(new_info, sample_info = NULL,
                         save_name = NULL, coloring = NULL){
  pdf(save_name, height = 10, width = 10)
  
  if(nrow(new_info) > 0){
    if(nrow(new_info) > 1){
      cluster_rows <- TRUE
    } else {
      cluster_rows <- FALSE
    }
    if(is_null(sample_info)){
      pheatmap(new_info, color = magma(n = 10),
               cluster_rows = cluster_rows)
      
      grid::grid.newpage()
      new_info <- t(scale(t(new_info)))
      
      pheatmap(new_info, color = magma(n = 10), cluster_rows = cluster_rows)    
    } else {
      pheatmap(new_info, color = magma(n = 10), annotation_col = sample_info,
               annotation_colors = coloring, cluster_rows = cluster_rows)
      
      grid::grid.newpage()
      new_info <- t(scale(t(new_info)))
      
      pheatmap(new_info, color = magma(n = 10), annotation_col = sample_info,
               annotation_colors = coloring, cluster_rows = cluster_rows)
    }
    
  } else {
    return(NULL)
  }
  
  
  dev.off()
  graphics.off()
  
}

for(i in names(test_circos)){
  directory <- test_circos[[i]]$directory
  subset_counts <- test_circos[[i]]$subset_counts
  
  heatmap_res <- make_heatmap_info(all_info_split, 
                                   select_cols = c("v_gene", "j_gene"),
                                   subset_counts = subset_counts)
  
  plot_heatmap(new_info = heatmap_res$heatmap_plot, sample_info = NULL,
               save_name = file.path(directory, "heavy_v_j_heatmap.pdf"),
               coloring = NULL)
  
  ## Isotype switched status ---------------------------------------------------
  isotypes_use <- c("IGHA", "IGHD", "IGHG", "IGHM")
  
  invisible(lapply(isotypes_use, function(x){
    test_data <- all_info_split %>%
      dplyr::filter(isotype == x)
    
    ## All samples
    heatmap_res <- make_heatmap_info(test_data, 
                                     select_cols = c("v_gene", "j_gene"),
                                     subset_counts = subset_counts)
    
    plot_heatmap(new_info = heatmap_res$heatmap_plot, sample_info = NULL,
                 save_name = file.path(directory,
                                       paste0(x, "heavy_v_j_heatmap.pdf")),
                 coloring = NULL)
  }))
  
}
graphics.off()

# Barplots ---------------------------------------------------------------------
# Here I want to make barplots that are the percent of the repertoire for 
# each VJ pair. 
# Steps
# 1. Count all vj pairs within each grouping (sample, tetramer, isotype, status)
# 2. Count all cells within each grouping (sample, tetramer, isotype, status)
# 3. Use these two to find percent
# 4. Plot percent -> try to keep colors consistent between all plots...?
# 5. Repeat this subsetting to sample (group tetramer, isotype, status),
# tetramer (group sample, isotype, status), isotype (group sample, tetramer,
# status), status (group sample, tetramer, isotype) 
# I'm not certain on the subsetting front at the moment
keep_chains <- c("IGH", "IGH;IGK", "IGH;IGK;IGK",
                 "IGH;IGK;IGL", "IGH;IGL", 
                 "IGH;IGL;IGL")

test_set <- all_info_split %>%
  dplyr::filter(chains %in% c("IGH")) %>% 
  dplyr::filter(all_chains %in% keep_chains) %>%
  dplyr::select(v_gene, j_gene, Status) %>%
  dplyr::mutate(vj_gene = paste(v_gene, j_gene, sep = "_")) %>%
  dplyr::group_by(vj_gene) %>%
  dplyr::add_count(name = "sample_vj_count") %>%
  dplyr::ungroup() %>%
  dplyr::add_count(name = "sample_count") %>%
  dplyr::ungroup() %>%
  dplyr::group_by(vj_gene) %>%
  dplyr::add_count(name = "full_vj_count") %>%
  dplyr::ungroup() %>%
  dplyr::distinct() %>%
  dplyr::mutate(percent_vj = sample_vj_count / sample_count * 100,
                percent_sample = sample_vj_count / full_vj_count * 100)


pdf(file.path(vdj_dir, "vj_barplots_test.pdf"),
    height = 10, width = 20)
# print(ggplot2::ggplot(test_set, ggplot2::aes(x = sample, y = percent_vj,
#                                        fill = vj_gene)) +
#   ggplot2::geom_bar(stat = "identity", position = "stack"))

print(ggplot2::ggplot(test_set, ggplot2::aes(x = vj_gene, y = sample_vj_count,
                                             fill = sample)) +
        ggplot2::geom_bar(stat = "identity", position = "stack") +
        ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
        ggplot2::ggtitle("VJ count colored by sample"))

print(ggplot2::ggplot(test_set, ggplot2::aes(x = vj_gene, y = sample_vj_count,
                                             fill = Status)) +
        ggplot2::geom_bar(stat = "identity", position = "stack") +
        ggplot2::scale_fill_manual(values = status_colors) +
        ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggplot2::ggtitle("VJ count colored by status"))

dev.off()


# Make a data frame
# Need
# 1. Percent of the vj gene (Here, adding across samples should be 100%)
# 2. Percent of the sample (Here adding across vj genes should be 100%)
# 3. Number of samples
# 4. Number of samples per status
# Do for V and VJ

vdj_data <- all_info_split %>%
  dplyr::filter(chains %in% c("IGH")) %>% 
  dplyr::filter(all_chains %in% keep_chains) %>%
  dplyr::select(v_gene, j_gene, Status) %>%
  dplyr::mutate(vj_gene = paste(v_gene, j_gene, sep = "_")) %>%
  dplyr::group_by(vj_gene) %>%
  dplyr::add_count(name = "sample_vj_count") %>%
  dplyr::ungroup() %>%
  dplyr::add_count(name = "sample_count") %>%
  dplyr::ungroup() %>%
  dplyr::group_by(vj_gene) %>%
  dplyr::add_count(name = "full_vj_count") %>%
  dplyr::ungroup() %>%
  dplyr::distinct() %>%
  dplyr::mutate(percent_vj = sample_vj_count / sample_count * 100,
                percent_sample = sample_vj_count / full_vj_count * 100)


# Want percent vj
# 
# 
# fig <- plot_ly(
#   type = 'scatterpolar',
#   mode = 'lines'
# ) 
# 
# for(i in colnames(test_data)){
#   theta_vals <- c(test_data$v_gene, test_data$v_gene[[1]])
#   if(i != "vj_gene"){
#     fig <- fig %>%
#       add_trace(
#         r = c(test_data[[i]], test_data[[i]][[1]]),
#         theta = theta_vals,
#         name = i
#       )
#   }
# }
