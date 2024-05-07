library(Seurat)
library(tidyverse)
library(here)
library(scAnalysisR)
library(djvdj)
library(UpSetR)
library(muscat)
library(SingleCellExperiment)

source(here("src/scripts/muscat_plotting_functions.R"))

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

sample <- "merged"

all_ref_dir <-
  "/beevol/home/wellskri/Analysis/references/single_cell_references"

HTO <- FALSE
ADT <- FALSE

if(normalization_method == "SCT"){
  SCT <- TRUE
  seurat_assay <- "SCT"
} else {
  SCT <- FALSE
  seurat_assay <- "RNA"
}

# Set directories
base_dir <- here()

save_dir <- file.path(base_dir, "results", "R_analysis", sample)
fig_dir <- file.path(save_dir, "images", "figures")

ifelse(!dir.exists(fig_dir), dir.create(fig_dir), FALSE)

# Read in the data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed_no_doublet.rds"))

plot_tetramers <- c("INS-tet", "TET-tet", "IA2-tet", "GAD-tet", "DNA-tet")

# Colors -----------------------------------------------------------------------
final_colors <- c("Resting_memory" = "#924bdb", # Resting memory
                  "Naive_1" = "#69ba3d", # Naive 1
                  "Naive_2" = "#9a43a4", # Naive 2
                  "Memory_IgE_IgG" = "#bf9b31", # Memory IgE/IgG1
                  "Naive_3" = "#6477ce", # Naive 3
                  "Memory_IgA" = "#d15131", # Memory IA
                  "Early_memory" = "#4c9e8e", # Early Memory
                  "BND2" = "#cc4570", #Bnd2
                  "DN2" = "#648d4f", # DN2
                  "Activated_memory" = "#985978", # Activated memory
                  "Activated_naive" = "#a06846", # Activated naive
                  "B.intermediate" = "#00008b",
                  "CD14.Mono" = "#e0205a",
                  "pDC" = "#ffb6d3",
                  "Plasmablast" = "#ffac14",
                  "CD8.TEM" = "#000000") 

tetramer_colors <- MetBrewer::met.brewer(name = "Juarez", n = 7,
                                         type = "continuous")
names(tetramer_colors) <- c("Negative", "Doublet", "INS-tet", "TET-tet",
                            "IA2-tet", "GAD-tet", "DNA-tet")

sample_colors <- MetBrewer::met.brewer(name = "Archambault", n = 16,
                                       type = "continuous")


all_samples <- unique(seurat_data$sample)

names(sample_colors) <- all_samples

antigen_colors <- MetBrewer::met.brewer(name = "Demuth", n = 4,
                                        type = "discrete")

names(antigen_colors) <- c("Negative", "Tet_antigen", "other",
                           "diabetes_antigen")


status_colors <- MetBrewer::met.brewer(name = "Egypt", n = 4,
                                       type = "discrete")

seurat_data$Status <- gsub(" ", "_", seurat_data$Status)

names(status_colors) <- c("nd", "no", "aab_stage_1", "aab_stage_2")
# Set levels -------------------------------------------------------------------
seurat_data$hash.ID <- factor(seurat_data$hash.ID,
                              levels = names(tetramer_colors))

celltype_levels <- c("Naive_1", "Naive_3",
                     "BND2", "B.intermediate",
                     "Memory_IgA", "Resting_memory", 
                     "Activated_memory",
                     "Plasmablast",
                     "CD14.Mono", "CD8.TEM")

seurat_data$RNA_combined_celltype <- factor(seurat_data$RNA_combined_celltype,
                                            levels = celltype_levels)

# Make seurat subsets ----------------------------------------------------------
seurat_no_dub <- subset(seurat_data, 
                        subset = hash.ID != "Doublet")

seurat_b_cell <- subset(seurat_data, 
                        subset = RNA_combined_celltype %in% celltype_levels[1:8])

seurat_b_cell$RNA_combined_celltype <- factor(seurat_b_cell$RNA_combined_celltype,
                                       levels = celltype_levels[1:8])

seurat_b_no_dub <- subset(seurat_b_cell, subset = hash.ID != "Doublet")

# Celltype plots ---------------------------------------------------------------
# UMAP of celltypes
celltype_plot <- plotDimRed(seurat_data, col_by = "RNA_combined_celltype",
                            color = final_colors, plot_type = "mnn.umap",
                            ggrastr = TRUE, size = 0.25)[[1]]

pdf(file.path(fig_dir, "celltype_umap.pdf"), height = 8, width = 8)

print(celltype_plot)

dev.off()

# UMAP of celltypes separated by individual
all_celltype_plots <- lapply(all_samples, function(x){
  return_plot <- plotDimRed(seurat_data, col_by = "RNA_combined_celltype",
                            color = final_colors, highlight_group = TRUE,
                            group = x, meta_data_col = "orig.ident",
                            plot_type = "mnn.umap", ggrastr = TRUE,
                            size = 0.25)[[1]]
}) 

pdf(file.path(fig_dir, "sample_celltype_umap.pdf"), height = 20, width = 10)

print(cowplot::plot_grid(plotlist = all_celltype_plots, nrow = 6, ncol = 2))

dev.off()

# Antigen plots ----------------------------------------------------------------
# Violin plots with the doublet
all_plots1 <- featDistPlot(seurat_data, geneset = plot_tetramers, 
                           combine = FALSE, sep_by = "hash.ID",
                           assay = "ADT_norm", color = tetramer_colors)

pdf(file.path(fig_dir, "antigen_violins.pdf"), height = 5, width = 10)

print(cowplot::plot_grid(plotlist = all_plots1, nrow = 2, ncol = 3))

dev.off()

# Violin plots without the doublet
all_plots1 <- featDistPlot(seurat_no_dub, geneset = plot_tetramers, 
                           combine = FALSE, sep_by = "hash.ID",
                           assay = "ADT_norm", color = tetramer_colors)

pdf(file.path(fig_dir, "antigen_violins_no_doublet.pdf"), height = 5, 
    width = 10)

print(cowplot::plot_grid(plotlist = all_plots1, nrow = 2, ncol = 3))

dev.off()


# Bar plots of binding by individual, status, cell type. Do cell type individually
# separated by status, individual. With and without doublet, non-b cells

make_barplots <- function(full_seurat, b_seurat, save_name = "all"){
  # Barplots separated by individual
  sample_barplot <- scAnalysisR::stacked_barplots(full_seurat,
                                                  meta_col = "hash.ID",
                                                  split_by = "sample", 
                                                  color = tetramer_colors) +
    ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  diabetes_barplot <- scAnalysisR::stacked_barplots(full_seurat,
                                                    meta_col = "hash.ID",
                                                    split_by = "Status", 
                                                    color = tetramer_colors) +
    ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  cluster_barplot <- scAnalysisR::stacked_barplots(full_seurat,
                                                   meta_col = "hash.ID",
                                                   split_by = "RNA_combined_celltype", 
                                                   color = tetramer_colors) +
    ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  bcell_barplot <- scAnalysisR::stacked_barplots(b_seurat,
                                                 meta_col = "hash.ID",
                                                 split_by = "RNA_combined_celltype", 
                                                 color = tetramer_colors) +
    ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  pdf(file.path(fig_dir, paste0("sample_antigen_barplot_",
                                save_name, ".pdf")), 
      width = 4, height = 4)
  print(sample_barplot)
  
  dev.off()
  
  pdf(file.path(fig_dir, paste0("diabetes_antigen_barplot_", 
                                save_name, ".pdf")), 
      width = 3, height = 4)
  
  print(diabetes_barplot)
  
  dev.off()
  
  pdf(file.path(fig_dir, paste0("all_celltype_antigen_barplot_", 
                                save_name, ".pdf")), 
      width = 4, height = 4)
  print(cluster_barplot)
  
  dev.off()
  
  pdf(file.path(fig_dir, paste0("b_celltype_antigen_barplot_",
                                save_name, ".pdf")), 
      width = 4, height = 4)
  print(bcell_barplot)
  
  dev.off()
  
  # Barplots separated by celltype
  celltype_barplots <- lapply(unique(full_seurat$RNA_combined_celltype), function(x){
    seurat_subset <- subset(full_seurat, subset = RNA_combined_celltype == x)
    save_barplot <- scAnalysisR::stacked_barplots(seurat_subset,
                                                  meta_col = "hash.ID",
                                                  split_by = "sample", 
                                                  color = tetramer_colors) +
      ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ggplot2::ggtitle(x)
  })
  
  pdf(file.path(fig_dir, paste0("individual_celltype_barplot_by_sample_", 
                                save_name, ".pdf")),
      width = 4, height = 4)
  
  print(celltype_barplots)
  
  dev.off()
  
  diabetes_barplots <- lapply(unique(full_seurat$RNA_combined_celltype), function(x){
    seurat_subset <- subset(full_seurat, subset = RNA_combined_celltype == x)
    save_barplot <- scAnalysisR::stacked_barplots(seurat_subset,
                                                  meta_col = "hash.ID",
                                                  split_by = "Status", 
                                                  color = tetramer_colors) +
      ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ggplot2::ggtitle(x)
  })
  
  pdf(file.path(fig_dir, paste0("individual_celltype_barplot_by_status_", 
                                save_name, ".pdf")),
      width = 3, height = 4)
  
  print(diabetes_barplots)
  
  dev.off()
}

make_barplots(full_seurat = seurat_data, b_seurat = seurat_b_cell,
              save_name = "all")

make_barplots(full_seurat = seurat_no_dub, b_seurat = seurat_b_no_dub,
              save_name = "no_doublet")

# Mutation plots ---------------------------------------------------------------
sep_columns <- c("chains", "cdr3", "cdr3_length",
                 "cdr3_nt_length", "v_gene", "d_gene", "j_gene", "c_gene",
                 "reads", "umis", "productive", "full_length",
                 "v_ins", "v_del", "v_mis", "d_ins", "d_del",
                 "d_mis", "j_ins", "j_del", "j_mis", "c_ins", "c_del", "c_mis",
                 "all_ins", "all_del", "all_mis", "vd_ins", "vd_del", "dj_ins",
                 "dj_del", "v_mis_freq", "d_mis_freq", "j_mis_freq",
                 "c_mis_freq", "all_mis_freq")
keep_columns <- c("isotype", "RNA_combined_celltype", "sample", "paired",
                  "clonotype_id", "Status", "hash.ID")

all_info <- seurat_b_cell[[]] %>%
  dplyr::select(dplyr::all_of(c(sep_columns, keep_columns)), n_chains)


# First try only with n_chains = 2
all_info <- all_info %>%
  dplyr::filter(n_chains == 2)

all_data <- lapply(sep_columns, function(x){
  sep_info <- all_info %>%
    dplyr::select(dplyr::all_of(x)) %>%
    dplyr::rename("cell" = dplyr::all_of(x)) %>%
    tidyr::separate(cell, c("chain1", "chain2"), sep = ";")
  
  colnames(sep_info) <- paste(x, colnames(sep_info), sep = "_")
  
  return(sep_info)
})

all_data <- do.call(cbind, all_data)

keep_cols <- all_info %>%
  dplyr::select(dplyr::all_of(keep_columns))

all_data <- cbind(keep_cols, all_data)


all_data$all_mis_freq_chain1 <- as.double(all_data$all_mis_freq_chain1)
all_data$all_mis_freq_chain2 <- as.double(all_data$all_mis_freq_chain2)

isotype_plot_h <- ggplot2::ggplot(all_data, ggplot2::aes(x = all_mis_freq_chain1,
                                                         y = isotype,
                                                         fill = isotype)) +
  ggridges::geom_density_ridges() +
  ggplot2::scale_fill_brewer(palette = "Set1") +
  ggplot2::ggtitle("Frequency of micmatches in heavy chain") +
  ggplot2::xlab("Frequency of mismatches")


celltype_plot_h <- ggplot2::ggplot(all_data, ggplot2::aes(x = all_mis_freq_chain1,
                                                          y = RNA_combined_celltype,
                                                          fill = RNA_combined_celltype)) +
  ggridges::geom_density_ridges() +
  ggplot2::scale_fill_manual(values = final_colors) +
  ggplot2::xlim(c(-0.02, 0.2)) +
  ggplot2::ggtitle("Frequency of mismatches in heavy chain")

sample_plot_h <- ggplot2::ggplot(all_data, ggplot2::aes(x = all_mis_freq_chain1,
                                                        y = sample,
                                                        fill = sample)) +
  ggridges::geom_density_ridges() +
  ggplot2::scale_fill_manual(values = sample_colors) +
  ggplot2::xlim(c(-0.02, 0.2)) +
  ggplot2::ggtitle("Frequency of mismatches in heavy chain")

status_plot_h <- ggplot2::ggplot(all_data, ggplot2::aes(x = all_mis_freq_chain1,
                                                        y = Status,
                                                        fill = Status)) +
  ggridges::geom_density_ridges() +
  ggplot2::scale_fill_manual(values = status_colors) +
  ggplot2::xlim(c(-0.02, 0.2)) +
  ggplot2::ggtitle("Frequency of mismatches in heavy chain")

tetramer_plot_h <- ggplot2::ggplot(all_data, ggplot2::aes(x = all_mis_freq_chain1,
                                                          y = hash.ID,
                                                          fill = hash.ID)) +
  ggridges::geom_density_ridges() +
  ggplot2::scale_fill_manual(values = tetramer_colors) +
  ggplot2::xlim(c(-0.02, 0.2)) +
  ggplot2::ggtitle("Frequency of mismatches in heavy chain")

isotype_plot_l <- ggplot2::ggplot(all_data, ggplot2::aes(x = all_mis_freq_chain2,
                                                         y = isotype,
                                                         fill = isotype)) +
  ggridges::geom_density_ridges() +
  ggplot2::scale_fill_brewer(palette = "Set1") +
  ggplot2::ggtitle("Frequency of micmatches in light chain") +
  ggplot2::xlab("Frequency of mismatches")


celltype_plot_l <- ggplot2::ggplot(all_data, ggplot2::aes(x = all_mis_freq_chain2,
                                                          y = RNA_combined_celltype,
                                                          fill = RNA_combined_celltype)) +
  ggridges::geom_density_ridges() +
  ggplot2::scale_fill_manual(values = final_colors) +
  ggplot2::xlim(c(-0.02, 0.2)) +
  ggplot2::ggtitle("Frequency of mismatches in light chain")

status_plot_l <- ggplot2::ggplot(all_data, ggplot2::aes(x = all_mis_freq_chain2,
                                                        y = Status,
                                                        fill = Status)) +
  ggridges::geom_density_ridges() +
  ggplot2::scale_fill_manual(values = status_colors) +
  ggplot2::xlim(c(-0.02, 0.2)) +
  ggplot2::ggtitle("Frequency of mismatches in light chain")

tetramer_plot_l <- ggplot2::ggplot(all_data, ggplot2::aes(x = all_mis_freq_chain2,
                                                          y = hash.ID,
                                                          fill = hash.ID)) +
  ggridges::geom_density_ridges() +
  ggplot2::scale_fill_manual(values = tetramer_colors) +
  ggplot2::xlim(c(-0.02, 0.2)) +
  ggplot2::ggtitle("Frequency of mismatches in light chain")

pdf(file.path(fig_dir, "mutation_frequency.pdf"))
print(isotype_plot_h)
print(celltype_plot_h)
print(sample_plot_h)
print(status_plot_h)
print(tetramer_plot_h)
print(isotype_plot_l)
print(celltype_plot_l)
print(status_plot_l)
print(tetramer_plot_l)
dev.off()

# 
# seurat_new <- seurat_b_cell
# 
# seurat_new <- AddMetaData(seurat_new, all_data)
# 
# featDistPlot(seurat_new, geneset = "all_mis_freq_chain1",
#              sep_by = "final_celltype",
#              col_by = "Status", combine = FALSE)

# Clone plots ------------------------------------------------------------------
# Clones by sample

# Clones by cell type (with NA, without NA, with monotytes, T cells,
# without monocytes T cells)

# Clones with at least 3 clones by Status and antigen --> percent of total cells

# Gene list plots --------------------------------------------------------------
gene_lists <- readRDS(here("files/mia_gene_lists.rds"))
plot_tetramers <- c("INS-tet", "TET-tet", "IA2-tet", "GAD-tet", "DNA-tet")


ifelse(!dir.exists(file.path(fig_dir, "gene_list_heatmaps")),
       dir.create(file.path(fig_dir, "gene_list_heatmaps")), FALSE)

# DE plots ---------------------------------------------------------------------
ifelse(!dir.exists(file.path(fig_dir, "de_heatmaps")),
       dir.create(file.path(fig_dir, "de_heatmaps")), FALSE)
# First remove all non-b cells
celltypes_keep <- c("Naive_1", "Naive_2", "Naive_3",
                    "BND_cluster", "B.intermediate",
                    "Early_memory", "Memory_IgE_IgG", "Resting_memory",
                    "Plasmablast")
seurat_data <- subset(seurat_data, subset = RNA_combined_celltype %in% celltypes_keep)

Idents(seurat_data) <- "test_id"


# Any diabetes antigen and no binding
name_mapping <- c("INS-tet" = "diabetes_antigen",
                  "GAD-tet" = "diabetes_antigen",
                  "IA2-tet" = "diabetes_antigen",
                  "TET-tet" = "Tet_antigen",
                  "Negative" = "Negative",
                  "Doublet" = "other",
                  "DNA-tet" = "other")

seurat_data$test_id <- name_mapping[as.character(seurat_data$hash.ID)]

make_de_plots <- function(all_res, save_name, vsd, seurat_object){
  genes_plot <- lapply(unique(all_res$contrast), function(x){
    return_res <- all_res %>%
      dplyr::filter(contrast == x & abs(logFC) > 0.25)
    
    return(return_res$gene)
  })
  
  names(genes_plot) <- unique(all_res$contrast)
  
  pdf(file.path(fig_dir, paste0("upset_plot_pseudobulk_", 
                                save_name, "_de.pdf")))
  print(upset(fromList(genes_plot), order.by = "freq"))
  dev.off()
  
  # Heatmap of DE genes
  # Change sample order
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
  
  seurat_object$meta_col <- paste(seurat_object$sample, seurat_object$test_id,
                                  sep = "_")
  meta_df <- seurat_object[[c("sample", "test_id", "Status", "meta_col")]]
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
                     "Status" = status_colors)
  
  
  pdf(file.path(fig_dir, "de_heatmaps", paste0(save_name, "_no_nd_all.pdf")))
  print(make_corrected_heatamp(vsd = vsd, meta_ave = meta_ave,
                               gene_list = genes_plot$`NO - ND`,
                               meta_col = "meta_col", plot_meta_col = FALSE,
                               max_val = 2.5, min_val = -2.5,
                               cluster_rows = TRUE, cluster_cols = FALSE,
                               coloring = color_list, plot_rownames = FALSE))
  dev.off()
  
  # new_meta <- meta_df %>%
  #   dplyr::select(sample, Status) %>%
  #   dplyr::distinct()
  # 
  # rownames(new_meta) <- new_meta$sample
  # 
  # new_meta$sample <- factor(new_meta$sample, levels = sample_order)
  # 
  # new_meta <- new_meta %>%
  #   arrange(sample)
  # 
  # plot_heatmap(seurat_object = seurat_object, gene_list = genes_plot$`NO - ND`,
  #              meta_df = new_meta, average_expression = TRUE, 
  #              meta_col = "sample", cluster_rows = TRUE, plot_rownames = FALSE,
  #              color_list = color_list)
  # dev.off()
  # plot_heatmap(seurat_object = seurat_object, gene_list = genes_plot$`NO - ND`,
  #              meta_df = new_meta, average_expression = TRUE, 
  #              meta_col = "sample", cluster_rows = TRUE, plot_rownames = FALSE,
  #              color_list = color_list, cluster_cols = TRUE)
  # 
  # plot only the diabetes antigen
  subset_meta <- meta_ave %>%
    dplyr::filter(test_id == "diabetes_antigen") 
  
  subset_meta$sample <- factor(subset_meta$sample, levels = sample_order)
  
  subset_meta <- subset_meta %>%
    dplyr::arrange(Status, sample)
  
  pdf(file.path(fig_dir, "de_heatmaps", paste0(save_name,
                                               "_no_nd_diabetes_antigen.pdf")))
  
  print(make_corrected_heatamp(vsd = vsd, meta_ave = subset_meta,
                               gene_list = genes_plot$`NO_diabetes_antigen - ND_diabetes_antigen`,
                               meta_col = "meta_col", plot_meta_col = FALSE,
                               max_val = 2.5, min_val = -2.5,
                               cluster_rows = TRUE, cluster_cols = FALSE,
                               coloring = color_list, plot_rownames = FALSE))
  
  dev.off()
  
  # plot only the diabetes antigen genes with tet
  subset_meta <- meta_ave %>%
    dplyr::filter(test_id %in% c("Tet_antigen", "diabetes_antigen"))
  
  subset_meta$sample <- factor(subset_meta$sample, levels = sample_order)
  subset_meta$test_id <- factor(subset_meta$test_id, levels = antigen_order)
  
  subset_meta <- subset_meta %>%
    dplyr::arrange(Status, test_id, sample)
  
  subset_meta2 <- subset_meta %>%
    dplyr::arrange(Status, sample, test_id)
  
  pdf(file.path(fig_dir, "de_heatmaps",
                paste0(save_name, "_no_nd_diabetes_antigen_with_tet.pdf")))
  
  print(make_corrected_heatamp(vsd = vsd, meta_ave = subset_meta,
                               gene_list = genes_plot$`NO_diabetes_antigen - ND_diabetes_antigen`,
                               meta_col = "meta_col", plot_meta_col = FALSE,
                               max_val = 2.5, min_val = -2.5,
                               cluster_rows = TRUE, cluster_cols = FALSE,
                               coloring = color_list, plot_rownames = FALSE))
  
  grid::grid.newpage()
  
  print(make_corrected_heatamp(vsd = vsd, meta_ave = subset_meta2,
                               gene_list = genes_plot$`NO_diabetes_antigen - ND_diabetes_antigen`,
                               meta_col = "meta_col", plot_meta_col = FALSE,
                               max_val = 2.5, min_val = -2.5,
                               cluster_rows = TRUE, cluster_cols = FALSE,
                               coloring = color_list, plot_rownames = FALSE))
  
  dev.off()
  
}

# De genes made in 07_de_analysis.R
vsd_bcells <- readRDS(file.path(save_dir, "rda_obj", 
                                "psuedobulk_batch_corrected_processed_date_object.rds"))
all_res_bcells <- read.csv(file.path(save_dir, "files", 
                                     "pseudobulk_all_bcell_de_processed_date.csv"))

make_de_plots(all_res = all_res_bcells, save_name = "b_cells", 
              vsd = vsd_bcells, seurat_object = seurat_data)


vsd_bcells <- readRDS(file.path(save_dir, "rda_obj", 
                                "psuedobulk_batch_corrected_object.rds"))
all_res_bcells <- read.csv(file.path(save_dir, "files", 
                                     "pseudobulk_all_bcell_de.csv"))

make_de_plots(all_res = all_res_bcells, save_name = "b_cells_collection_bin", 
              vsd = vsd_bcells, seurat_object = seurat_data)

# Plots for 230524 (stacked barplots) ------------------------------------------
# antigen_colors<- c("Negative" = "#7fc97f",
#                    "INS-tet" = "#beaed4",
#                    "TET-tet" = "#fdc086",
#                    "IA2-tet" = "#ffff99",
#                    "GAD-tet" = "#386cb0",
#                    "DNA-tet" = "#f0027f")
# 
# # Keep only cell types of interest and non-doublets
# seurat_sub <- subset(seurat_nd_no, subset = hash.ID != "Doublet" &
#                        final_celltype %in% c("Naive_1",
#                                              "Naive_2",
#                                              "Naive_3",
#                                              "Resting_memory"))
# 
# seurat_sub$mem_naive <- ifelse(seurat_sub$final_celltype %in%
#                                  c("Naive_1", "Naive_2", "Naive_3"), 
#                                "Naive", "Resting_memory")
# 
# seurat_sub$celltype_status <- paste(seurat_sub$mem_naive, seurat_sub$Status,
#                                     sep = "_")
# 
# # Make barplots
# new_barplots <- scAnalysisR::stacked_barplots(seurat_object = seurat_sub, 
#                                               meta_col = "hash.ID",
#                                               split_by = "celltype_status",
#                                               return_values = TRUE,
#                                               color = antigen_colors)
# 
# # Find total of numbers of cell per condition
# total_counts <- new_barplots$data %>%
#   dplyr::group_by(split_by) %>%
#   dplyr::mutate(total_count = paste(sum(Freq), "cells", sep = " ")) %>%
#   dplyr::select(split_by, total_count) %>%
#   dplyr::distinct()
# 
# # Make a handful of plots to pick from
# plot_1 <- new_barplots$barplot +
#   ggplot2::theme(axis.text.x = element_text(angle = 45, hjust=1, size = 10)) +
#   ggplot2::geom_text(
#     data = total_counts, 
#     ggplot2::aes(x = split_by, y = 50, label = total_count),
#     show.legend = F,
#     inherit.aes = F,
#     color = "black", 
#     angle = 90,
#     size = 8
#   ) +
#   ggplot2::labs(x = "", fill = "Antigen") 
# 
# plot_2 <- new_barplots$barplot +
#   ggplot2::theme(axis.text.x = element_text(angle = 45, hjust=1, size = 10)) +
#   ggplot2::geom_text(
#     data = total_counts, 
#     ggplot2::aes(x = split_by, y = 50, label = total_count),
#     show.legend = F,
#     inherit.aes = F,
#     color = "black", 
#     angle = 90,
#     size = 5
#   ) +
#   ggplot2::labs(x = "", fill = "Antigen") 
# 
# plot_3 <- new_barplots$barplot +
#   ggplot2::theme(axis.text.x = element_text(angle = 45, hjust=1, size = 10)) +
#   ggplot2::labs(x = "", fill = "Antigen") 
# 
# pdf(file.path(fig_dir, "celltype_percents_new_colors.pdf"), width = 8,
#     height = 8)
# print(plot_1)
# print(plot_2)
# print(plot_3)
# dev.off()
