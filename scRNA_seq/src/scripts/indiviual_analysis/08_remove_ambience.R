library(Seurat)
library(tidyverse)
library(cowplot)
library(here)
library(scAnalysisR)
library(DropletUtils)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log


cell_types <- "RNA_celltype"
clusters <- "RNA_cluster"
pval <- 0.05
logfc <- 0.5

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

args <- commandArgs(trailingOnly = TRUE)

sample <- args[[1]]
sample <- gsub("__.*", "", sample)
#sample <- "JH310-12_AP"

results_dir <- args[[2]]
#results_dir <- here("results")

sample_info <- args[[4]]
#sample_info <- here("files/sample_info.tsv")

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

# Read in data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))

DefaultAssay(seurat_data) <- "RNA"
# Background contamination -----------------------------------------------------
raw_path <- file.path(results_dir, sample, "outs", "multi",
                      "count", "raw_feature_bc_matrix")

raw_sce <- DropletUtils::read10xCounts(raw_path, col.names = TRUE)

# Remove the antibodies
raw_sce <- raw_sce[grepl("ENSG", rownames(raw_sce)), ]

# Rename genes based on gene name not ens id
gene_mapping <- read.table(file.path(raw_path, "features.tsv.gz"),
                           sep = "\t")

colnames(gene_mapping) <- c("ens_id", "gene_id", "type")

# Seurat deals with gene name duplicates with the "make names" argument
gene_mapping <- gene_mapping %>%
  dplyr::filter(type == "Gene Expression") %>%
  dplyr::mutate(gene_id = make.unique(gene_id)) %>%
  dplyr::filter(gene_id %in% rownames(seurat_data)) 

# Subset to the cells and genes kept by seurat
raw_sce <- raw_sce[rownames(raw_sce) %in% gene_mapping$ens_id, ]
filtered_sce <- raw_sce[ , colnames(raw_sce) %in% colnames(seurat_data)]

if(!identical(colnames(filtered_sce), colnames(seurat_data))){
  filtered_sce <- filtered_sce[ , order(match(colnames(filtered_sce),
                                              colnames(seurat_data)))]
}

colLabels(filtered_sce) <- seurat_data[[clusters]][[1]]

set.seed(100)
e.out <- emptyDrops(counts(raw_sce))

amb <- metadata(e.out)$ambient[,1]
stripped <- filtered_sce[names(amb),]

out <- removeAmbience(counts(stripped), ambient=amb, groups=colLabels(stripped))

# Make sure something happened
differences <- out != counts(stripped)
sums <- rowSums(differences)
changed <- sums[sums > 0]

counts(stripped, withDimnames=FALSE) <- out

final_counts <- counts(stripped)

gene_mapping <- gene_mapping %>%
  dplyr::filter(ens_id %in% rownames(final_counts))

changed_ids <- gene_mapping %>%
  dplyr::filter(ens_id %in% names(changed))

changed_ids$count <- changed[changed_ids$ens_id]

# Now order this to be the same as the counts
final_counts <- final_counts[rownames(final_counts) %in% gene_mapping$ens_id, ]

final_counts <- final_counts[order(match(rownames(final_counts),
                                         gene_mapping$ens_id)), ]

rownames(final_counts) <- gene_mapping$gene_id

seurat_data[["AMBRNA"]] <- CreateAssayObject(counts = final_counts)

seurat_data <- NormalizeData(seurat_data, assay = "AMBRNA")

saveRDS(seurat_data, file.path(save_dir, "rda_obj", "seurat_processed.rds"))

# 
# avg_expr <- AverageExpression(seurat_data, group.by = "RNA_celltype")
# 
# avg_expr_first <- data.frame(avg_expr$RNA)
# avg_expr_second <- data.frame(avg_expr$AMBRNA)
# 
# correlation <- cor(avg_expr_first$Activated_memory, avg_expr_second$Activated_memory)
# 
# new_plot <- data.frame("uncorrected" = avg_expr_first$CD14.Mono,
#                        "corrected" = avg_expr_second$CD14.Mono)
# rownames(new_plot) <- rownames(avg_expr_first)
# 
# ggplot2::ggplot(new_plot, ggplot2::aes(x = uncorrected, y = corrected)) +
#   ggplot2::geom_point()
# 
# 
# new_plot <- new_plot %>%
#   dplyr::arrange(desc(uncorrected))
# 
# differences2 <- GetAssayData(seurat_data, assay = "AMBRNA", slot = "counts") != 
#   GetAssayData(seurat_data, assay = "RNA", slot = "counts")
# sums <- rowSums(differences2)
# changed <- sums[sums > 0]
# 
# 
# group1 <- GetAssayData(seurat_data, assay = "AMBRNA", slot = "counts")["FP236383.3",]
# 
# group2 <- GetAssayData(seurat_data, assay = "RNA", slot = "counts")["FP236383.3",]
# 
# 
# p1 <- featDistPlot(seurat_data, geneset = "FP236383.3", sep_by = clusters, combine = FALSE,
#                    assay = "RNA")[[1]]
# 
# p2 <- featDistPlot(seurat_data, geneset = "FP236383.3", sep_by = clusters, combine = FALSE,
#                    assay = "AMBRNA")[[1]]
# 
# cowplot::plot_grid(p1, p2)
# 
# 
# p1 <- featDistPlot(seurat_data, geneset = "IGKV1-16", sep_by = clusters, combine = FALSE,
#                    assay = "RNA")[[1]]
# 
# p2 <- featDistPlot(seurat_data, geneset = "IGKV1-16", sep_by = clusters, combine = FALSE,
#                    assay = "AMBRNA")[[1]]
# 
# cowplot::plot_grid(p1, p2)
# 
# 
# p1 <- featDistPlot(seurat_data, geneset = "MT-CO1", sep_by = clusters, combine = FALSE,
#                    assay = "RNA")[[1]]
# 
# p2 <- featDistPlot(seurat_data, geneset = "MT-CO1", sep_by = clusters, combine = FALSE,
#                    assay = "AMBRNA")[[1]]
# 
# cowplot::plot_grid(p1, p2)
