library(Seurat)
library(tidyverse)
library(here)
library(scAnalysisR)
library(viridis)
library(clustifyr)
library(singlecellmethods)
library(batchelor)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

cor_cutoff <- 0.5

all_ref_dir <-
  "/beevol/home/wellskri/Analysis/references/single_cell_references"

normalization_method <- "log" # can be SCT or log
# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

args <- commandArgs(trailingOnly = TRUE)

sample <- args[[1]]
sample <- gsub("__.*", "", sample)
#sample <- "merged"

results_dir <- args[[2]]
#results_dir <- here("results")

sample_info <- args[[4]]
#sample_info <- here("files/sample_info.tsv")

sample_info <- read.table(sample_info, fill = TRUE, header = TRUE)

sample_info <- sample_info[sample_info$sample == sample,]
vars.to.regress <- NULL
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
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))

# Remove B cells with more than n = 2 chains -----------------------------------
# The chains we want to remove are:
# IGH;IGH;IGK
# IGH;IGH;IGK;IGK
# IGH;IGH;IGK;IGL
# IGH;IGH;IGL
# IGH;IGH;IGL;IGL
# IGH;IGK;IGK
# IGH;IGL;IGL
cm <- confusionMatrix(seurat_data$chains, seurat_data$Doublet_finder)
cm <- cm / rowSums(cm)

pheatmap::pheatmap(cm)

cm <- confusionMatrix(seurat_data$Doublet_finder, seurat_data$chains)
cm <- cm / rowSums(cm)

pheatmap::pheatmap(cm)

# So remove all the chains with multiple heavy, remove "doublets" from
# all non-plasmablasts

remove_chains <- c("IGH;IGH;IGK", "IGH;IGH;IGK;IGK", "IGH;IGH;IGK;IGL",
                   "IGH;IGH;IGL", "IGH;IGH;IGL;IGL")

remove_cells_chain <- seurat_data[[]] %>%
  dplyr::filter(chains %in% remove_chains)

remove_cells_chain <- rownames(remove_cells_chain)

remove_cells_doublet <- seurat_data[[]] %>%
  dplyr::filter(Doublet_finder == "Doublet" & 
                  RNA_combined_celltype != "Plasmablast")

remove_cells_doublet <- rownames(remove_cells_doublet)

remove_cells <- unique(c(remove_cells_chain, remove_cells_doublet))

keep_cells <- colnames(seurat_data)[!colnames(seurat_data) %in% remove_cells]

seurat_data <- subset(seurat_data, cells = keep_cells)

# Reprocess data ---------------------------------------------------------------
seurat_data <- seurat_data %>%
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = vars.to.regress)

# Update variable features
VariableFeatures(seurat_data) <- VariableFeatures(seurat_data)[!grepl("IG[H|L|K]",
                                                                      VariableFeatures(seurat_data))]

## PCA -------------------------------------------------------------------------

# PCA of gene expression, weighted by cell number per sample
seurat_data <- singlecellmethods::RunBalancedPCA(obj = seurat_data, 
                                                 weight.by = "orig.ident",
                                                 npcs = 50)
## MNN -------------------------------------------------------------------------
sce_data <- as.SingleCellExperiment(seurat_data)

set.seed(0)
# Run fastMNN on batch
corrected_data <- fastMNN(sce_data, batch = sce_data$sample,
                          subset.row = VariableFeatures(seurat_data))

# Check reduced dims
if(!identical(rownames(SingleCellExperiment::reducedDim(x = corrected_data)), 
              colnames(seurat_data))){
  stop("names are not the same between the object and fastmnn output!!!")
}

# Change names
colnames(SingleCellExperiment::reducedDim(x = corrected_data)) <-
  paste0("mnn_", 1:ncol(SingleCellExperiment::reducedDim(x = corrected_data)))

# Add to seurat object
seurat_data$mnn <- CreateDimReducObject(
  embeddings = SingleCellExperiment::reducedDim(x = corrected_data),
  loadings = as.matrix(SingleCellExperiment::rowData(x = corrected_data)),
  assay = DefaultAssay(object = seurat_data),
  key = "mnn_"
)

rm(corrected_data)

## UMAP ------------------------------------------------------------------------

seurat_data$RNA_cluster <- NULL

# Remove previous clustering
remove_cols <- colnames(seurat_data[[]])[grepl("res\\.[0-9]",
                                               colnames(seurat_data[[]]))]

for (i in remove_cols){
  seurat_data[[i]] <- NULL
}

# UMAP of gene expression
set.seed(0)
umap_data <- group_cells(seurat_data, sample, save_dir, nPCs = RNA_pcs,
                         resolution = resolution, assay = seurat_assay,
                         HTO = HTO,
                         reduction = batch_correction)
seurat_data <- umap_data$object

plot_type <- paste0(batch_correction, ".umap")

seurat_data$corrected_cluster <- seurat_data$RNA_cluster


## Name cell types -------------------------------------------------------------

###############
# Seurat pbmc #
###############

# Information for cell mapping
ref_dir <- file.path(all_ref_dir, "pbmc/seurat")

ref_mat <- read.csv(file.path(ref_dir, "clustifyr_reference_l2.csv"),
                    header = TRUE, row.names = 1)

pdf(file.path(save_dir, "images", "no_doublet_celltype_mapping.pdf"))

cluster_res <- name_clusters(seurat_data, ref_mat,
                             save_dir = save_dir,
                             save_name = "comb_celltype_seurat", ADT = FALSE,
                             assay = "RNA",
                             nfeatures = 1000, clusters = "RNA_cluster",
                             plot_type = plot_type,
                             features = VariableFeatures(seurat_data))

seurat_data <- cluster_res$object

seurat_res_seurat <- cluster_res$RNA

grid::grid.newpage()

print(plotDimRed(seurat_data, col_by = "RNA_comb_celltype_seurat",
                 plot_type = plot_type, ggrastr = TRUE))

cm <- confusionMatrix(seurat_data$RNA_comb_celltype_seurat,
                      seurat_data$RNA_celltype_seurat)
cm <- cm / rowSums(cm)

pheatmap::pheatmap(cm)

write.csv(seurat_res_seurat, 
          file.path(save_dir, "files/celltype_mapping_seurat.csv"))



####################
# Published object #
####################

ref_dir <- here("files/210825_object")


ref_mat <- read.csv(file.path(ref_dir, "clustifyr_reference.csv"),
                    header = TRUE, row.names = 1)

DefaultAssay(seurat_data) <- "RNA"

grid::grid.newpage()

cluster_res <- name_clusters(seurat_data, ref_mat,
                             save_dir = save_dir,
                             save_name = "comb_celltype_bnd", ADT = FALSE,
                             assay = "RNA",
                             features = VariableFeatures(seurat_data),
                             clusters = "RNA_cluster",
                             plot_type = plot_type)

seurat_data <- cluster_res$object

seurat_res_bnd <- cluster_res$RNA

grid::grid.newpage()

new_colors <- c("Resting_memory" = "#924bdb", # Resting memory
                "Naive_1" = "#69ba3d", # Naive 1
                "Naive_2" = "#9a43a4", # Naive 2
                "Memory_IgE_IgG" = "#bf9b31", # Memory IgE/IgG1
                "Naive_3" = "#6477ce", # Naive 3
                "Memory_IgA" = "#d15131", # Memory IA
                "Early_memory" = "#4c9e8e", # Early Memory
                "BND2" = "#cc4570", #Bnd2
                "DN2" = "#648d4f", # DN2
                "Activated_memory" = "#985978", # Activated memory
                "Activated_naive" = "#a06846") # Activated naive

print(plotDimRed(seurat_data, col_by = "RNA_comb_celltype_bnd",
                 plot_type = plot_type, color = new_colors,
                 ggrastr = TRUE))

cm <- confusionMatrix(seurat_data$RNA_comb_celltype_bnd,
                      seurat_data$RNA_celltype_bnd)
cm <- cm / rowSums(cm)

pheatmap::pheatmap(cm)


write.csv(seurat_res_bnd,
          file.path(save_dir, "files/celltype_mapping_bnd.csv"))


celltype_colors <- readRDS(file = file.path(save_dir, "color_palette.rds"))

# Merge dfs

seurat_res <- cbind(seurat_res_seurat, seurat_res_bnd)

seurat_cluster <- cor_to_call(seurat_res) %>% 
  mutate(type = ifelse(r < cor_cutoff, "undetermined", type))


new_clusters <- seurat_cluster$type
names(new_clusters) <- seurat_cluster$cluster
seurat_data$RNA_combined_celltype <- new_clusters[seurat_data$RNA_cluster]

print(plotDimRed(seurat_data, col_by = "RNA_combined_celltype",
                 plot_type = plot_type,
                 color = celltype_colors, ggrastr = TRUE))

cm <- confusionMatrix(seurat_data$RNA_combined_celltype,
                      seurat_data$RNA_celltype)
cm <- cm / rowSums(cm)

pheatmap::pheatmap(cm)

dev.off()

saveRDS(seurat_data, file.path(save_dir, "rda_obj", "seurat_processed_no_doublet.rds"))