library(Seurat)
library(tidyverse)
library(cowplot)
library(harmony)
library(here)
library(scAnalysisR)
library(viridis)
library(clustifyr)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

cor_cutoff <- 0.5

all_ref_dir <-
  "/pl/active/Anschutz_BDC/resources/single_cell_references"

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

args <- commandArgs(trailingOnly = TRUE)

sample <- args[[1]]
sample <- gsub("__.*", "", sample)
#sample <- "Npod6456_PLN"

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
save_dir <- file.path(results_dir, "R_analysis", sample)

# Read in data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))

# Work around because I'm too lazy to update my code
seurat_data@reductions$adt.umap <- seurat_data@reductions$pdsb.umap
#-------------------------------------------------------------------------------

# NOTE update the section below to include any references pertanent to your
# sample.

################
# Seurat pbmc #
###############

# Information for cell mapping
ref_dir <- file.path(all_ref_dir, "pbmc/seurat")

ref_mat <- read.csv(file.path(ref_dir, "clustifyr_reference_l2.csv"),
                    header = TRUE, row.names = 1)

pdf(file.path(save_dir, "images", "celltype_mapping.pdf"))

cluster_res <- name_clusters(seurat_data, ref_mat,
                             save_dir = save_dir,
                             save_name = "celltype_seurat", ADT = TRUE,
                             assay = "RNA",
                             nfeatures = 1000, clusters = "RNA_cluster",
                             plot_type = "rna.umap",
                             features = VariableFeatures(seurat_data))

seurat_data <- cluster_res$object

seurat_res_seurat_rna <- cluster_res$RNA
seurat_res_seurat_adt <- cluster_res$ADT
seurat_res_seurat_wnn <- cluster_res$WNN


grid::grid.newpage()

all_celltypes <- c("RNA_celltype_seurat",
                   "ADT_celltype_seurat",
                   "combined_celltype_seurat")

print(plotDimRed(seurat_data, col_by = all_celltypes,
                 plot_type = "rna.umap"))
# 
# cm <- confusionMatrix(seurat_data$RNA_celltype_seurat,
#                       seurat_data$ADT_celltype_seurat)
# 
# cm <- cm /rowSums(cm)
# pheatmap::pheatmap(cm)

DefaultAssay(seurat_data) <- "ADT"

all_adts <- rownames(seurat_data)

DefaultAssay(seurat_data) <- "RNA"


all_violins <- featDistPlot(seurat_data, geneset = all_adts,
                            col_by = "RNA_celltype_seurat",
                            sep_by = "RNA_celltype_seurat",
                            combine = FALSE)

print(all_violins)

all_violins <- featDistPlot(seurat_data, geneset = all_adts,
                            col_by = "ADT_celltype_seurat",
                            sep_by = "ADT_celltype_seurat",
                            combine = FALSE)

print(all_violins)

all_violins <- featDistPlot(seurat_data, geneset = all_adts,
                            col_by = "combined_celltype_seurat",
                            sep_by = "combined_celltype_seurat",
                            combine = FALSE)

print(all_violins)


mapping_wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb = mapping_wb, sheetName = "seurat_rna")
openxlsx::writeData(wb = mapping_wb, sheet = "seurat_rna",
                    x = seurat_res_seurat_rna)
openxlsx::addWorksheet(wb = mapping_wb, sheetName = "seurat_adt")
openxlsx::writeData(wb = mapping_wb, sheet = "seurat_adt",
                    x = seurat_res_seurat_adt)
openxlsx::addWorksheet(wb = mapping_wb, sheetName = "seurat_combined")
openxlsx::writeData(wb = mapping_wb, sheet = "seurat_combined",
                    x = seurat_res_seurat_wnn)



#-------------------------------------------------------------------------------

####################
# Published object #
####################

ref_dir <- "/beevol/home/wellskri/Analysis/Mia_Smith/Catherine_Nicolas/full_antigen_pos_data/files/210825_object"


ref_mat <- read.csv(file.path(ref_dir, "clustifyr_reference.csv"),
                    header = TRUE, row.names = 1)

DefaultAssay(seurat_data) <- "RNA"

grid::grid.newpage()

cluster_res <- name_clusters(seurat_data, ref_mat,
                             save_dir = save_dir,
                             save_name = "celltype_bnd", ADT = TRUE,
                             assay = "RNA",
                             features = VariableFeatures(seurat_data),
                             clusters = "RNA_cluster",
                             plot_type = "rna.umap")

seurat_data <- cluster_res$object

seurat_res_bnd_rna <- cluster_res$RNA
seurat_res_bnd_adt <- cluster_res$ADT
seurat_res_bnd_wnn <- cluster_res$WNN


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

all_celltypes <- c("RNA_celltype_bnd",
                   "ADT_celltype_bnd",
                   "combined_celltype_bnd")

print(plotDimRed(seurat_data, col_by = all_celltypes,
                 plot_type = "rna.umap", color = new_colors))

openxlsx::addWorksheet(wb = mapping_wb, sheetName = "bnd_rna")
openxlsx::writeData(wb = mapping_wb, sheet = "bnd_rna",
                    x = seurat_res_bnd_rna)
openxlsx::addWorksheet(wb = mapping_wb, sheetName = "bnd_adt")
openxlsx::writeData(wb = mapping_wb, sheet = "bnd_adt",
                    x = seurat_res_bnd_adt)
openxlsx::addWorksheet(wb = mapping_wb, sheetName = "bnd_combined")
openxlsx::writeData(wb = mapping_wb, sheet = "bnd_combined",
                    x = seurat_res_bnd_wnn)

#-------------------------------------------------------------------------------

# Merge dfs

merge_dfs <- function(res1, res2, seurat_object, save_name, cluster_col){
  seurat_res <- cbind(res1, res2)
  
  seurat_cluster <- cor_to_call(seurat_res) %>% 
    mutate(type = ifelse(r < cor_cutoff, "undetermined", type))
  
  
  new_clusters <- seurat_cluster$type
  names(new_clusters) <- seurat_cluster$cluster
  seurat_object[[save_name]] <- new_clusters[seurat_data[[cluster_col]][[1]]]
  return(seurat_object)
}

seurat_data <- merge_dfs(seurat_object = seurat_data,
                         res1 = seurat_res_seurat_rna,
                         res2 = seurat_res_bnd_rna,
                         save_name = "RNA_celltype",
                         cluster_col = "RNA_cluster")

seurat_data <- merge_dfs(seurat_object = seurat_data,
                         res1 = seurat_res_seurat_adt,
                         res2 = seurat_res_bnd_adt,
                         save_name = "ADT_celltype",
                         cluster_col = "ADT_cluster")

seurat_data <- merge_dfs(seurat_object = seurat_data,
                         res1 = seurat_res_seurat_wnn,
                         res2 = seurat_res_bnd_wnn,
                         save_name = "combined_celltype",
                         cluster_col = "combined_cluster")

all_celltypes <- c("RNA_celltype",
                   "ADT_celltype",
                   "combined_celltype")

print(plotDimRed(seurat_data, col_by = all_celltypes,
                 plot_type = "rna.umap"))



all_violins <- featDistPlot(seurat_data, geneset = all_adts,
                            col_by = "RNA_celltype",
                            sep_by = "RNA_celltype",
                            combine = FALSE)

print(all_violins)

all_violins <- featDistPlot(seurat_data, geneset = all_adts,
                            col_by = "ADT_celltype",
                            sep_by = "ADT_celltype",
                            combine = FALSE)

print(all_violins)

all_violins <- featDistPlot(seurat_data, geneset = all_adts,
                            col_by = "combined_celltype",
                            sep_by = "combined_celltype",
                            combine = FALSE)

print(all_violins)

all_violins <- featDistPlot(seurat_data, geneset = all_adts,
                            col_by = "RNA_celltype",
                            sep_by = "RNA_cluster",
                            combine = FALSE)

print(all_violins)

all_violins <- featDistPlot(seurat_data, geneset = all_adts,
                            col_by = "ADT_celltype",
                            sep_by = "ADT_cluster",
                            combine = FALSE)

print(all_violins)

all_violins <- featDistPlot(seurat_data, geneset = all_adts,
                            col_by = "combined_celltype",
                            sep_by = "combined_cluster",
                            combine = FALSE)

print(all_violins)

dev.off()

openxlsx::saveWorkbook(wb = mapping_wb,
                       file = file.path(save_dir, "files/celltype_mapping.xlsx"),
                       overwrite = TRUE)

saveRDS(seurat_data, file.path(save_dir, "rda_obj", "seurat_processed.rds"))