library(ggplot2)
library(Seurat)
library(tidyr)
library(cowplot)
library(dplyr)
library(pathview)
library(openxlsx)
library(gprofiler2)
library(scAnalysisR)
library(here)

source(here("src/scripts/common_setup.R"))

# Update this with any pathways you would like to visiualize
path_id_list <- c(NFKB_sig_path = "04064", T_cell_receptor = "04660",
                  cytokine_receptor_interaction = "04060",
                  cell_adhesion = "04514", Hematopoitic_lineage = "04640",
                  Type1_diabetes = "04940", NK_cytoxicity = "04650",
                  th1_th2_differentiation = "04658", TNF_signaling = "04668",
                  immunodeficiency = "05340", Th17_differentiation = "04659")

gen_id <- "hsa"

# Set to TRUE if you needed to set cell types manually, also
# name the manual cell types "final_celltype"
manual_cell_types <- FALSE

# First, determine the clusters
if (ADT & run_adt_umap) {
  clusters <- "combined_cluster"
} else {
  clusters <- "RNA_cluster"
}

# Then, determine the cell types
if (manual_cell_types) {
  cell_types <- "final_celltype"
} else if (ADT & run_adt_umap) {
  cell_types <- "combined_celltype"
} else {
  cell_types <- "RNA_celltype"
}

pval <- 0.05
logfc <- 0.5

# Read in data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))

save_dir <- here("results", "R_analysis", sample)

out_dir <- file.path(save_dir, "images", "pathways")

out_dir_go <- file.path(save_dir, "images", "GSEA")
out_dir_text <- file.path(save_dir, "files", "GSEA")

# Create output directory
out_dir %>%
  dir.create(showWarnings = F)

out_dir_go %>%
  dir.create(showWarnings = F)

out_dir_text %>%
  dir.create(showWarnings = F)


# Read in the data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))


# Cell type DE -----------------------------------------------------------------

# Find markers
marker_genes <- find_markers(seurat_object = seurat_data,
                             seurat_assay = seurat_assay,
                             test_idents = cell_types)


# Run gost
gost_output <- run_gost(seurat_de = marker_genes,
                        sources = c("GO:BP", "KEGG", "REAC", "TF"))

# Save gost output
save_gost(gost_output,
          save_dir_plots = file.path(out_dir_go, "celltype"),
          save_dir_text = file.path(out_dir_text, "celltype"))



de_to_pathview(seurat_de = marker_genes,
               path_id_list = path_id_list,
               out_dir = file.path(out_dir, "celltype"),
               seurat_assay = seurat_assay,
               test_idents = cell_types,
               gen_id = gen_id)

# cluster DE -------------------------------------------------------------------

# Find markers
marker_genes <- find_markers(seurat_object = seurat_data,
                             seurat_assay = seurat_assay,
                             test_idents = clusters)


# Run gost
gost_output <- run_gost(seurat_de = marker_genes,
                        sources = c("GO:BP", "KEGG", "REAC", "TF"))

# Save gost output
save_gost(gost_output,
          save_dir_plots = file.path(out_dir_go, "cluster"),
          save_dir_text = file.path(out_dir_text, "cluster"))


de_to_pathview(seurat_de = marker_genes,
               path_id_list = path_id_list,
               out_dir = file.path(out_dir, "cluster"),
               seurat_assay = seurat_assay,
               test_idents = clusters,
               gen_id = gen_id)

