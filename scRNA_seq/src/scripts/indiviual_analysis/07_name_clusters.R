library(Seurat)
library(tidyverse)
library(cowplot)
library(harmony)
library(here)
library(scAnalysisR)
library(viridis)
library(clustifyr)

source(here("src/scripts/common_setup.R"))

# Read in data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))

reference_configs <- list(
  list(
    ref_dir = file.path(all_ref_dir, "pbmc/seurat"),
    ref_file = "clustifyr_reference_l2.csv",
    save_name = "celltype_seurat"
  ),
  list(
    ref_dir = "/beevol/home/wellskri/Analysis/Mia_Smith/Catherine_Nicolas/full_antigen_pos_data/files/210825_object",
    ref_file = "clustifyr_reference.csv",
    save_name = "celltype_bnd"
  )
  # Add more reference configurations as needed
)

if (ADT){
  # Work around because I'm too lazy to update my code
  seurat_data@reductions$adt.umap <- seurat_data@reductions$pdsb.umap
}

cluster_data <- function(seurat_data, ref_mat, save_dir, save_name,
                         ADT, run_adt_umap, mapping_wb, VDJ = FALSE){

  if (VDJ) {
    features <- VariableFeatures(seurat_data)
  } else {
    features <- NULL
  }
  # Name the clusters
  cluster_res <- name_clusters(seurat_data, ref_mat,
                               save_dir = save_dir,
                               save_name = save_name, ADT = ADT,
                               assay = "RNA",
                               nfeatures = 1000, clusters = "RNA_cluster",
                               plot_type = "rna.umap",
                               features = features)

  # Pull out the data
  seurat_data <- cluster_res$object
  res_rna <- cluster_res$RNA

  return_list <- list(object = seurat_data, rna = res_rna)

  # Figure out the meta data column
  all_celltypes <- c(paste0("RNA_", save_name))

  # write to excel
  write_excel(wb = mapping_wb, sheet_name = paste0(save_name, "_rna"), data = res_rna)

  # Repeat for ADT
  if (ADT & run_adt_umap) {

    # Pull out cell types determined based on ADT and WNN clustering
    res_adt <- cluster_res$ADT
    res_wnn <- cluster_res$WNN

    return_list <- c(return_list, list(adt = res_adt, wnn = res_wnn))

    plot_violins(save_name, seurat_data)

    # Add ADT data to the excel file
    write_excel(wb = mapping_wb, sheet_name = paste0(save_name, "_adt"), data = res_adt)
    write_excel(wb = mapping_wb, sheet_name = paste0(save_name, "_combined"), data = res_wnn)
  }

  grid::grid.newpage()

  print(plotDimRed(seurat_data, col_by = all_celltypes,
                  plot_type = "rna.umap"))

  return(return_list)

}

plot_violins <- function(save_name, seurat_data) {
  # Find meta data columns
  all_celltypes <- c(all_celltypes,
                    c(paste0("ADT_", save_name),
                      paste0("combined_", save_name)))

  # Find all ADTs
  DefaultAssay(seurat_data) <- "ADT"

  all_adts <- rownames(seurat_data)

  DefaultAssay(seurat_data) <- "RNA"

  # Plot ADTs across all types of the cell type calls
  all_violins <- featDistPlot(seurat_data, geneset = all_adts,
                              col_by = paste("RNA_", save_name),
                              sep_by = paste("RNA_", save_name),
                              combine = FALSE)

  print(all_violins)

  all_violins <- featDistPlot(seurat_data, geneset = all_adts,
                              col_by = paste("ADT_", save_name),
                              sep_by = paste("ADT_", save_name),
                              combine = FALSE)

  print(all_violins)

  all_violins <- featDistPlot(seurat_data, geneset = all_adts,
                              col_by = paste("combined_", save_name),
                              sep_by = paste("combined_", save_name),
                              combine = FALSE)
  print(all_violins)

}

write_excel <- function(wb, sheet_name, data) {
  openxlsx::addWorksheet(wb = wb, sheetName = sheet_name)
  openxlsx::writeData(wb = wb, sheet = sheet_name,
                      x = data)
}

process_all_references <- function(seurat_data, reference_configs, save_dir,
                                   ADT, run_adt_umap, VDJ, mapping_wb) {
  results <- list()
  
  for (config in reference_configs) {
    ref_mat <- read.csv(file.path(config$ref_dir, config$ref_file), header = TRUE, row.names = 1)
    
    res <- cluster_data(seurat_data = seurat_data, 
                        ref_mat = ref_mat, 
                        save_dir = save_dir,
                        save_name = config$save_name, 
                        ADT = ADT, 
                        run_adt_umap = run_adt_umap,
                        VDJ = VDJ, 
                        mapping_wb = mapping_wb)

    seurat_data <- res$object
    res$object <- NULL
    results[[config$save_name]] <- res
  }
  
  return(list(results = results, final_object = seurat_data))
}

# References ----------------------------------------------------------------
# Create a workbook
mapping_wb <- openxlsx::createWorkbook()

# Open a PDF
pdf(file.path(save_dir, "images", "celltype_mapping.pdf"))

if (VDJ_B | VDJ_T) {
  VDJ <- TRUE
} else {
  VDJ <- FALSE
}

all_results <- process_all_references(seurat_data, reference_configs, save_dir,
                                      ADT, run_adt_umap, VDJ, mapping_wb)

seurat_data <- all_results$final_object
## merged ---------------------------------------------------------------------------

# Merge dfs

merge_dfs <- function(res_list, seurat_object, save_name, cluster_col){
  full_res <- do.call(cbind, res_list)
  
  full_celltype <- cor_to_call(full_res) %>% 
    mutate(type = ifelse(r < cor_cutoff, "undetermined", type))
  
  
  new_clusters <- full_celltype$type
  names(new_clusters) <- full_celltype$cluster
  seurat_object[[save_name]] <- new_clusters[seurat_data[[cluster_col]][[1]]]
  return(seurat_object)
}

merge_rna_results <- function(all_results, seurat_object) {
  rna_results <- lapply(all_results$results, function(x) x$rna)
  merged_rna <- merge_dfs(rna_results, seurat_object, "RNA_celltype", "RNA_cluster")
  return(merged_rna)
}

merge_adt_results <- function(all_results, seurat_object) {
  adt_results <- lapply(all_results$results, function(x) x$adt)
  merged_adt <- merge_dfs(adt_results, seurat_object, "ADT_celltype", "ADT_cluster")
  return(merged_adt)
}

merge_wnn_results <- function(all_results, seurat_object) {
  wnn_results <- lapply(all_results$results, function(x) x$wnn)
  merged_wnn <- merge_dfs(wnn_results, seurat_object, "combined_celltype", "combined_cluster")
  return(merged_wnn)
}

final_seurat_object <- merge_rna_results(all_results, final_seurat_object)
if (ADT) {
  final_seurat_object <- merge_adt_results(all_results, final_seurat_object)
  final_seurat_object <- merge_wnn_results(all_results, final_seurat_object)
}

plot_violins(save_name = "celltype", seurat_data = seurat_data)

dev.off()

openxlsx::saveWorkbook(wb = mapping_wb,
                       file = file.path(save_dir, "files/celltype_mapping.xlsx"),
                       overwrite = TRUE)

saveRDS(seurat_data, file.path(save_dir, "rda_obj", "seurat_processed.rds"))