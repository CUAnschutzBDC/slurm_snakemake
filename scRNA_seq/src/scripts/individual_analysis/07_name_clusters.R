library(Seurat)
library(tidyverse)
library(cowplot)
library(harmony)
library(here)
library(scAnalysisR)
library(viridis)
library(clustifyr)
source(here("src", "scripts", "common_setup.R"))

cor_cutoff <- 0.5

all_ref_dir <-
  "/pl/active/Anschutz_BDC/resources/single_cell_references"

# Read in data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))

res_list_rna <- list()

if(ADT){
  res_list_adt <- list()
  res_list_combined <- list()
}
#-------------------------------------------------------------------------------

# NOTE update the section below to include any references pertinent to your
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

res_list_rna$seurat_rna <- cluster_res$RNA

grid::grid.newpage()

all_celltypes <- c("RNA_celltype_seurat")

if(ADT){
  res_list_adt$seurat_adt <- cluster_res$ADT
  res_list_combined$seurat_wnn <- cluster_res$WNN
  all_celltypes <- c(all_celltypes,
                     "ADT_celltype_seurat",
                     "combined_celltype_seurat")
  
  DefaultAssay(seurat_data) <- "ADT"
  
  all_adts <- rownames(seurat_data)
  
  DefaultAssay(seurat_data) <- "RNA"
}

print(plotDimRed(seurat_data, col_by = all_celltypes,
                 plot_type = "rna.umap"))

mapping_wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb = mapping_wb, sheetName = "seurat_rna")
openxlsx::writeData(wb = mapping_wb, sheet = "seurat_rna",
                    x = res_list_rna$seurat_rna)

if (ADT){
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
  
  openxlsx::addWorksheet(wb = mapping_wb, sheetName = "seurat_adt")
  openxlsx::writeData(wb = mapping_wb, sheet = "seurat_adt",
                      x = res_list_adt$seurat_adt)
  openxlsx::addWorksheet(wb = mapping_wb, sheetName = "seurat_combined")
  openxlsx::writeData(wb = mapping_wb, sheet = "seurat_combined",
                      x = res_list_combined$seurat_wnn)
}


#-------------------------------------------------------------------------------

####################
# Published object #
####################

ref_dir <- file.path(all_ref_dir, "pancreas/Baron_2016")

ref_mat <- read.csv(file.path(ref_dir, "clustifyr_human_reference.csv"),
                    header = TRUE, row.names = 1)

DefaultAssay(seurat_data) <- "RNA"

grid::grid.newpage()

cluster_res <- name_clusters(seurat_data, ref_mat,
                             save_dir = save_dir,
                             save_name = "celltype_baron", ADT = TRUE,
                             assay = "RNA",
                             features = VariableFeatures(seurat_data),
                             clusters = "RNA_cluster",
                             plot_type = "rna.umap")

seurat_data <- cluster_res$object

res_list_rna$baron_rna <- cluster_res$RNA

grid::grid.newpage()

all_celltypes <- c("RNA_celltype_baron")

openxlsx::addWorksheet(wb = mapping_wb, sheetName = "baron_rna")
openxlsx::writeData(wb = mapping_wb, sheet = "baron_rna",
                    x = res_list_rna$baron_rna)

if(ADT){
  res_list_adt$baron_adt <- cluster_res$ADT
  res_list_combined$baron_wnn <- cluster_res$WNN
  all_celltypes <- c(all_celltypes,
                     "ADT_celltype_baron",
                     "combined_celltype_baron")
  
  DefaultAssay(seurat_data) <- "ADT"
  
  all_adts <- rownames(seurat_data)
  
  DefaultAssay(seurat_data) <- "RNA"
  
  openxlsx::addWorksheet(wb = mapping_wb, sheetName = "baron_adt")
  openxlsx::writeData(wb = mapping_wb, sheet = "baron_adt",
                      x = res_list_adt$baron_adt)
  openxlsx::addWorksheet(wb = mapping_wb, sheetName = "baron_combined")
  openxlsx::writeData(wb = mapping_wb, sheet = "baron_combined",
                      x = res_list_combined$baron_wnn)
}

print(plotDimRed(seurat_data, col_by = all_celltypes,
                 plot_type = "rna.umap"))

#-------------------------------------------------------------------------------

# Merge dfs

merge_dfs <- function(res_list, seurat_object, save_name, cluster_col){
  seurat_res <- do.call(cbind, res_list)
  
  seurat_cluster <- cor_to_call(seurat_res) %>% 
    mutate(type = ifelse(r < cor_cutoff, "undetermined", type))
  
  
  new_clusters <- seurat_cluster$type
  names(new_clusters) <- seurat_cluster$cluster
  meta_data <- seurat_object[[]]
  meta_data[[save_name]] <- new_clusters[meta_data[[cluster_col]]]
  meta_data <- meta_data[, save_name, drop = FALSE]
  
  seurat_object <- AddMetaData(seurat_object, metadata = meta_data)
  return(seurat_object)
}

all_celltypes <- c("RNA_celltype")

seurat_data <- merge_dfs(
  seurat_object = seurat_data,
  res_list = res_list_rna,
  save_name = "RNA_celltype",
  cluster_col = "RNA_cluster"
)

if(ADT){
  seurat_data <- merge_dfs(
    seurat_object = seurat_data,
    res_list = res_list_adt,
    save_name = "ADT_celltype",
    cluster_col = "ADT_cluster"
  )
  
  seurat_data <- merge_dfs(
    seurat_object = seurat_data,
    res_list = res_list_combined,
    save_name = "combined_celltype",
    cluster_col = "combined_cluster"
  )  
  
  all_celltypes <- c(all_celltypes,
                     "ADT_celltype",
                     "combined_celltype")
}



print(plotDimRed(seurat_data, col_by = all_celltypes,
                 plot_type = "rna.umap"))


if(ADT){
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
  
}

dev.off()

openxlsx::saveWorkbook(wb = mapping_wb,
                       file = file.path(save_dir, "files/celltype_mapping.xlsx"),
                       overwrite = TRUE)

saveRDS(seurat_data, file.path(save_dir, "rda_obj", "seurat_processed.rds"))
