library(Seurat)
library(tidyverse)
library(cowplot)
library(openxlsx)
library(here)
library(scAnalysisR)


# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

source(here("src/scripts/functions.R"))

remove_ambience <- FALSE

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
#sample <- "pancreas"

sample_info <- args[[3]]
#sample_info <- here("files/sample_info.tsv")

results_dir <- args[[2]]
#results_dir <- here("results")

sample_info <- read.table(sample_info, fill = TRUE, header = TRUE)

samples_use <- sample_info[sample_info$sample != sample, ]$sample

sample_info <- sample_info[sample_info$sample == sample,]

HTO <- sample_info$HTO
ADT <- sample_info$ADT
hash_ident <- sample_info$hash_ident

cell_types <- "RNA_combined_celltype"
clusters <- "RNA_cluster"
pval <- 0.05
logfc <- 0.5

# ADT_pca <- sample_info$adt_pca
# run_adt_umap <- sample_info$adt_umap

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

mapping_file <- read.csv(here("files/species_mapping_file.csv"))

gene_lists <- read.table(file.path("/beevol", "home", "wellskri", "Analysis", 
                                   "references", "single_cell_references",
                                   "gene_lists", 
                                   "PanglaoDB_markers_27_Mar_2020.tsv.gz"),
                         sep = "\t", header = TRUE)

# Cell type DE -----------------------------------------------------------------


# RNA cluster DE ---------------------------------------------------------------

seurat_data$cluster_celltype <- paste0(seurat_data[[clusters]][[1]], "_",
                                       seurat_data[[cell_types]][[1]])

marker_list <- find_write_markers_orthologs(seurat_object = seurat_data,
                                            meta_col = "cluster_celltype",
                                            pval = pval,
                                            logfc = logfc,
                                            assay = "RNA",
                                            save_dir = save_dir,
                                            mapping_file = mapping_file,
                                            mapping_gene_col = "gene_id",
                                            mapping_ortholog_col = c("Mouse.gene.name",
                                                                     "Human.gene.name",
                                                                     "Dog.gene.name",
                                                                     "Pig.gene.name"),
                                            gene_lists = NULL)

write.csv(marker_list, file.path(save_dir, "files", "all_markers.csv"))


gene_list <- c("INS", "SST", "GCG", "MKI67", "CFTR")

plots <- featDistPlot(seurat_data, geneset = gene_list, sep_by = "RNA_cluster",
                      combine = FALSE)

pdf(file.path(save_dir, "images", "marker_gene_violin.pdf"),
    height = 6, width = 15)
print(plots)
dev.off()


library(enrichR)
websiteLive <- getOption("enrichR.live")
if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes   
}

if (websiteLive) dbs <- listEnrichrDbs()

dbs_use <- c("Human_Gene_Atlas", "Mouse_Gene_Atlas",
             "PanglaoDB_Augmented_2021", "Tabula_Muris",
             "Tabula_Sapiens")

wb <- openxlsx::createWorkbook()

human <- lapply(unique(marker_list$cluster), function(x){
  markers_use <- marker_list %>%
    dplyr::filter(cluster == x, p_val_adj < 0.05)
  
  enriched <- enrichr(unique(markers_use$Human.gene.name), dbs_use)
  
  top_enriched <- lapply(enriched, function(enriched_df){
    enriched_df <- enriched_df %>%
      dplyr::arrange(desc(Combined.Score))
    
    top_term <- enriched_df[1,]
  })
  top_enriched <- do.call(rbind, top_enriched)
  
  top_enriched$database <- dbs_use
  
  x <- as.character(x)
  # Save as excel file
  openxlsx::addWorksheet(wb = wb, sheetName = x)
  openxlsx::writeData(wb = wb, sheet = x, x = top_enriched)
  
  return(top_enriched)
})

openxlsx::saveWorkbook(wb = wb, file = file.path(save_dir, "files",
                                                 "enrichr_human_gene.xlsx"),
                       overwrite = TRUE)

# Try with module score --------------------------------------------------------
celltypes_use <- c("Acinar cells", "Alpha cells", "B cells",
                   "Beta cells", "Chromaffin cells", "Delta cells",
                   "Ductal cells", "Endothelial cells", "Gamma (PP) cells",
                   "Macrophages", "Monocytes")

all_plots <- lapply(celltypes_use, function(x){
  celltype_genes <- gene_lists %>%
    dplyr::filter(cell.type == x)
  
  all_genes <- celltype_genes$official.gene.symbol
  
  mapping_genes <- mapping_file %>%
    dplyr::filter(Human.gene.name %in% all_genes) %>%
    dplyr::filter(!is.na(Human.gene.name)) %>%
    dplyr::filter(Human.gene.name != "")
  
  human_genes <- unique(mapping_genes$gene_id)
  
  ferret_genes <- setdiff(all_genes, unique(mapping_genes$Human.gene.name))
  
  x <- make.names(x)
  
  so <- AddModuleScore(seurat_data,
                       features = list(c(human_genes, ferret_genes)),
                       name = x)
  
  return_plot <- featDistPlot(so, geneset = paste0(x, "1"),
                              sep_by = "cluster_celltype",
                              combine = FALSE, col_by = "cluster_celltype")
  
  return(return_plot)
})

pdf(file.path(save_dir, "images", "panglao_db_module_violin.pdf"),
    height = 6, width = 15)
print(all_plots)
dev.off()

marker_genes <- c("INS", # Beta
                  "GCG", # Alpha
                  "SST", # Delta
                  "NKX2-2", "NKX6-1", "PAX6", # Endocrine
                  "CFTR", # ductal
                  "PDGFRA", # fibroblast
                  "CD3D", # T cells
                  "ITGAM", # Myeloid
                  "SOX2", "NOTCH1", # Neurons
                  "ALB", # Hepatocyte
                  "CD19", # B cell
                  "IL2RA", # T cell
                  "NES", # Stellate cell
                  "CDH5", "PECAM1", # Endothelail
                  "OLFM4", # epithelail
                  "FEV" # fev progenitor
)

marker_plot <- featDistPlot(seurat_data, geneset = marker_genes,
                            sep_by = "cluster_celltype",
                            combine = FALSE, col_by = "cluster_celltype")

pdf(file.path(save_dir, "images", "RNA_markers_violin.pdf"),
    height = 6, width = 15)
print(marker_plot)
dev.off()




rename_celltype <- c("0" = "acinar", 
                     "1" = "ductal", 
                     "2" = "ductal", 
                     "3" = "ductal", 
                     "4" = "beta", 
                     "5" = "ductal", 
                     "6" = "acinar", 
                     "7" = "ductal", 
                     "8" = "fev_endocrine_progenitor",
                     "9" = "acinar", 
                     "10" = "delta", 
                     "11" = "ductal", 
                     "12" = "ductal", 
                     "13" = "fev_endocrine_progenitor", 
                     "14" = "beta", 
                     "15" = "ductal", 
                     "16" = "alpha", 
                     "17" = "ductal", 
                     "18" = "acinar", 
                     "19" = "ductal", 
                     "20" = "ductal", 
                     "21" = "ductal", 
                     "22" = "beta", 
                     "23" = "ductal") 


seurat_data$final_combined_celltype <- rename_celltype[seurat_data$RNA_cluster]

all_plots <- lapply(unique(seurat_data$final_combined_celltype), function(x){
  plotDimRed(seurat_data, col_by = "RNA_cluster", plot_type = "rna.umap",
             highlight_group = TRUE, group = x, meta_data_col = "final_combined_celltype",
             ggrastr = TRUE)
})

pdf(file.path(save_dir, "images", "UMAP_celltype_single.pdf"),
    height = 6, width = 8)
print(all_plots)
dev.off()

seurat_data$cluster_celltype <- paste0(seurat_data[[clusters]][[1]], "_",
                                       seurat_data$final_combined_celltype)

marker_plot <- featDistPlot(seurat_data, geneset = marker_genes,
                            sep_by = "cluster_celltype",
                            combine = FALSE, col_by = "cluster_celltype")

pdf(file.path(save_dir, "images", "RNA_markers_violin_final.pdf"),
    height = 6, width = 15)
print(marker_plot)
dev.off()


pdf(file.path(save_dir, "images", "final_celltype_umap.pdf"))
print(plotDimRed(seurat_data, col_by = "final_combined_celltype", plot_type = "rna.umap"))
dev.off()


saveRDS(seurat_data, file.path(save_dir, "rda_obj", "seurat_processed.rds"))
