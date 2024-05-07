library(Seurat)
library(tidyverse)
library(cowplot)
library(harmony)
library(here)
library(scAnalysisR)
library(clustree)

# Look at comment in github for why I made this decision
remove_df_doublets <- TRUE

vars.to.regress <- NULL


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

adt_PCA <- sample_info$adt_pca
run_adt_umap <- sample_info$adt_umap


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
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_adt.rds"))

if(remove_df_doublets){
  # Remove doublet finder doublets
  Idents(seurat_data) <- "Doublet_finder"
  seurat_data <- subset(x = seurat_data, idents = "Singlet")  
}


# Remove HTO doublets
if(HTO){
  Idents(seurat_data) <- "HTO_classification.global"
  seurat_data <- subset(x = seurat_data, idents = "Singlet")
}

seurat_data <- seurat_data %>%
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = vars.to.regress)

# Update variable features
# Remove B cell chains
VariableFeatures(seurat_data) <- VariableFeatures(seurat_data)[!grepl("IG[H|L|K]",
                                                                      VariableFeatures(seurat_data))]

# Remove T cell chains
VariableFeatures(seurat_data) <- VariableFeatures(seurat_data)[!grepl("TR[A|B][V|J|C]",
                                                                      VariableFeatures(seurat_data))]



# PCA --------------------------------------------------------------------------

# PCA of gene expression
seurat_data <- PCA_dimRed(seurat_data, assay = seurat_assay)

RNA_plots <- plot_PCA(HTO = HTO, assay = seurat_assay,
                      sample_object = seurat_data)

pdf(file.path(save_dir, "images", "RNA_pca.pdf"))
print(RNA_plots)
dev.off()

if(ADT & run_adt_umap){
  if(adt_PCA){
    # PCA of surface protein
    seurat_data <- PCA_dimRed(seurat_data, assay = "ADT")
    
    ADT_plots <- plot_PCA(HTO = HTO, assay = "ADT", sample_object = seurat_data)
    
  } else {
    # set up dsb values to use in WNN analysis 
    DefaultAssay(seurat_data) <- "ADT"
    
    # hack seurat to use normalized protein values as a dimensionality reduction object.
    VariableFeatures(seurat_data) <- rownames(seurat_data)
    
    # run true pca to initialize dr pca slot for WNN 
    seurat_data <- ScaleData(seurat_data, verbose = FALSE) %>%
      RunPCA(reduction.name = "pdsb",
             features = VariableFeatures(seurat_data),
             verbose = FALSE)
    
    # make matrix of norm values to add as dr embeddings
    pseudo <- t(GetAssayData(seurat_data, slot = "data"))
    pseudo_colnames <- paste('pseudo', 1:ncol(pseudo), sep = "_")
    colnames(pseudo) <- pseudo_colnames
    # add to object 
    seurat_data@reductions$pdsb@cell.embeddings = pseudo
    
    ADT_plots <- plotDimRed(seurat_data,
                            col_by = c("orig.ident", "percent.mt",
                                       "nFeature_RNA", "nCount_RNA",
                                       "nFeature_ADT", "nCount_ADT"),
                            plot_type = "pdsb")
  }
  
  pdf(file.path(save_dir, "images", "ADT_pca.pdf"))
  print(ADT_plots)
  dev.off()
}

saveRDS(seurat_data, file.path(save_dir, "rda_obj", "seurat_processed.rds"))

