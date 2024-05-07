library(scAnalysisR)
library(tidyverse)
library(Seurat)
library(here)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

sample <- "merged"

sample_info <- here("files/sample_info.tsv")

results_dir <- here("results")

sample_info <- read.table(sample_info, fill = TRUE, header = TRUE)

samples_use <- sample_info[sample_info$sample != sample, ]$sample

sample_info <- sample_info[sample_info$sample == sample,]

HTO <- sample_info$HTO
ADT <- sample_info$ADT
hash_ident <- sample_info$hash_ident

if(normalization_method == "SCT"){
  SCT <- TRUE
  seurat_assay <- "SCT"
} else {
  SCT <- FALSE
  seurat_assay <- "RNA"
}

save_dir <- file.path(results_dir, "R_analysis", sample)

seurat_data <- readRDS(file = file.path(save_dir, "rda_obj",
                                        "seurat_adtnorm.rds"))


# MBC CD66 = CEACAM1
gene_lists <- list("Naive" = c("BACH2", "ZBTB16", "APBB2", "SPRY1", "TCL1A",
                               "IKZF2"),
                   "MBC" = c("CD27", "CEACAM1", "RASSF6", "TOX", "TRERF1",
                             "TRPV3", "POU2AF1", "RORA", "TNFRSF13B",
                             "CD60", "FCRL5"),
                   "Class_switched" = c("GDPD5", "BAIAP3", "TGM2", "MUC16"),
                   "ASC" = c("PRDM1", "MANF", "XBP1", "IL6R", "BCL6", "IRF4",
                             "TNFRSF17"),
                   "GC_emigrant" = c("NT5E", "MKI67", "CD40", "CD83", "MAP3KB",
                                     "MAP3K1", "FAS"),
                   "Marginal_zone" = c("SOX7", "SLA", "TMEM71", "HS3ST3B1",
                                       "RPS6KA2", "RHOBTB3", "AEBP1", "DMD"),
                   "B1" = c("SPN", "MYO1D", "PLSCR1", "PSTPIP2", "AHR", 
                            "CD300LF", "LYSMD2", "FAM160A1", "CYP11A1",
                            "ZBTB32", "BHLHE41", "TBC1D9", "IZUMO1R",
                            "GPR55", "TNFSF13", "MYD88"))

all_meta <- lapply(names(gene_lists), function(x){
  new_seurat <- AddModuleScore(seurat_data, features = gene_lists[x])
  new_meta <- new_seurat[[]] %>%
    dplyr::select("Cluster1")
  
  colnames(new_meta) <- x
  
  return(new_meta)
})

all_meta <- do.call(cbind, all_meta)

seurat_data <- AddMetaData(seurat_data, all_meta)

first_violin <- featDistPlot(seurat_data, geneset = names(all_meta),
                             sep_by = "RNA_celltype", combine = FALSE)

names(first_violin) <- names(all_meta)


pdf(file.path(save_dir, "images", "covid_study_modules.pdf"))

print(first_violin)

dev.off()
