library(Seurat)
library(tidyverse)
library(cowplot)
library(here)
library(scAnalysisR)
library(scuttle)
library(here)
library(openxlsx)
library(djvdj)


# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

args <- commandArgs(trailingOnly = TRUE)

sample <- args[[1]]
sample <- gsub("__.*", "", sample)
#sample <- "Npod6553_PLN"

sample_info <- args[[3]]
#sample_info <- here("files/sample_info.tsv")

results_dir <- args[[2]]
#results_dir <- here("results")

sample_info <- read.table(sample_info, fill = TRUE, header = TRUE)

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

vars.to.regress <- NULL

# Set directories
save_dir <- file.path(results_dir, "R_analysis", sample)

# Make directories
ifelse(!dir.exists(save_dir), dir.create(save_dir, recursive = TRUE), FALSE)

ifelse(!dir.exists(file.path(save_dir, "images")),
       dir.create(file.path(save_dir, "images")), FALSE)

ifelse(!dir.exists(file.path(save_dir, "files")),
       dir.create(file.path(save_dir, "files")), FALSE)

ifelse(!dir.exists(file.path(save_dir, "rda_obj")),
       dir.create(file.path(save_dir, "rda_obj")), FALSE)

mt_pattern <- "^MT-" # "^MT-" for human, "^mt-" for mice

# Create seurat object ---------------------------------------------------------
seurat_object <- create_seurat_object(sample = sample,
                                      count_path = results_dir,
                                      ADT = ADT, hashtag = HTO,
                                      tenx_structure = "multi7",
                                      min_features = 200
                                      )

# Add mitochondrial percent
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object,
                                                      pattern = mt_pattern)

# Add in scar counts -----------------------------------------------------------
scar_counts <- read.csv(file.path(save_dir,
                                  "files", "scar_denoised.csv"),
                        row.names = 1) %>%
  t()

scar_counts <- scar_counts[ , colnames(scar_counts) %in%
                             colnames(seurat_object)]

seurat_object[["SCAR_ADT"]] <- CreateAssayObject(counts = scar_counts)

seurat_object <- NormalizeData(seurat_object, assay = "SCAR_ADT",
                               normalization.method = "LogNormalize")
  
# Use scuttle for cutoffs ------------------------------------------------------
se <- as.SingleCellExperiment(seurat_object)
is.mito <- grep(mt_pattern, rownames(se))
per.cell <- perCellQCMetrics(se, subsets=list(Mito=is.mito))

qc.stats <- perCellQCFilters(per.cell,
                             sum.field = "sum",
                             detected.field = "detected",
                             sub.fields = "subsets_Mito_percent")

all_info <- cbind(per.cell, qc.stats)
colnames(all_info) <- paste0("cell_qc_", colnames(all_info))
all_info <- data.frame(all_info)
seurat_object <- AddMetaData(seurat_object, metadata = all_info)

rna_qual <- featDistPlot(seurat_object,
                          geneset = c("nFeature_RNA", "nCount_RNA",
                                      "percent.mt"),
                          plot_type = "violin", col_by = "cell_qc_discard")

if(ADT){
  adt_qual <- featDistPlot(seurat_object,
                            geneset = c("nFeature_ADT", "nCount_ADT",
                                        "percent.mt"),
                            plot_type = "violin", col_by = "cell_qc_discard")
}

pdf(file.path(save_dir, "images", "quality_plots.pdf"))
plot(rna_qual)
if(ADT){
  plot(adt_qual)
}
dev.off()

# Save before moving on
saveRDS(seurat_object, file = file.path(save_dir, "rda_obj",
                                        "seurat_unfilt.rds"))


seurat_object <- subset(seurat_object, subset = cell_qc_discard == FALSE)

# Remove ambeince --------------------------------------------------------------

# Normalization ----------------------------------------------------------------
# Single Cell Transform normalization
seurat_object <- SCTransform(seurat_object, vars.to.regress = vars.to.regress,
                             verbose = FALSE)

# Default normalization
DefaultAssay(seurat_object) <- "RNA"
seurat_object <- NormalizeData(seurat_object) %>% 
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = vars.to.regress)

# Add in cell cycle
seurat_object <- CellCycleScoring(seurat_object, 
                                  g2m.features = cc.genes$g2m.genes,
                                  s.features = cc.genes$s.genes)


# VDJ --------------------------------------------------------------------------
# Read in VDJ data
vdj_dir_b <- file.path(results_dir, sample, "outs/per_sample_outs",
                       sample, "vdj_b")

vdj_dir_t <- file.path(results_dir, sample, "outs/per_sample_outs",
                       sample, "vdj_t")

seurat_object <- import_vdj(input = seurat_object,
                            vdj_dir = vdj_dir_b,
                            prefix = "bcr_",
                            filter_paired = FALSE,
                            include_mutations = TRUE)


seurat_object <- import_vdj(input = seurat_object,
                            vdj_dir = vdj_dir_t,
                            prefix = "tcr_",
                            filter_paired = FALSE,
                            include_mutations = TRUE)

saveRDS(seurat_object, file = file.path(save_dir, "rda_obj",
                                        "seurat_start.rds"))
