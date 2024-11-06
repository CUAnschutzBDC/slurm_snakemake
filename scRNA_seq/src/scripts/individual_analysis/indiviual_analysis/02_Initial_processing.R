library(Seurat)
library(tidyverse)
library(cowplot)
library(here)
library(scAnalysisR)
library(scuttle)
library(here)
library(openxlsx)
library(djvdj)
library(scran)

source(here("src", "scripts", "common_setup.R"))

vars.to.regress <- NULL

# Make directories
ifelse(!dir.exists(save_dir), dir.create(save_dir, recursive = TRUE), FALSE)

ifelse(!dir.exists(file.path(save_dir, "images")),
       dir.create(file.path(save_dir, "images")), FALSE)

ifelse(!dir.exists(file.path(save_dir, "files")),
       dir.create(file.path(save_dir, "files")), FALSE)

ifelse(!dir.exists(file.path(save_dir, "rda_obj")),
       dir.create(file.path(save_dir, "rda_obj")), FALSE)

mt_pattern <- "^MT-" # "^MT-" for human, "^mt-" for mice
rb_pattern <- "^RP[SL]" # ^RP[SL] for human ^Rp[s|l] for mice

# Create seurat object ---------------------------------------------------------
if(VDJ_T | VDJ_B){
  tenx_structure <- "multi7"
} else {
  tenx_structure <- "count"
}
seurat_object <- create_seurat_object(
  sample = sample,
  count_path = results_dir,
  ADT = ADT, hashtag = HTO,
  tenx_structure = tenx_structure,
  min_features = 200
)

# Add mitochondrial percent
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object,
                                                      pattern = mt_pattern)

seurat_object[["percent.ribo"]] <- PercentageFeatureSet(seurat_object,
                                                        pattern = rb_pattern)

# Add in scar counts -----------------------------------------------------------
if (ADT) {
  scar_counts <- read.csv(file.path(save_dir,
                                    "files", "scar_denoised.csv"),
                          row.names = 1) %>%
    t()
  
  scar_counts <- scar_counts[ , colnames(scar_counts) %in%
                                colnames(seurat_object)]
  
  seurat_object[["SCAR_ADT"]] <- CreateAssayObject(counts = scar_counts)
  
  seurat_object <- NormalizeData(
    seurat_object, assay = "SCAR_ADT",
    normalization.method = "CLR",
    margin = 2
  )  
}

  
# Use scuttle for cutoffs ------------------------------------------------------
se <- suppressWarnings(as.SingleCellExperiment(seurat_object))
is.mito <- grep(mt_pattern, rownames(se))
per.cell <- perCellQCMetrics(se, subsets=list(Mito=is.mito))

qc.stats <- perCellQCFilters(
  per.cell,
  sum.field = "sum",
  detected.field = "detected",
  sub.fields = "subsets_Mito_percent"
)

all_info <- cbind(per.cell, qc.stats)
colnames(all_info) <- paste0("cell_qc_", colnames(all_info))
all_info <- data.frame(all_info)
seurat_object <- AddMetaData(seurat_object, metadata = all_info)

rna_qual <- featDistPlot(seurat_object,
                          geneset = c("nFeature_RNA", "nCount_RNA",
                                      "percent.mt", "percent.ribo"),
                          plot_type = "violin", col_by = "cell_qc_discard")

if(ADT){
  adt_qual <- featDistPlot(seurat_object,
                            geneset = c("nFeature_ADT", "nCount_ADT",
                                        "percent.mt", "percent.ribo"),
                            plot_type = "violin", col_by = "cell_qc_discard")
}

pdf(file.path(save_dir, "images", "quality_plots.pdf"),
    height = 10, width = 8)
plot(rna_qual)
if(ADT){
  plot(adt_qual)
}
dev.off()

# Save before moving on
saveRDS(seurat_object, file = file.path(save_dir, "rda_obj",
                                        "seurat_unfilt.rds"))


seurat_object <- subset(seurat_object, subset = cell_qc_discard == FALSE)

# Normalization ----------------------------------------------------------------
# First move the RNA assay to be the seurat norm
seurat_object[["RNASEU"]] <- seurat_object[["RNA"]]
# Scran normalization
# Need to get to natural log 
# https://github.com/LTLA/scuttle/issues/12#issuecomment-871196592
se <- se[, colnames(se) %in% colnames(seurat_object)]
clusters <- quickCluster(se)
se <- computeSumFactors(se, clusters = clusters)
se <- logNormCounts(se)

# Pull out normalized counts
norm_counts <- assay(se, "logcounts")

# Change to ln by multiplying by log(2)
norm_counts <- norm_counts * log(2)

seurat_object[["RNA"]] <- CreateAssay5Object(
  counts = GetAssayData(object = seurat_object, 
                        assay = "RNASEU", 
                        layer = "counts"), 
  data = norm_counts
)

seurat_object <- seurat_object %>%
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = vars.to.regress)

# Single Cell Transform normalization
# Workaround because of this weird globals error
# Takes too long and too much memory, I don't use it, scran has been shown
# to be better.
# options(future.globals.maxSize = 1000 * 1024^2)
# seurat_object <- SCTransform(seurat_object, 
#                              vars.to.regress = vars.to.regress,
#                              verbose = FALSE)

# Default normalization
DefaultAssay(seurat_object) <- "RNASEU"
seurat_object <- NormalizeData(seurat_object) %>% 
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = vars.to.regress)

DefaultAssay(seurat_object) <- "RNA"
# Add in cell cycle
seurat_object <- CellCycleScoring(seurat_object, 
                                  g2m.features = cc.genes$g2m.genes,
                                  s.features = cc.genes$s.genes)


# VDJ --------------------------------------------------------------------------
if (VDJ_T) {
  vdj_dir_t <- file.path(results_dir, sample, "outs", "per_sample_outs",
                         sample, "vdj_t")
  
  seurat_object <- import_vdj(
    input = seurat_object,
    vdj_dir = vdj_dir_t,
    prefix = "tcr_",
    filter_paired = FALSE,
    include_mutations = FALSE)
}
if (VDJ_B){
  vdj_dir_b <- file.path(results_dir, sample, "outs", "per_sample_outs",
                         sample, "vdj_b")
  
  seurat_object <- import_vdj(
    input = seurat_object,
    vdj_dir = vdj_dir_b,
    prefix = "bcr_",
    filter_paired = FALSE,
    include_mutations = TRUE)
  
}
saveRDS(seurat_object, file = file.path(save_dir, "rda_obj",
                                        "seurat_start.rds"))
