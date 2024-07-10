library(Seurat)
library(tidyverse)
library(cowplot)
library(here)
library(scAnalysisR)
library(scuttle)
library(here)
library(openxlsx)
library(djvdj)

source(here("src/scripts/common_setup.R"))

# TODO add in fix here for human/mouse
mt_pattern <- "^MT-" # "^MT-" for human, "^mt-" for mice


# Make directories
ifelse(!dir.exists(save_dir), dir.create(save_dir, recursive = TRUE), FALSE)

ifelse(!dir.exists(file.path(save_dir, "images")),
       dir.create(file.path(save_dir, "images")), FALSE)

ifelse(!dir.exists(file.path(save_dir, "files")),
       dir.create(file.path(save_dir, "files")), FALSE)

ifelse(!dir.exists(file.path(save_dir, "rda_obj")),
       dir.create(file.path(save_dir, "rda_obj")), FALSE)

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
  # Add in scar counts ---------------------------------------------------------
  scar_counts <- read.csv(file.path(save_dir,
                                    "files", "scar_denoised.csv"),
                          row.names = 1) %>%
    t()

  scar_counts <- scar_counts[ , colnames(scar_counts) %in%
                              colnames(seurat_object)]

  seurat_object[["SCAR_ADT"]] <- CreateAssayObject(counts = scar_counts)

  seurat_object <- NormalizeData(seurat_object, assay = "SCAR_ADT",
                                normalization.method = "CLR",
                                margin = 2)

  adt_qual <- featDistPlot(seurat_object,
                            geneset = c("nFeature_ADT", "nCount_ADT"),
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
# TODO add in fix here for human/mouse
seurat_object <- CellCycleScoring(seurat_object, 
                                  g2m.features = cc.genes$g2m.genes,
                                  s.features = cc.genes$s.genes)


# VDJ --------------------------------------------------------------------------
add_vdj <- function(results_dir, sample, seurat_object, vdj_type){
  if (vdj_type == "VDJ_B"){
    out_dir <- "vdj_b"
    prefix <- "bcr_"
  } else if (vdj_type == "VDJ_T"){
    out_dir <- "vdj_t"
    prefix <- "tcr_"
  } else {
    stop("vdj_type should be VDJ_B or VDJ_T")
  }
  vdj_dir <- file.path(results_dir, sample, "outs", "per_sample_outs",
                       sample, out_dir)

  seurat_object <- import_vdj(input = seurat_object,
                              vdj_dir = vdj_dir,
                              prefix = prefix,
                              filter_paired = FALSE,
                              include_mutations = TRUE)

  return(seurat_object)
}
# Read in VDJ data
if (VDJ_B){
  seurat_object <- add_vdj(results_dir = results_dir,
                           sample = sample,
                           seurat_object = seurat_object,
                           vdj_type = "VDJ_B")
}

if (VDJ_T){
  seurat_object <- add_vdj(results_dir = results_dir,
                           sample = sample,
                           seurat_object = seurat_object,
                           vdj_type = "VDJ_T")
}

saveRDS(seurat_object, file = file.path(save_dir, "rda_obj",
                                        "seurat_start.rds"))
