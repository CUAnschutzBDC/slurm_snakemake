# Following this tutorial
# https://bioinformatics-core-shared-training.github.io/UnivCambridge_ScRnaSeq_Nov2021/Markdowns/10_MultiSplComp.html
library(Seurat)
library(tidyverse)
library(cowplot)
library(here)
library(scAnalysisR)
library(scran)
library(scuttle)
library(UpSetR)
library(edgeR)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

#args <- commandArgs(trailingOnly = TRUE)

#sample <- args[[1]]
#sample <- gsub("__.*", "", sample)
sample <- "merged"

#sample_info <- args[[4]]
sample_info <- here("files/sample_info.tsv")

#results_dir <- args[[2]]
results_dir <- here("results")

#sample_metadata <- args[[3]]
sample_metadata <- here("files/Deidentified_donor_metadata.xlsx")
sample_metadata <- openxlsx::readWorkbook(sample_metadata, detectDates = TRUE,
                                          sheet = "New_metadata") %>%
  dplyr::select(Sample.Name, millions.of.cells.frozen, dplyr::contains("num_"))

sample_metadata$Sample.Name <- gsub("JH313-15_JBM",
                                    "JH313-15_JB",
                                    sample_metadata$Sample.Name)

sample_info <- read.table(sample_info, fill = TRUE, header = TRUE)

samples_use <- sample_info[sample_info$sample != sample, ]$sample

sample_info <- sample_info[sample_info$sample == sample,]

HTO <- sample_info$HTO
ADT <- sample_info$ADT
hash_ident <- sample_info$hash_ident


ADT_pca <- sample_info$adt_pca
run_adt_umap <- sample_info$adt_umap

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
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed_no_doublet.rds"))

seurat_data$age_bin <- ifelse(seurat_data$Age.at.Collection..years. <= 10,
                              "less_than_10",
                              ifelse(seurat_data$Age.at.Collection..years. <= 20,
                                     "10_to_20", "more_than_20"))

sample_metadata <- sample_metadata %>%
  dplyr::filter(Sample.Name %in% seurat_data$Sample.Name)


full_seurat_meta <- seurat_data[[]] %>%
  dplyr::select(Sample.Name) %>%
  tibble::rownames_to_column("barcode")

sample_metadata <- merge(sample_metadata, full_seurat_meta, all.x = FALSE,
                         all.y = TRUE)

sample_metadata <- sample_metadata %>%
  tibble::column_to_rownames("barcode") %>%
  dplyr::select(!Sample.Name)

seurat_data <- AddMetaData(seurat_data, metadata = sample_metadata)

# Colors -----------------------------------------------------------------------
all_samples <- unique(seurat_data$sample)
sample_colors <- MetBrewer::met.brewer(name = "Archambault", 
                                       n = length(all_samples),
                                       type = "continuous")

names(sample_colors) <- all_samples

antigen_colors <- MetBrewer::met.brewer(name = "Demuth", n = 4,
                                        type = "discrete")

names(antigen_colors) <- c("Negative", "Tet_antigen", "other",
                           "diabetes_antigen")

status_colors <- MetBrewer::met.brewer(name = "Egypt", n = 3,
                                       type = "discrete")

names(status_colors) <- c("nd", "no", "aab")


status_colors2 <- MetBrewer::met.brewer(name = "Egypt", n = 4,
                                        type = "discrete")

names(status_colors2) <- c("nd", "no", "aab stage 1", "aab stage 2")

celltype_colors <- readRDS(file.path(save_dir, "color_palette.rds"))

# Update seurat object ---------------------------------------------------------

# First remove all non-b cells
celltypes_keep <- c("Naive_1", "Naive_2", "Naive_3",
                    "BND_cluster", "B.intermediate",
                    "Early_memory", "Memory_IgE_IgG", "Resting_memory",
                    "Plasmablast", "DN2", "Early_memory", "BND2",
                    "Memory_IgA", "Activated_memory", "Activated_naive")

seurat_data <- subset(seurat_data, subset = RNA_celltype %in% celltypes_keep)

name_mapping <- c("INS-tet" = "diabetes_antigen",
                  "GAD-tet" = "diabetes_antigen",
                  "IA2-tet" = "diabetes_antigen",
                  "TET-tet" = "Tet_antigen",
                  "Negative" = "Negative",
                  "Doublet" = "other",
                  "DNA-tet" = "other")

seurat_data$test_id <- name_mapping[as.character(seurat_data$hash.ID)]

Idents(seurat_data) <- "test_id"

seurat_data$test_full <- paste(seurat_data$Status, seurat_data$test_id,
                               sep = "_")

seurat_data$new_status <- gsub("aab stage [1|2]", "aab", seurat_data$Status)

seurat_data$test_full_new_status <- paste(seurat_data$new_status, seurat_data$test_id,
                                          sep = "_")

sce <- as.SingleCellExperiment(seurat_data)

summed <- aggregateAcrossCells(sce,
                               id = DataFrame(label = sce$test_id,
                                              sample = sce$sample))

# Subset to just one -----------------------------------------------------------
# Think about this, should I do it on everything?
labelToGet <- "diabetes_antigen"
current <- summed[,summed$label==labelToGet]

countsToUse <- counts(current)
colnames(countsToUse) <- colData(current)$sample

# Preprocess -------------------------------------------------------------------
y <- DGEList(countsToUse, samples = colData(current))

discarded <- current$ncells < 20
y <- y[,!discarded]
summary(discarded)

keep <- filterByExpr(y, group=current$Status)
y <- y[keep,]
summary(keep)

y <- calcNormFactors(y)
y$samples

par(mfrow=c(2,4))
for (i in seq_len(ncol(y))) {
  plotMD(y, column=i)
}

mapping <- c("nd" = "blue",
             "aab stage 1" = "darkgreen",
             "aab stage 2" = "red",
             "no" = "black")

y$samples$color <- mapping[y$samples$Status]

par(mfrow = c(1,1))
plotMDS(cpm(y, log=TRUE), 
        col = y$samples$color
)

# Modeling ---------------------------------------------------------------------
design <- model.matrix(~date.processed.for.scSeq + factor(Status), y$samples)
design

y <- estimateDisp(y, design)
summary(y$trended.dispersion)

plotBCV(y)

fit <- glmQLFit(y, design, robust=TRUE)
summary(fit$var.prior)

summary(fit$df.prior)

plotQLDisp(fit)

res <- glmQLFTest(fit, coef=ncol(design))
summary(decideTests(res))
