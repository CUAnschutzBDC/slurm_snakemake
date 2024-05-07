library(Seurat)
library(tidyverse)
library(cowplot)
library(here)
library(scAnalysisR)
library(muscat)
library(scran)
library(DESeq2)
library(UpSetR)

source(here("src/scripts/muscat_plotting_functions.R"))

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
                    "DN2", "Early_memory", "BND2",
                    "Memory_IgA", "Activated_memory", "Activated_naive")

seurat_data <- subset(seurat_data, subset = RNA_celltype %in% celltypes_keep)
seurat_data <- subset(seurat_data, subset = Status %in% c("nd", "no"))

# Any diabetes antigen and no binding
seurat_data$test_id <- ifelse(seurat_data$hash.ID %in% c("INS-tet",
                                                         "GAD-tet",
                                                         "IA2-tet"),
                              "diabetes_antigen", ifelse(seurat_data$hash.ID 
                                                         == "Negative",
                                                         "Negative", "other"))

seurat_data$frozen_cell_bin <- ifelse(seurat_data$millions.of.cells.frozen < 11,
                                      "10_or_less",
                                      ifelse(is.na(seurat_data$millions.of.cells.frozen),
                                             "unknown", "10_or_more"))

Idents(seurat_data) <- "test_id"

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

# All B cells ------------------------------------------------------------------

## Rename columns for compatability --------------------------------------------

# Name columns for muscat sample_id, group_id, cluster_id
seurat_data$sample_id <- paste(seurat_data$sample, seurat_data$test_id,
                               sep = "_")
seurat_data$group_id <- seurat_data$new_status

# The cluster id I'm not sure about, I'll do once with cluster and once
# with full groups
seurat_data$celltype <- "bcells"
seurat_data$cluster_id <- seurat_data$celltype


# Change to single cell experiment
sce <- as.SingleCellExperiment(seurat_data)
sce <- prepSCE(sce,
               kid = "cluster_id",
               gid = "group_id",
               sid = "sample_id",
               drop = FALSE)

# Aggregate data to pseudotbulk using default parameters from tutorial
pb <- aggregateData(sce, assay = "counts",
                    fun = "sum",
                    by = c("cluster_id", "sample_id"))


pb_mds <- pbMDS(pb)
pb_mds

### Check for batch ------------------------------------------------------------
plot_data <- pb_mds$data
# fix the dates
colData(pb)$Collection_Date <- zoo::as.Date(colData(pb)$Collection.Date)
colData(pb)$Collection_Date <- gsub("-", "_", colData(pb)$Collection_Date)
colData(pb)$processed_date <- zoo::as.Date(as.numeric(colData(pb)$date.processed.for.scSeq))
colData(pb)$processed_date <- gsub("-", "_", colData(pb)$processed_date)


colData(pb)$Collection_Date <- ifelse(colData(pb)$Collection_Date == 
                                        "2022_02_08", "2023_02_08",
                                      colData(pb)$Collection_Date)


### Run batch correction -------------------------------------------------------


pb$date_age <- paste(pb$processed_date, pb$age_bin,
                     sep = "_")

pb$date_age <- factor(pb$date_age)
pb$test_full <- gsub(" ", "_", pb$test_full)
pb$test_full <- factor(pb$test_full)
pb$test_full_new_status <- factor(pb$test_full_new_status)

date_mapping <- c("2023_01_20" = 1,
                  "2023_02_01" = 2,
                  "2023_02_08" = 3,
                  "2023_02_10" = 4,
                  "2023_02_23" = 5,
                  "2023_03_02"= 6,
                  "2023_03_09" = 7,
                  "2023_03_27" = 8,
                  "2023_03_28" = 9,
                  "2023_04_04" = 10,
                  "2023_04_10" = 11, 
                  "2023_04_11" = 12)

pb$new_date <- date_mapping[pb$Collection_Date]

date_bin <- c("2023_01_20" = 1,
              "2023_02_01" = 2,
              "2023_02_08" = 2,
              "2023_02_10" = 2,
              "2023_02_23" = 2,
              "2023_03_02"= 3,
              "2023_03_09" = 3,
              "2023_03_27" = 3,
              "2023_03_28" = 3,
              "2023_04_04" = 4,
              "2023_04_10" = 4, 
              "2023_04_11" = 4)

pb$date_bin <- date_bin[pb$Collection_Date]

pb$date_bin <- factor(pb$date_bin)

pb$full_date <- paste(pb$processed_date, pb$date_bin,
                      sep = "_")

pb$full_date <- factor(pb$full_date)
# Definitely grouping by individual, somewhat by family, can't correct by
dds <- DESeqDataSet(pb, design = ~Collection_Date + test_full_new_status)

vsd_one <- vst(dds, blind = F)

vsd_one$frozen_cell_bin <- factor(vsd_one$frozen_cell_bin)

vsd_one$age_bin <- factor(vsd_one$age_bin)

mat <- assay(vsd_one)

# Day is the batch, the model.matrix should include everything but the batch
mat <- limma::removeBatchEffect(mat, vsd_one$Collection.Date,
                                design = model.matrix(~vsd_one$test_full))

vsd_two <- vsd_one
assay(vsd_two) <- mat

mod_mat <- model.matrix(design(dds), colData(dds))


ifelse(!dir.exists(file.path(save_dir, "images", "pseudobulk")),
       dir.create(file.path(save_dir, "images", "pseudobulk")),
       FALSE)
pdf(file.path(save_dir, "images", "pseudobulk", "all_b_cells_no_aab.pdf"))


print(plot_pca(vsd_one, group_by = "frozen_cell_bin") +
        ggplot2::ggtitle("uncorrected"))

print(plot_pca(vsd_one, group_by = "age_bin") +
        ggplot2::ggtitle("uncorrected"))

vsd_one$num_Pe_neg_sorted_cells_loaded_on_10X <- factor(vsd_one$num_Pe_neg_sorted_cells_loaded_on_10X)
vsd_one$num_PE_pos_sorted_cells_loaded_on_10X <- factor(vsd_one$num_PE_pos_sorted_cells_loaded_on_10X)


color_one <- viridis::magma(n = 256)

color_scale <- grDevices::colorRampPalette(color_one)

neg_scale <- color_scale(length(levels(vsd_one$num_Pe_neg_sorted_cells_loaded_on_10X)))
names(neg_scale) <- levels(vsd_one$num_Pe_neg_sorted_cells_loaded_on_10X)

vsd_one$Collection_Date <- factor(vsd_one$Collection_Date)
collection_scale <- color_scale(length(levels(vsd_one$Collection_Date)))
names(collection_scale) <- levels(vsd_one$Collection_Date)


print(plot_pca(vsd_one, group_by = "num_Pe_neg_sorted_cells_loaded_on_10X",
               color_palette = neg_scale) +
        ggplot2::ggtitle("uncorrected"))

pos_scale <- color_scale(length(levels(vsd_one$num_PE_pos_sorted_cells_loaded_on_10X)))
names(pos_scale) <- levels(vsd_one$num_PE_pos_sorted_cells_loaded_on_10X)

print(plot_pca(vsd_one, group_by = "num_PE_pos_sorted_cells_loaded_on_10X",
               color_palette = pos_scale) +
        ggplot2::ggtitle("uncorrected"))

print(plot_pca(vsd_one, group_by = "Collection_Date",
               color_palette = collection_scale) +
        ggplot2::ggtitle("uncorected"))

vsd_one$millions.of.cells.frozen <- factor(vsd_one$millions.of.cells.frozen)

print(plot_pca(vsd_one, group_by = "millions.of.cells.frozen") + 
        ggplot2::ggtitle("uncorrected"))

print(plot_pca(vsd_one, group_by = "new_status", color_palette = status_colors) +
        ggplot2::ggtitle("uncorrected"))

print(plot_pca(vsd_one, group_by = "sample", color_palette = sample_colors) +
        ggplot2::ggtitle("uncorrected"))

print(plot_pca(vsd_one, group_by = "test_id", color_palette = antigen_colors) +
        ggplot2::ggtitle("uncorrected"))
print(plot_pca(vsd_one, group_by = "processed_date") +
        ggplot2::ggtitle("uncorrected"))

print(plot_pca(vsd_two, group_by = "new_status", color_palette = status_colors) +
        ggplot2::ggtitle("corrected"))

print(plot_pca(vsd_two, group_by = "Status", color_palette = status_colors2) +
        ggplot2::ggtitle("corrected"))

print(plot_pca(vsd_two, group_by = "sample", color_palette = sample_colors) +
        ggplot2::ggtitle("corrected"))

print(plot_pca(vsd_two, group_by = "test_id", color_palette = antigen_colors) +
        ggplot2::ggtitle("corrected"))

dev.off()

### Run deseq with muscat and the new batch term -------------------------------

formula <- ~full_date + test_full
cd <- as.data.frame(colData(pb))
design <- model.matrix(formula, cd)

ND <- colMeans(design[which(cd$group_id == "nd"), ])
NO <- colMeans(design[which(cd$group_id == "no"), ])
AAB <- colMeans(design[which(cd$group_id == "aab"), ])
diabetes_antigen <- colMeans(design[which(cd$test_id == "diabetes_antigen"), ])
Negative <- colMeans(design[which(cd$test_id == "Negative"), ])
tet <- colMeans(design[which(cd$test_id == "Tet_antigen"), ])
other <- colMeans(design[which(cd$test_id == "other"), ])
ND_diabetes_antigen <- colMeans(design[which(cd$test_id == "diabetes_antigen"
                                             & cd$group_id == "nd") , ])
ND_Negative <- colMeans(design[which(cd$test_id == "Negative" 
                                     & cd$group_id == "nd"), ])
ND_tet <- colMeans(design[which(cd$test_id == "Tet_antigen"
                                & cd$group_id == "nd"), ])
ND_other <- colMeans(design[which(cd$test_id == "other"
                                  & cd$group_id == "nd"), ])

NO_diabetes_antigen <- colMeans(design[which(cd$test_id == "diabetes_antigen"
                                             & cd$group_id == "no") , ])
NO_Negative <- colMeans(design[which(cd$test_id == "Negative" 
                                     & cd$group_id == "no"), ])
NO_tet <- colMeans(design[which(cd$test_id == "Tet_antigen"
                                & cd$group_id == "no"), ])
NO_other <- colMeans(design[which(cd$test_id == "other"
                                  & cd$group_id == "no"), ])

AAB_diabetes_antigen <- colMeans(design[which(cd$test_id == "diabetes_antigen"
                                              & cd$group_id == "aab") , ])
AAB_Negative <- colMeans(design[which(cd$test_id == "Negative" 
                                      & cd$group_id == "aab"), ])
AAB_tet <- colMeans(design[which(cd$test_id == "Tet_antigen"
                                 & cd$group_id == "aab"), ])
AAB_other <- colMeans(design[which(cd$test_id == "other"
                                   & cd$group_id == "aab"), ])

all_contrasts <- limma::makeContrasts(NO - ND,
                                      diabetes_antigen - Negative,
                                      diabetes_antigen - tet,
                                      NO_diabetes_antigen - ND_diabetes_antigen,
                                      NO_Negative - ND_Negative,
                                      NO_tet - ND_tet, 
                                      NO_other - ND_other,
                                      levels = design)

de_genes_all <- muscat::pbDS(pb, method = "DESeq2",
                             design = design, contrast = all_contrasts)

all_res <- lapply(de_genes_all$table, function(x){
  return(x$bcells)
})

all_res <- do.call(rbind, all_res)

all_res <- all_res %>%
  dplyr::filter(p_adj.glb < 0.05)

# Find percent of cells expressing genes
freqs <- calcExprFreqs(sce, assay = "counts", th = 0)

all_freqs <- assay(freqs)

# Test against a cutoff
gids <- levels(sce$group_id)
frq10 <- vapply(as.list(assays(freqs)), 
                function(u) apply(u[, gids] > 0.1, 1, any), 
                logical(nrow(sce)))

frq10 <- data.frame(frq10)

frq10 <- frq10[frq10$bcells,, drop = FALSE]

all_res <- all_res[all_res$gene %in% rownames(frq10),]

# Save the results
write.csv(all_res, file.path(save_dir, "files",
                             "pseudobulk_all_bcell_de_no_aab.csv"))

de_wb <- openxlsx::createWorkbook()

for(i in unique(all_res$contrast)){
  subset_res <- all_res %>%
    dplyr::filter(contrast == i)
  
  if(grepl("diabetes_antigen", i)){
    i <- gsub("diabetes_antigen", "da", i)
  }
  
  openxlsx::addWorksheet(wb = de_wb, sheetName = i)
  openxlsx::writeData(wb = de_wb, sheet = i, x = subset_res)
}

openxlsx::saveWorkbook(wb = de_wb,
                       file = file.path(save_dir, "files",
                                        "pseudobulk_all_bcell_de_no_aab.xlsx"),
                       overwrite = TRUE)

saveRDS(pb, file.path(save_dir, "rda_obj", "muscat_object_no_aab.rds"))
saveRDS(vsd_two, file.path(save_dir, "rda_obj", 
                           "psuedobulk_batch_corrected_object_no_aab.rds"))


