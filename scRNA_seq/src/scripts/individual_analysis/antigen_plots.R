library(Seurat)
library(tidyverse)
library(cowplot)
library(here)
library(scAnalysisR)
source(here("src", "scripts", "common_setup.R"))

# Read in data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))

table(seurat_data$RNA_celltype)

table(seurat_data$combined_celltype, seurat_data$HTO_classification)

table(seurat_data$combined_celltype, seurat_data$manual_hash_id)
cm <- confusionMatrix(seurat_data$RNA_celltype, seurat_data$manual_hash_id)
cm <- cm /rowSums(cm)

cm <- data.frame(cm)
original_colnames <- colnames(cm)
cm$cell_type <- rownames(cm)
cm <- cm[,c("cell_type", original_colnames)]

write.table(cm , file.path(save_dir, "files", "oligo_determination.csv"),
            sep = ",", row.names = FALSE)


cm <- confusionMatrix(seurat_data$combined_celltype, seurat_data$manual_hash_id)

cm <- data.frame(cm)
original_colnames <- colnames(cm)
cm$cell_type <- rownames(cm)
cm <- cm[,c("cell_type", original_colnames)]

write.table(cm , file.path(save_dir, "files", "oligo_determination_counts.csv"),
            sep = ",", row.names = FALSE)

tet <- GetAssayData(seurat_data, assay = "SCAR_ADT")
tet <- tet[!grepl("^anti", rownames(tet)),]

pdf(file.path(save_dir, "images", "antigen_by_celltype.pdf"))
print(featDistPlot(seurat_data, geneset = rownames(tet), assay = "SCAR_ADT",
             combine = FALSE, sep_by = "combined_celltype"))

dev.off()


pdf(file.path(save_dir, "images", "antigen_by_antigen_call.pdf"))
print(featDistPlot(seurat_data, geneset = rownames(tet), assay = "SCAR_ADT",
                   combine = FALSE, sep_by = "manual_hash_id"))

dev.off()