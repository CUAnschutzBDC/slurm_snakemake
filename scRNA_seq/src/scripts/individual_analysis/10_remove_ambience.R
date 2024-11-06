library(Seurat)
library(tidyverse)
library(cowplot)
library(here)
library(scAnalysisR)
library(DropletUtils)
source(here("src", "scripts", "common_setup.R"))

# Read in data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))

DefaultAssay(seurat_data) <- "RNA"
# Background contamination -----------------------------------------------------
if(VDJ_T | VDJ_B){
  raw_path <- file.path(results_dir, sample, "outs", "multi",
                        "count", "raw_feature_bc_matrix")  
} else {
  raw_path <- file.path(results_dir, sample, "outs",
                        "raw_feature_bc_matrix")  
}


raw_sce <- DropletUtils::read10xCounts(raw_path, col.names = TRUE)

# Remove the antibodies
raw_sce <- raw_sce[grepl("ENS", rownames(raw_sce)), ]

# Rename genes based on gene name not ens id
gene_mapping <- read.table(file.path(raw_path, "features.tsv.gz"),
                           sep = "\t")

colnames(gene_mapping) <- c("ens_id", "gene_id", "type")

# Seurat deals with gene name duplicates with the "make names" argument
gene_mapping <- gene_mapping %>%
  dplyr::filter(type == "Gene Expression") %>%
  dplyr::mutate(gene_id = make.unique(gene_id)) %>%
  dplyr::filter(gene_id %in% rownames(seurat_data)) 

# Subset to the cells and genes kept by seurat
raw_sce <- raw_sce[rownames(raw_sce) %in% gene_mapping$ens_id, ]
filtered_sce <- raw_sce[ , colnames(raw_sce) %in% colnames(seurat_data)]

if(!identical(colnames(filtered_sce), colnames(seurat_data))){
  filtered_sce <- filtered_sce[ , order(match(colnames(filtered_sce),
                                              colnames(seurat_data)))]
}

colLabels(filtered_sce) <- seurat_data[[clusters]][[1]]

set.seed(100)
e.out <- emptyDrops(counts(raw_sce))

amb <- metadata(e.out)$ambient[,1]
stripped <- filtered_sce[names(amb),]

out <- removeAmbience(counts(stripped), ambient=amb, groups=colLabels(stripped))

# Make sure something happened
differences <- out != counts(stripped)
sums <- rowSums(differences)
changed <- sums[sums > 0]

counts(stripped, withDimnames=FALSE) <- out

final_counts <- counts(stripped)

gene_mapping <- gene_mapping %>%
  dplyr::filter(ens_id %in% rownames(final_counts))

changed_ids <- gene_mapping %>%
  dplyr::filter(ens_id %in% names(changed))

changed_ids$count <- changed[changed_ids$ens_id]

# Now order this to be the same as the counts
final_counts <- final_counts[rownames(final_counts) %in% gene_mapping$ens_id, ]

final_counts <- final_counts[order(match(rownames(final_counts),
                                         gene_mapping$ens_id)), ]

rownames(final_counts) <- gene_mapping$gene_id

seurat_data[["AMBRNA"]] <- CreateAssayObject(counts = final_counts)

seurat_data <- NormalizeData(seurat_data, assay = "AMBRNA")

saveRDS(seurat_data, file.path(save_dir, "rda_obj", "seurat_processed.rds"))
