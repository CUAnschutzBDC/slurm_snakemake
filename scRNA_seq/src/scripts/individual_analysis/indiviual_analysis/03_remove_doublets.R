library(DoubletFinder)
library(Seurat)
library(tidyverse)
library(here)
library(scAnalysisR)

source(here("src", "scripts", "common_setup.R"))

# Read in data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_start.rds"))

# Remove "negatives"
if(HTO){
  Idents(seurat_data) <- "HTO_classification.global"
  seurat_data <- subset(x = seurat_data, idents = c("Singlet", "Doublet"))
}

# PCA and UMAP
set.seed(0)
DefaultAssay(seurat_data) <- seurat_assay
seurat_data <- PCA_dimRed(seurat_data, assay = seurat_assay)
npcs <- 20

umap_data <- group_cells(seurat_data, assay = seurat_assay, nPCs = npcs)

seurat_data <- umap_data[[1]]

# Run Doublet finder -----------------------------------------------------------
set.seed(0)
sweep.res.seurat_data <- paramSweep(seurat_data, PCs = 1:npcs, sct = SCT)

if(HTO){
  ## pK Identification (ground-truth)
  gt.calls <- seurat_data@meta.data[rownames(sweep.res.seurat_data[[1]]),
                                    "HTO_classification.global"]
  gt.calls <- as.factor(gt.calls)
  sweep.stats_sample <- summarizeSweep(sweep.res.seurat_data, GT = TRUE,
                                       GT.calls = gt.calls)
  pdf(file.path(save_dir, "images", "doublet_finder_peak.pdf"))
  
  bcmvn_sample <- find.pK(sweep.stats_sample)
  dev.off()
  max_AUC <- max(bcmvn_sample$MeanAUC)
  max_BC <- max(bcmvn_sample$MeanBC)
  max_BCmetric <- max(bcmvn_sample$BCmetric)
} else {
  ## pK Identification (no ground-truth)
  sweep.stats_sample <- summarizeSweep(sweep.res.seurat_data, GT = FALSE)
  
  pdf(file.path(save_dir, "images", "doublet_finder_peak.pdf"))
  
  bcmvn_sample <- find.pK(sweep.stats_sample)
  dev.off()
  max_BC <- max(bcmvn_sample$MeanBC)
  max_BCmetric <- max(bcmvn_sample$BCmetric)
}


# Pull out best pK
best_pk <- bcmvn_sample[bcmvn_sample$BCmetric == max_BCmetric,]$pK

pK <- as.double(as.character(best_pk))

## Homotypic Doublet Proportion Estimate ---------------------------------------
annotations <- seurat_data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) 

# Find expected proportion -----------------------------------------------------
## CHANGE BASED ON CELL NUMBERS
# Multiplet rate	# cells loaded	# cells recovered
# 0.40%	825	500
# 0.80%	1,650	1,000
# 1.60%	3,300	2,000
# 2.40%	4,950	3,000
# 3.20%	6,600	4,000
# 4.00%	8,250	5,000
# 4.80%	9,900	6,000
# 5.60%	11,550	7,000
# 6.40%	13,200	8,000
# 7.20%	14,850	9,000
# 8.00%	16,500	10,000

multiplet_mapping <- c("500" = 0.004,
                       "1000" = 0.008,
                       "2000" = 0.016,
                       "3000" = 0.024,
                       "4000" = 0.032,
                       "5000" = 0.04,
                       "6000" = 0.048,
                       "7000" = 0.056,
                       "8000" = 0.064,
                       "9000" = 0.072,
                       "10000" = 0.08)

# Find number of cells
ncells <- ncol(seurat_data)

differences <- abs(as.numeric(names(multiplet_mapping)) - ncells)

# Find the index of the element with the minimum difference
closest_index <- which.min(differences)

# Retrieve the closest number from the vector
expected_proportion <- as.numeric(multiplet_mapping[closest_index])

nExp_poi <- round(0.032*nrow(seurat_data@meta.data))  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


## Find pK from bcmvn output, pN selection is less important--------------------

seurat_data <- doubletFinder(
  seurat_data, PCs = 1:npcs, 
  pN = 0.25,
  pK = pK, nExp = nExp_poi.adj,
  reuse.pANN = FALSE,
  sct = SCT
)

# Rename to be a consistent name between all samples
df_column <- colnames(seurat_data[[]])[grep("DF.classifications",
                                            colnames(seurat_data[[]]))]


seurat_data$Doublet_finder <- seurat_data[[df_column]]

# Remove the original name
seurat_data[[df_column]] <- NULL

pdf(file.path(save_dir, "images", "doublet_finder.pdf"))

print(plotDimRed(seurat_data, col_by = "Doublet_finder",
                 plot_type = "rna.umap"))
print(plotDimRed(seurat_data, col_by = "JCHAIN", plot_type = "rna.umap"))

print(plotDimRed(seurat_data, col_by = "CD3D", plot_type = "rna.umap"))

dev.off()

saveRDS(seurat_data, file.path(save_dir, "rda_obj", "seurat_doublet.rds"))