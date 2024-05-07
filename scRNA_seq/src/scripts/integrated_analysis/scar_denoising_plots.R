library(here)
library(scAnalysisR)
library(tidyverse)
library(Seurat)
source(here("src/scripts/muscat_plotting_functions.R"))


# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

# Set directories
base_dir <- here()
sample <- "merged"


save_dir <- file.path(base_dir, "results", "R_analysis", sample)

ifelse(!dir.exists(file.path(save_dir, "images", "scar_denoising")),
       dir.create(file.path(save_dir, "images", "scar_denoising")),
       FALSE)

# Read in the data
seurat_data <- readRDS(file.path(save_dir, "rda_obj",
                                 "seurat_processed_no_doublet.rds"))

all_samples <- unique(seurat_data$sample)

all_plots <- lapply(all_samples, function(test_sample){
  
  scar_counts <- read.csv(file.path(base_dir, "results", "R_analysis", test_sample,
                                    "files", "scar_denoised.csv"),
                          row.names = 1) 
  
  seurat_data_sub <- subset(seurat_data, subset = sample == test_sample)
  
  raw_counts <- GetAssayData(seurat_data_sub, assay = "ADT", slot = "counts")
  
  colnames(raw_counts) <- gsub("_[0-9]+", "", colnames(raw_counts))
  
  raw_counts <- t(as.matrix(raw_counts))
  
  # First remove the incorrect insulin and other unnecessary ADTs
  all_sums <- colSums(scar_counts)
  
  # 1000 is arbitrary but should work to identify the actual ADTs. This should
  # be checked and adjusted.
  keep_adts <- all_sums[all_sums > 1000] 
  
  scar_counts <- scar_counts[,names(keep_adts)]
  
  colnames(scar_counts) <- gsub("INS_tet[a|b]", "INS_tet", colnames(scar_counts))
  
  scar_counts <- scar_counts[rownames(scar_counts) %in% rownames(raw_counts),]
  
  colnames(scar_counts) = gsub("_", ".", colnames(scar_counts))
  
  scar_counts <- scar_counts %>%
    tidyr::pivot_longer(cols = colnames(scar_counts), names_to = "tetramer",
                        values_to = "count") %>%
    dplyr::mutate(type = "scar")
  
  raw_counts <- raw_counts %>%
    data.frame %>%
    tidyr::pivot_longer(cols = colnames(.), names_to = "tetramer",
                        values_to = "count") %>%
    dplyr::mutate(type = "raw")
  
  all_counts <- rbind(scar_counts, raw_counts) %>%
    dplyr::mutate(count = log1p(count))
  
  scar_plot <- ggplot2::ggplot(all_counts, ggplot2::aes(x = tetramer,
                                                        y = count,
                                                        fill = type)) + 
    ggplot2::geom_violin() +
    ggplot2::geom_violin(scale = "width",
                         position = ggplot2::position_dodge(1),
                         width = 0.9) +
    ggplot2::stat_summary(fun = median, geom = "point", size = 2,
                          position = ggplot2::position_dodge(1)) +
    ggplot2::scale_fill_brewer(palette = "Set1") +
    ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  
  pdf(file.path(save_dir, "images", "scar_denoising",
                paste0(test_sample, ".pdf")))
  
  print(scar_plot)
  
  dev.off()  
})
