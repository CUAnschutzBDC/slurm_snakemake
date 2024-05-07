library(Seurat)
library(tidyverse)
library(cowplot)
library(openxlsx)
library(here)
library(scAnalysisR)
library(clustifyr)
library(splitstackshape)
library(data.table)


# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log


cell_types <- "combined_celltype"
clusters <- "combined_cluster"
pval <- 0.05
logfc <- 0.5

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

normalization_method <- "log" # can be SCT or log

args <- commandArgs(trailingOnly = TRUE)

sample <- args[[1]]
sample <- gsub("__.*", "", sample)
#sample <- "Npod6456_PLN"
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

ADT_pca <- sample_info$adt_pca
run_adt_umap <- sample_info$adt_umap

RNA_pcs <- sample_info$PCs
resolution <- sample_info$resolution

if(normalization_method == "SCT"){
  SCT <- TRUE
  seurat_assay <- "SCT"
} else {
  SCT <- FALSE
  seurat_assay <- "RNA"
}

# Set directories
save_dir <- file.path(results_dir, "R_analysis")

# Read in data
cd4_save_dir <- file.path(save_dir, paste0(sample, "_cd4"))
cd8_save_dir <- file.path(save_dir, paste0(sample, "_cd8"))
seurat_cd4 <- readRDS(file.path(cd4_save_dir, "rda_obj", "seurat_processed.rds"))
seurat_cd8 <- readRDS(file.path(cd8_save_dir, "rda_obj", "seurat_processed.rds"))

find_clones <- function(seurat_object){
  seurat_object$cell_chains <- seurat_object$tcr_chains
  sep_columns <- c("chains", "cdr3", "cdr3_length",
                   "cdr3_nt_length", "v_gene", "d_gene", "j_gene", "c_gene",
                   "reads", "umis", "productive", "full_length",
                   "v_ins", "v_del", "v_mis", "d_ins", "d_del",
                   "d_mis", "j_ins", "j_del", "j_mis", "c_ins", "c_del", "c_mis",
                   "all_ins", "all_del", "all_mis", "vd_ins", "vd_del", "dj_ins",
                   "dj_del", "v_mis_freq", "d_mis_freq", "j_mis_freq",
                   "c_mis_freq", "all_mis_freq")
  sep_columns <- paste0("tcr_", sep_columns)
  keep_columns <- c("RNA_celltype_seurat", "RNA_cluster", "cell_chains")
  
  all_info <- seurat_object[[]] %>%
    dplyr::select(dplyr::all_of(c(sep_columns, keep_columns)), tcr_n_chains) %>%
    tibble::rownames_to_column("barcode") %>%
    dplyr::filter(!is.na(tcr_n_chains))
  
  all_info_split <- cSplit(all_info, sep_columns, sep = ";", direction = "long") %>%
    dplyr::filter(!is.na(tcr_chains))

  
  return_info_alpha <- find_clones_chain(all_info_split, chain = "TRA")
  return_info_beta <- find_clones_chain(all_info_split, chain = "TRB")
  
  return(list("alpha" = return_info_alpha,
              "beta" = return_info_beta))
}

determine_second_chain <- function(all_info, group_barcodes,
                                   new_group_barcodes, all_barcode_clones,
                                   group, new_group){
  # If the second chain is the same and doesn't have a different first
  # chain, we are fine. Otherwise, raise an error
  
  # So select all barcoes that have been found
  new_all_barcodes <- all_info %>%
    dplyr::filter(barcode %in% c(new_group_barcodes, group_barcodes)) %>%
    dplyr::group_by(barcode) %>%
    dplyr::select(barcode, group_name) %>%
    nest()
  
  new_all_barcode_clones <- unlist(new_all_barcodes$data)
  
  # If this is true, there are only two chains total. Name based on
  # the most common, but ask Maki
  if(all(unique(new_all_barcode_clones) %in% unique(all_barcode_clones))){
    if(length(group_barcodes) > length(new_group_barcodes)){
      renaming_list <- c(group, group)
      names(renaming_list) <- c(group, new_group)
    } else if(length(group_barcodes) < length(new_group_barcodes)){
      renaming_list <- c(new_group, new_group)
      names(renaming_list) <- c(group, new_group)
    } else {
      # If the number is the same, pick the smaller clone number
      new_group_num <- gsub("group", "", new_group)
      group_num <- gsub("group", "", group)
      if(as.integer(new_group_num) < as.integer(group_num)){
        renaming_list <- c(new_group, new_group)
        names(renaming_list) <- c(group, new_group)
      } else {
        renaming_list <- c(group, group)
        names(renaming_list) <- c(group, new_group)
      }
    }
  } else {
    # Mark as undetermined
    renaming_list <- c("undetermined", "undetermined")
    names(renaming_list) <- c(group, new_group)
  }
  return(renaming_list)
}

find_clones_chain <- function(all_info_split, chain){
  
  # First do just Alpha alone
  all_info <- all_info_split %>%
    dplyr::filter(tcr_chains == chain)
  
  
  # Right now I have one line per chain. What I want is 
  # 1. Find clone step 1 based on any clones that overlap
  # 2. Within each of those clones? Check if there are any cells with more than 1 
  # chain
  # 3. If there are cells with more than one chain, check if those second
  # chains match
  
  ## Find clones based on one chain only -----------------------------------------
  
  # Convert the data frame to a data.table
  setDT(all_info)
  
  # Add a new column with group names
  all_info[, group_name := paste0("group", .GRP), 
                 by = .(tcr_cdr3, tcr_v_gene, tcr_j_gene)]
  
  all_info <- all_info %>%
    dplyr::group_by(group_name) %>%
    dplyr::add_count(name = "group_count")
  
  # Find if there is a second chain ----------------------------------------------
  # First get names of clones that have more than one cell in them
  all_clones <- all_info %>%
    dplyr::filter(group_count > 1)
  
  # Goal
  # Each barcode should only be associated with clones
  # This means that any barcodes with two chains should be corrected to only
  # one group name.
  # If either of these chains are associated with a group, the two
  # groups should be merged
  # Current issues
  # Some of these will come up twice using this strategy
  clone_info <- lapply(unique(all_info$group_name), function(group){
    #print(group)
    # First subset to just the group
    group_clones <- all_info %>%
      dplyr::filter(group_name == group)
    
    # Then pull out all barcodes
    group_barcodes <- group_clones$barcode
    
    # Now pull out all information for those barcodes and see if there are
    # two chains in any
    all_barcodes <- all_info %>%
      dplyr::filter(barcode %in% group_barcodes) %>%
      dplyr::group_by(barcode) %>%
      dplyr::select(barcode, group_name) %>%
      nest()
    
    all_barcode_clones <- unlist(all_barcodes$data)
    
    # If there is another group found, then
    # 1. Make sure all cells in the second group were part of the first group
    # or all cells in the first group were part of the second (depending on
    # which is larger)
    # 2. If they were, then see which group has more cells and return the mapping
    # of renaming so all cells have the same group for their chains, if they
    # aren't, throw an error
    if(length(unique(all_barcode_clones)) == 2){
      new_group <- unique(all_barcode_clones)[unique(all_barcode_clones != group)]
      # First subset to just the group
      new_group_clones <- all_info %>%
        dplyr::filter(group_name == new_group)
      
      # Then pull out all barcodes
      new_group_barcodes <- new_group_clones$barcode
      
      if(length(group_barcodes) > length(new_group_barcodes)){
        if(all(new_group_barcodes %in% group_barcodes)){
          renaming_list <- c(group, group)
          names(renaming_list) <- c(group, new_group)
        } else {
          renaming_list <- determine_second_chain(all_info, group_barcodes,
                                                  new_group_barcodes,
                                                  all_barcode_clones,
                                                  group, new_group)
        }
      } else if(length(new_group_barcodes) > length(group_barcodes)){
        if(all(group_barcodes %in% new_group_barcodes)){
          renaming_list <- c(new_group, new_group)
          names(renaming_list) <- c(group, new_group)
        } else {
          renaming_list <- determine_second_chain(all_info, group_barcodes,
                                                  new_group_barcodes, 
                                                  all_barcode_clones,
                                                  group, new_group)
        }
      } else {
        if(all(group_barcodes %in% new_group_barcodes)){
          # If the number is the same, pick the smaller clone number
          new_group_num <- gsub("group", "", new_group)
          group_num <- gsub("group", "", group)
          if(as.integer(new_group_num) < as.integer(group_num)){
            renaming_list <- c(new_group, new_group)
            names(renaming_list) <- c(group, new_group)
          } else {
            renaming_list <- c(group, group)
            names(renaming_list) <- c(group, new_group)
          } 
        } else {
          renaming_list <- determine_second_chain(all_info, group_barcodes,
                                                  new_group_barcodes,
                                                  all_barcode_clones,
                                                  group, new_group)
        }
      }
    } else if(length(unique(all_barcode_clones)) == 1){
      renaming_list <- c(group)
      names(renaming_list) <- group
    } else {
      # Mark as undetermined
      group_names <- unique(all_barcode_clones)
      renaming_list <- rep("undetermined_first_second_chain_clash", times = 3)
      names(renaming_list) <- group_names
    }
    
    return(renaming_list)
  })
  
  
  clone_info <- unlist(clone_info)
  
  clone_info <- clone_info[!duplicated(names(clone_info))]
  
  all_info$group_name_updated <- clone_info[all_info$group_name]
  
  return(all_info)
  
}

cd4_clones <- find_clones(seurat_cd4)
cd8_clones <- find_clones(seurat_cd8)

# Subset based on UMI counts for the tcr reads
if(sample == "Npod6553_PLN"){
  # Remove any clones with low UMI count as they could be background
  cd4_clones <- lapply(cd4_clones, function(x){
    x <- x[x$tcr_umis > 1,]
    return(x)
  })
  
  cd8_clones <- lapply(cd8_clones, function(x){
    x <- x[x$tcr_umis > 1,]
    return(x)
  })
}

# Save data
save_wb <- openxlsx::createWorkbook()

cd4_tra <- cd4_clones$alpha %>%
  dplyr::select(barcode, tcr_cdr3, tcr_cdr3_length,
                tcr_v_gene, tcr_j_gene, tcr_c_gene,
                tcr_productive, RNA_celltype_seurat,
                group_name, group_name_updated)

cd4_trb <- cd4_clones$beta %>%
  dplyr::select(barcode, tcr_cdr3, tcr_cdr3_length,
                tcr_v_gene, tcr_j_gene, tcr_c_gene,
                tcr_productive, RNA_celltype_seurat,
                group_name, group_name_updated)

cd8_tra <- cd8_clones$alpha %>%
  dplyr::select(barcode, tcr_cdr3, tcr_cdr3_length,
                tcr_v_gene, tcr_j_gene, tcr_c_gene,
                tcr_productive, RNA_celltype_seurat,
                group_name, group_name_updated)

cd8_trb <- cd8_clones$beta %>%
  dplyr::select(barcode, tcr_cdr3, tcr_cdr3_length,
                tcr_v_gene, tcr_j_gene, tcr_c_gene,
                tcr_productive, RNA_celltype_seurat,
                group_name, group_name_updated)

openxlsx::addWorksheet(wb = save_wb, sheetName = "cd8_alpha")
openxlsx::writeData(wb = save_wb, sheet = "cd8_alpha", x = cd8_tra)

openxlsx::addWorksheet(wb = save_wb, sheetName = "cd8_beta")
openxlsx::writeData(wb = save_wb, sheet = "cd8_beta", x = cd8_trb)

openxlsx::addWorksheet(wb = save_wb, sheetName = "cd4_alpha")
openxlsx::writeData(wb = save_wb, sheet = "cd4_alpha", x = cd4_tra)

openxlsx::addWorksheet(wb = save_wb, sheetName = "cd4_beta")
openxlsx::writeData(wb = save_wb, sheet = "cd4_beta", x = cd4_trb)

openxlsx::saveWorkbook(wb = save_wb, file = file.path(save_dir, sample, "files",
                                                      "cd4_cd8_clones.xlsx"),
                       overwrite = TRUE)

