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

#sample <- args[[1]]
#sample <- gsub("__.*", "", sample)
sample <- "cd8"

#sample_info <- args[[3]]
sample_info <- here("files/sample_info.tsv")

#results_dir <- args[[2]]
results_dir <- here("results")

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
save_dir <- file.path(results_dir, "R_analysis", sample)

# Read in data
seurat_data <- readRDS(file.path(save_dir, "rda_obj", "seurat_processed.rds"))


seurat_data$cell_chains <- seurat_data$tcr_chains
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

all_info <- seurat_data[[]] %>%
  dplyr::select(dplyr::all_of(c(sep_columns, keep_columns)), tcr_n_chains) %>%
  tibble::rownames_to_column("barcode") %>%
  dplyr::filter(!is.na(tcr_n_chains))
# 
# 
# # First try only with n_chains = 2
# all_info <- all_info %>%
#   dplyr::filter(n_chains == 2)
# 
# all_data <- lapply(sep_columns, function(x){
#   sep_info <- all_info %>%
#     dplyr::select(dplyr::all_of(x)) %>%
#     dplyr::rename("cell" = dplyr::all_of(x)) %>%
#     tidyr::separate(cell, c("chain1", "chain2"), sep = ";")
#   
#   colnames(sep_info) <- paste(x, colnames(sep_info), sep = "_")
#   
#   return(sep_info)
# })


all_info_split <- cSplit(all_info, sep_columns, sep = ";", direction = "long") %>%
  dplyr::filter(!is.na(tcr_chains))

# First do just Alpha alone
all_info_alpha <- all_info_split %>%
  dplyr::filter(tcr_chains == "TRA")


# Right now I have one line per chain. What I want is 
# 1. Find clone step 1 based on any clones that overlap
# 2. Within each of those clones? Check if there are any cells with more than 1 
# chain
# 3. If there are cells with more than one chain, check if those second
# chains match

## Find clones based on one chain only -----------------------------------------

# Convert the data frame to a data.table
setDT(all_info_alpha)

# Add a new column with group names
all_info_alpha[, group_name := paste0("group", .GRP), 
               by = .(tcr_cdr3, tcr_v_gene, tcr_j_gene)]

all_info_alpha <- all_info_alpha %>%
  dplyr::group_by(group_name) %>%
  dplyr::add_count(name = "group_count")

# Find if there is a second chain ----------------------------------------------
# First get names of clones that have more than one cell in them
all_clones <- all_info_alpha %>%
  dplyr::filter(group_count > 1)

determine_second_chain <- function(all_info_alpha, group_barcodes,
                                   new_group_barcodes, all_barcode_clones,
                                   group, new_group){
  # If the second chain is the same and doesn't have a different first
  # chain, we are fine. Otherwise, raise an error
  
  # So select all barcoes that have been found
  new_all_barcodes <- all_info_alpha %>%
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

# Goal
# Each barcode should only be associated with clones
# This means that any barcodes with two chains should be corrected to only
# one group name.
# If either of these chains are associated with a group, the two
# groups should be merged
# Current issues
# Some of these will come up twice using this strategy
clone_info <- lapply(unique(all_info_alpha$group_name), function(group){
  #print(group)
  # First subset to just the group
  group_clones <- all_info_alpha %>%
    dplyr::filter(group_name == group)
  
  # Then pull out all barcodes
  group_barcodes <- group_clones$barcode
  
  # Now pull out all information for those barcodes and see if there are
  # two chains in any
  all_barcodes <- all_info_alpha %>%
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
    new_group_clones <- all_info_alpha %>%
      dplyr::filter(group_name == new_group)
    
    # Then pull out all barcodes
    new_group_barcodes <- new_group_clones$barcode
    
    if(length(group_barcodes) > length(new_group_barcodes)){
      if(all(new_group_barcodes %in% group_barcodes)){
        renaming_list <- c(group, group)
        names(renaming_list) <- c(group, new_group)
      } else {
        renaming_list <- determine_second_chain(all_info_alpha, group_barcodes,
                                                new_group_barcodes,
                                                all_barcode_clones,
                                                group, new_group)
      }
    } else if(length(new_group_barcodes) > length(group_barcodes)){
      if(all(group_barcodes %in% new_group_barcodes)){
        renaming_list <- c(new_group, new_group)
        names(renaming_list) <- c(group, new_group)
      } else {
        renaming_list <- determine_second_chain(all_info_alpha, group_barcodes,
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
        renaming_list <- determine_second_chain(all_info_alpha, group_barcodes,
                                                new_group_barcodes,
                                                all_barcode_clones,
                                                group, new_group)
      }
    }
  } else if(length(unique(all_barcode_clones)) == 1){
    renaming_list <- c(group)
    names(renaming_list) <- group
  } else {
    stop("Need to write a rule to deal with 3 groups found for a set of barcodes")
  }

  return(renaming_list)
})


clone_info <- unlist(clone_info)

clone_info <- clone_info[!duplicated(names(clone_info))]

all_info_alpha$group_name_updated <- clone_info[all_info_alpha$group_name]


# Things to check
# 1. The ones that changed
# * Did they change in a way I would expect? 
  # * expected patterns - should be places where the clone size was larger
  #   than 1 or that more than one alpha chain was identified.
# 2. The ones that didn't change
# * Should they have changed?
  # * expected patterns - should be places where the clone size was 1
  #   or that one alpha chain was identified or that the clone was present
  #   in the first set.

# Checking updated
updated <- all_info_alpha %>%
  dplyr::filter(group_name != group_name_updated)

# All but 1 has multiple TRA chains
table(updated$cell_chains)

table(updated$group_count)

not_updated <- all_info_alpha %>%
  dplyr::filter(group_name == group_name_updated)

# The number of multiplets here is the same
table(not_updated$cell_chains)


# This does get rid of any cells with multiple ra chains
not_updated_singles <- not_updated %>%
  dplyr::filter(!barcode %in% updated$barcode)

table(not_updated_singles$cell_chains)

table(not_updated_singles$group_count)

greater_group <- not_updated %>%
  dplyr::filter(group_count > 1)

# Only 4 have more than one TRA
table(greater_group$cell_chains)

greater_group %>%
  dplyr::filter(grepl("TRA;TRA", cell_chains)) %>%
  dplyr::select(group_name, group_name_updated, barcode)

updated %>%
  dplyr::filter(group_count > 1) %>%
  dplyr::select(group_name, group_name_updated, barcode)

# Check that these groups are no longer present
grep("group610", all_info_alpha$group_name_updated)
grep("group930", all_info_alpha$group_name_updated)

# Check the one barcode that is alone
all_info_alpha %>%
  dplyr::filter(barcode == "TGAGGGAGTTAAGAAC-1") %>%
  dplyr::select(barcode, group_name, group_name_updated)


# Everything looked good to me so I am somewhat confident that this is working
