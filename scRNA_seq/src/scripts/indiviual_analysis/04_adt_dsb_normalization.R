library(Seurat)
library(dsb)
library(tidyverse)
library(here)
library(scAnalysisR)

args <- commandArgs(trailingOnly = TRUE)

sample <- args[[1]]
sample <- gsub("__.*", "", sample)
#sample <- "Npod6456_PLN"

sample_info <- args[[3]]
#sample_info <- here("files/sample_info.tsv")

results_dir <- args[[2]]
#results_dir <- here("results")

sample_info <- read.table(sample_info, fill = TRUE, header = TRUE)

sample_info <- sample_info[sample_info$sample == sample,]

HTO <- sample_info$HTO
ADT <- sample_info$ADT
hash_ident <- sample_info$hash_ident

# Set directories
save_dir <- file.path(results_dir, "R_analysis", sample)

# Load in raw data
sample_path <- file.path(results_dir, sample, "outs", "per_sample_outs", 
		                         sample, "count", "raw_feature_bc_matrix")

sample_data <- Read10X(data.dir = sample_path)
if(ADT){
  if(HTO){
    
    # Load in HTOs and demultiplex
    protein_data <- sample_data[["Antibody Capture"]]
    hashtag_data <- protein_data[grepl("Hashtag", rownames(protein_data)), ]
    ADT_data <- protein_data[!grepl("Hashtag", rownames(protein_data)), ]
    sample_object <- CreateSeuratObject(counts = hashtag_data,
                                        assay = "HTO")
    
    sample_object[["ADT"]] <- CreateAssayObject(counts = ADT_data)
    
    sample_object <- NormalizeData(sample_object, assay = "HTO",
                                   normalization.method = "CLR")
    
    sample_object <- NormalizeData(sample_object, assay = "ADT",
                                   normalization.method = "CLR")
    
    sample_object_negative <- subset(sample_object, nCount_HTO == 0)
    
    sample_object_test <- subset(sample_object, nCount_HTO > 0)
    
    
    
    # Demultiplex
    sample_object_test <- HTODemux(sample_object_test, assay = "HTO",
                                   positive.quantile = 0.99, verbose = TRUE)
    
    sample_object_negative$HTO_classification.global <- "Negative"
    
    # This gave WAY too many positive cells so I'm just using the seurat object as the positive
    print(table(sample_object_test$HTO_classification))
    
    sample_object <- merge(sample_object_negative, sample_object_test, merge.data = TRUE)
    
    Idents(sample_object) <- "HTO_classification.global"
    
    # Create an object of just negative cells
    negative_object <- subset(sample_object, idents = "Negative")
    
    # Read in my exisitng object as the "positive"
    positive_object <- readRDS(file.path(save_dir, "rda_obj", "seurat_doublet.rds"))
    
    # Make sure there are no overlapping cells
    length(intersect(rownames(negative_object), rownames(positive_object)))
    
    neg_adt_matrix <- GetAssayData(negative_object, assay = "ADT",
                                   slot = 'counts') %>% as.matrix()
    positive_adt_matrix <- GetAssayData(positive_object, assay = "ADT",
                                        slot = 'counts') %>% as.matrix()
    
  } else {
    # Read in my exisitng object as the "positive"
    starting_object <- readRDS(file.path(save_dir, "rda_obj", "seurat_unfilt.rds"))
    
    # define a vector of cell-containing barcodes and remove them from unfiltered data 
    stained_cells <- colnames(starting_object)
    background <- setdiff(colnames(sample_data$`Gene Expression`), stained_cells)
    
    # split the data into separate matrices per assay 
    prot <- sample_data$`Antibody Capture`
    rna <- sample_data$`Gene Expression`
    
    # create metadata of droplet QC stats used in standard scRNAseq processing
    rna_size <- log10(Matrix::colSums(rna))
    prot_size <- log10(Matrix::colSums(prot))
    ngene <- Matrix::colSums(rna > 0)
    mtgene <- grep(pattern = "^MT-", rownames(rna), value = TRUE)
    propmt <- Matrix::colSums(rna[mtgene, ]) / Matrix::colSums(rna)
    md <- as.data.frame(cbind(propmt, rna_size, ngene, prot_size))
    md$bc <- rownames(md)
    md$droplet_class <- ifelse(test = md$bc %in% stained_cells, yes = 'cell',
                               no = 'background')
    
    # filter barcodes to only include those with data for both assays 
    md <- md %>% dplyr::filter(rna_size > 0 & prot_size > 0 )
    
    # ggplot(md, aes(y = rna_size, x = prot_size )) +
    #   theme_bw() +
    #   geom_bin2d(bins = 300) +
    #   scale_fill_viridis_c(option = "C") +
    #   facet_wrap(~droplet_class)
    
    background_only <- md %>% 
      dplyr::filter(droplet_class == "background") 
    upper <- quantile(background_only$prot_size, probs = 0.95)
    lower <- quantile(background_only$prot_size, probs = 0.05)
    
    ## CHANGE BASED ON THE ABOVE PLOT!!! ##
    
    # define a vector of background droplet barcodes based on protein library size and mRNA content
    # background_drops <- md[md$prot_size > 0.5 &
    #                          md$prot_size < 1.5 & md$rna_size < 2.75, ]$bc
    
    background_drops <- md[md$prot_size < upper & md$prot_size > lower,]$bc
    
    neg_adt_matrix <- as.matrix(prot[ , background_drops])
    
    # Get positive cells
    positive_object <- readRDS(file.path(save_dir, "rda_obj/seurat_doublet.rds"))
    
    positive_cells <- colnames(positive_object)
    positive_adt_matrix <- as.matrix(prot[ , positive_cells])
    
  }
  
  # Retain the original object
  positive_object[["CLR_ADT"]] <- positive_object[["ADT"]]
  
  # calculate quantiles of the raw protein matrix 
  d1 <- data.frame(pmax = apply(positive_adt_matrix, 1, max)) %>% 
    rownames_to_column('prot') %>% arrange(pmax)
  
  remove_genes <- d1[d1$pmax < 100,]$prot
  
  # CHECK D1 AND DETERMINE IF ANYTHING NEEDS TO BE REMOVED!!!
  # Think about removing CD27.1?
  prot_names <- rownames(positive_adt_matrix)
  positive_adt_matrix <- positive_adt_matrix[!prot_names %in% remove_genes, ]
  neg_adt_matrix <- neg_adt_matrix[!prot_names %in% remove_genes, ]
  
  # Fix insulin names
  rownames(neg_adt_matrix) <- gsub("INS_tet[a|b]", "INS_tet", 
                                   rownames(neg_adt_matrix))
  
  rownames(positive_adt_matrix) <- gsub("INS_tet[a|b]", "INS_tet", 
                                   rownames(positive_adt_matrix))
  
  # Run the normalization
  dsb_norm_prot <- DSBNormalizeProtein(
    cell_protein_matrix = positive_adt_matrix,
    empty_drop_matrix = neg_adt_matrix,
    denoise.counts = FALSE, # Set to TRUE if isotype controls are present
    use.isotype.control = FALSE)
  
  # Add to object
  positive_object[["DSB_ADT"]] <- CreateAssayObject(data = dsb_norm_prot)
  
  DefaultAssay(positive_object) <- "DSB_ADT"
  plot_adts <- rownames(positive_object)
  violins <- featDistPlot(positive_object, sep_by = "RNA_cluster",
                          geneset = plot_adts, assay = "DSB_ADT",
                          combine = FALSE) 

  DefaultAssay(positive_object) <- "RNA"  
} else {
  # Read in data
  positive_object <- readRDS(file.path(save_dir, "rda_obj", "seurat_doublet.rds"))
}

# Save
saveRDS(positive_object, file.path(save_dir, "rda_obj", "seurat_adt.rds"))
