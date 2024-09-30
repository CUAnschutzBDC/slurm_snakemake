# Goals
# Map the groupx to sample name
# Barplot with different splicing changes for each comparison
# Gene ontology for each comparison... look for something cool

library(tidyverse)
library(here)
library(gprofiler2)
library(KEGGREST)
library(org.Hs.eg.db)
library(ggraph)
library(ggdendro)

# Get group to sample mapping
mapping_file <- read.table(here("files/rmats_sample_info.tsv"), header = TRUE)

mapping_file <- mapping_file %>%
  dplyr::select(group, name) %>%
  dplyr::distinct()

# Get all combined files
all_dirs <- list.dirs(here("results/rmats_cutadapt_trim"), recursive = FALSE)
file_names <- c("combined_JCEC.txt", "combined_JC.txt")

make_full_file <- function(file_name, mapping_file, all_dirs){
  full_file <- lapply(all_dirs, function(full_dir){
    dir_name <- basename(full_dir)
    full_path <- file.path(full_dir, file_name)
    
    main_file <- read.table(full_path)
    # Split this to pull out the group names
    samples <- strsplit(dir_name, split = "_")
    # Use this and the mapping to get sample names
    # Add the sample names to the object
    main_file$sample_one <- mapping_file[
      mapping_file$group == samples[[1]][[1]], "name"
    ]
    
    main_file$sample_two <- mapping_file[
      mapping_file$group == samples[[1]][[2]], "name"
    ]
    
    main_file$comparison <- paste(main_file$sample_one, main_file$sample_two,
                                  sep = "_")
    
    return(main_file)
  })
  full_file <- do.call(rbind, full_file)
  return(full_file)
}

gost_plots <- function(results_table, source, title){
  source_results <- results_table[grep(source, results_table$source), ]
  if (nrow(source_results) > 0) {
    source_results <- source_results[order(source_results$precision), ]
    if(nrow(source_results) > 40){
      source_results <- source_results[1:40, ]
    }
    source_results$term_name <- factor(source_results$term_name, levels = 
                                         unique(source_results$term_name))
    source_results$log_padj <- -log10(source_results$p_value)
    go_plot <-ggplot(source_results, aes(x = precision,
                                         y = term_name,
                                         color = log_padj,
                                         size = intersection_size)) +
      geom_point() +
      theme_classic() +
      scale_size(name = "Intersection",
                 range(c(1,max(source_results$intersection_size)))) +
      viridis::scale_color_viridis() +
      theme(text = ggplot2::element_text(size = 10)) +
      ggtitle(paste0(source, ": ", title)) +
      labs(color = "-log10(p-value)") +
      xlab("Precision (proportion of genes)") +
      ylab("Term")
    
    
    
    return(go_plot)
  }
}

run_gost <- function(rmats_file,
                     comparison,
                     custom_bg = FALSE, 
                     correction_method = "gSCS",
                     exclude_iea = FALSE,
                     save_dir = NULL,
                     plot_dir = NULL,
                     evcodes = FALSE) {
  
  subset_file <- rmats_file[rmats_file$comparison == comparison,]
  gene_list <- unique(subset_file$geneSymbol)
  gost_res <- gost(query = gene_list,
                   organism = "hsapiens",
                   ordered_query = TRUE,
                   exclude_iea = exclude_iea,
                   user_threshold = 0.05,
                   correction_method = correction_method,
                   custom_bg = custom_bg,
                   evcodes = evcodes)
  
  if(evcodes){
    sig_gost_res <- gost_res$result[gost_res$result$significant == TRUE, ] %>%
      dplyr::select(!c(intersection, evidence_codes))
  } else {
    sig_gost_res <- gost_res$result[gost_res$result$significant == TRUE, ]
  }
  
  sig_gost_res$expected_num <- (sig_gost_res$query_size * 
                                  sig_gost_res$term_size)/ length(custom_bg)
  sig_gost_res$representation <- sig_gost_res$intersection_size/sig_gost_res$expected_num
  
  if (!is.null(save_dir)){
    save_name_csv <- paste0(comparison, ".csv")
    save_rds <- paste0(comparison, ".rds")
    save_excel <- paste0(comparison, ".xlsx")
    saveRDS(gost_res, file = file.path(save_dir, save_rds))
    if(nrow(sig_gost_res) > 1){
      sig_res_csv <- apply(sig_gost_res, 2, as.character)
      write.csv(sig_res_csv, file = file.path(save_dir, save_name_csv))
      
      excel_wb <- openxlsx::createWorkbook()
      for(type in unique(sig_gost_res$source)){
        write_data <- sig_gost_res[sig_gost_res$source == type,]
        sheet_name <- gsub(":", "_", type)
        openxlsx::addWorksheet(wb = excel_wb, sheetName = sheet_name)
        openxlsx::writeData(wb = excel_wb, sheet = sheet_name, x = write_data)
      }
      openxlsx::saveWorkbook(wb = excel_wb,
                             file = file.path(save_dir, save_excel),
                             overwrite = TRUE)
    }
  }
  # These make the plots
  plots <- list()
  plots$GOBP <- gost_plots(sig_gost_res, "GO:BP", comparison)
  plots$GOMF <- gost_plots(sig_gost_res, "GO:MF", comparison)
  plots$GOCC <- gost_plots(sig_gost_res, "GO:CC", comparison)
  plots$KEGG <- gost_plots(sig_gost_res, "KEGG", comparison)
  plots$TF <- gost_plots(sig_gost_res, "TF", comparison)
  
  if (!is.null(plot_dir)){
    save_name <- paste0(comparison, ".pdf")
    
    # This opens up a pdf file. We will save many images into this file
    pdf(file.path(plot_dir, save_name))
    
    plots_list <- lapply(plots, function(x){
      if (!is.null(x)){
        plot(x)
      }
    })
    
    # This closes the pdf
    dev.off()
    
    
  }
  
  return(sig_gost_res)
}

p_val <- 0.05
difference <- 0.1
save_dir <- here("results/R_analysis/images/rmats/")

full_files <- lapply(file_names, function(file){
  return_file <- make_full_file(
    file_name = file,
    mapping_file = mapping_file,
    all_dirs = all_dirs)
  
  return_file <- return_file[return_file$FDR < p_val & 
                               abs(return_file$IncLevelDifference) > difference,]
  
  save_name <- gsub("\\.txt", "", file)
  write.csv(return_file,
            file = file.path(save_dir, save_name, "rmats_res.csv"))
  
  return(return_file)
})

list_names <- gsub("\\.txt", "", file_names)

names(full_files) <- list_names

# Think about cutoffs, FDR < 0.05, IncLevelDifference?? 0.01?


# Make a barplot
plot_file <- full_files[[2]]

plot_file <- plot_file %>%
  dplyr::select(splice_type, comparison) %>%
  dplyr::group_by(splice_type, comparison) %>%
  dplyr::add_count(name = "splice_counts") %>%
  dplyr::ungroup() %>%
  dplyr::group_by(comparison) %>%
  dplyr::add_count(name = "comparison_counts") %>%
  dplyr::ungroup() %>%
  dplyr::distinct() %>%
  dplyr::mutate(percent_splice_type = splice_counts / comparison_counts * 100)

colors <- RColorBrewer::brewer.pal(n = 5, name = "Set1")
names(colors) <- unique(plot_file$splice_type)

# Count
plot_one <- ggplot2::ggplot(
  plot_file,
  ggplot2::aes(x = comparison, y = splice_counts, fill = splice_type)
) +
  ggplot2::geom_bar(position = "stack", stat = "identity") +
  ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
  ggplot2::scale_fill_manual(values = colors)

# Percent
plot_two <- ggplot2::ggplot(
  plot_file,
  ggplot2::aes(x = comparison, y = percent_splice_type, fill = splice_type)
) +
  ggplot2::geom_bar(position = "stack", stat = "identity") +
  ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
  ggplot2::scale_fill_manual(values = colors)

save_dir <- here("results/R_analysis/images/rmats/combined_JC/")

pdf(file.path(save_dir, "splice_type_counts.pdf"),
    width = 6, height = 6)

print(plot_one)

dev.off()

pdf(file.path(save_dir, "splice_type_percents.pdf"),
    width = 6, height = 6)

print(plot_two)

dev.off()


# Write a function to run GO on all 
# Count table
count_path <- here("results/nkx2.2_timecourse_countTable_cutadapt_trim.txt")
count_table <- read.table(count_path, sep = "\t",
                          row.names = 1, header = T)
keep <- rowSums(count_table) >= 10
count_table <- count_table[keep,]

# Get all genes expressed in the dataset
background_genes <- gsub("_ENS.*", "", rownames(count_table))

all_res <- lapply(names(full_files), function(file_type){
  save_dir <- here("results/R_analysis/images/rmats", file_type)
  
  ifelse(!dir.exists(save_dir), dir.create(save_dir, recursive = TRUE), FALSE)
  use_file <- full_files[[file_type]]
  all_gost <- lapply(unique(use_file$comparison), function(comparison){
    run_gost(rmats_file = use_file,
             comparison = comparison,
             custom_bg = background_genes, 
             correction_method = "gSCS",
             exclude_iea = FALSE,
             save_dir = save_dir,
             plot_dir = save_dir,
             evcodes = FALSE) 
  })

  names(all_gost) <- unique(use_file$comparison)
  return(all_gost)
  
})

names(all_res) <- list_names

test_file <- all_res$combined_JC
names(test_file) <- unique(full_files[[2]]$comparison)

new_test <- test_file$S5d1_S5d4

all_go <- lapply(names(test_file), function(comparision_name){
  new_test <- test_file[[comparision_name]]
  new_test <- new_test %>%
    dplyr::filter(source == "GO:BP") %>%
    dplyr::filter(representation > 2) %>%
    dplyr::filter(intersection_size > 5) %>%
    dplyr::arrange(desc(representation))  %>%
    dplyr::mutate(comparison = comparision_name)
})

all_kegg <- lapply(names(test_file), function(comparision_name){
  new_test <- test_file[[comparision_name]]
  
  new_test <- new_test %>%
    dplyr::filter(source == "KEGG") %>%
    dplyr::filter(representation > 2) %>%
    dplyr::filter(intersection_size > 5) %>%
    dplyr::arrange(desc(representation))  %>%
    dplyr::mutate(comparison = comparision_name)
})

all_go <- do.call(rbind, all_go)
all_kegg <- do.call(rbind, all_kegg)

terms <- c("neuron", "presynapse", "vesicle", "clathrin")
terms <- paste(terms, collapse = "|")

interesting_go <- all_go[grepl(terms, all_go$term_name),]

interesting_go <- interesting_go %>% 
  dplyr::select(p_value, intersection_size, term_name, representation, comparison) %>%
  dplyr::filter(comparison == "S5d1_S5d4") %>%
  dplyr::filter(representation > 2.9)

interesting_go$term_name <- factor(interesting_go$term_name,
                                   levels = interesting_go$term_name)

ggplot2::ggplot(interesting_go, ggplot2::aes(x = representation,
                                             y = term_name,
                                             size = intersection_size,
                                             color = -log10(p_value))) +
  ggplot2::geom_point() +
  ggplot2::scale_color_viridis_c(option = "inferno") +
  ggplot2::theme_classic() +
  ggplot2::xlab("Fold Enrichment") +
  ggplot2::xlim(0, 8) +
  ggplot2::ggtitle("S5d1 vs S5d4")


p1 <- ggplot2::ggplot(interesting_go, ggplot2::aes(x = representation,
                                             y = term_name,
                                             fill = -log10(p_value))) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::scale_fill_viridis_c(option = "inferno") +
  ggplot2::theme_classic() +
  ggplot2::xlab("Fold Enrichment") +
  ggplot2::xlim(0, 8) +
  ggplot2::ggtitle("S5d1 vs S5d4")

save_dir <- here("results/R_analysis/images/rmats/combined_JC/")
pdf(file.path(save_dir, "vesicle_transport_go.pdf"),
    width = 10, height = 8)
print(p1)
dev.off()

top_go <- all_go  %>% 
  dplyr::select(p_value, intersection_size, term_name, representation, comparison) %>%
  dplyr::filter(comparison == "S5d1_S5d4") %>%
  dplyr::filter(representation > 2.9)


#### Kegg clustering -----------------------------------------------------------

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}


make_histogram <- function(query, gse_res, p_val_cutoff = 0.01){
  
  gse_res <- gse_res[gse_res$source == "KEGG" & 
                       gse_res$query == query,]
  
  gse_list <- gse_res$term_id
  names(gse_list) <- gse_res$term_name
  
  path <- keggLink("pathway", "hsa")
  
  
  all_genes <- lapply(gse_list, function(x){
    search_term <- gsub("KEGG:", "", x)
    
    tryCatch({
      genes_in_term <- path[grepl(search_term, path)]
      
      all_gene_id <- keggConv("ncbi-geneid", names(genes_in_term))
      
      all_gene_id <- gsub("ncbi-geneid:", "", all_gene_id)
      
      gene_symbols <- mapIds(org.Hs.eg.db, keys = all_gene_id, 
                             keytype = "ENTREZID", column = "SYMBOL")  
      
    }, error = function(e) {
      cat("Error:", conditionMessage(e), "\n")
      gene_symbols <- "NA"
    })
    
    return(as.character(gene_symbols))
  })
  
  total_values <- length(all_genes) ^ 2
  final_jaccard_matrix <- Matrix::Matrix(data = rep(0, total_values),
                                         nrow = length(all_genes),
                                         ncol = length(all_genes),
                                         sparse = TRUE)
  
  for (x in 1:length(all_genes)){
    for (y in 1:length(all_genes)){
      
      save_val <- jaccard(a = all_genes[[x]],
                          b = all_genes[[y]])
      final_jaccard_matrix[x, y] <- save_val
      final_jaccard_matrix[y, x] <- save_val
    }
  }
  
  colnames(final_jaccard_matrix) <- names(all_genes)
  rownames(final_jaccard_matrix) <- names(all_genes)
  
  hc_full <- stats::hclust(stats::as.dist(1- final_jaccard_matrix),
                           method = "ward.D2")
  
  gse_small <- gse_res[gse_res$p_value < p_val_cutoff,]
  
  jaccard_small <- final_jaccard_matrix[rownames(final_jaccard_matrix)
                                        %in% gse_small$term_name,
                                        colnames(final_jaccard_matrix)
                                        %in% gse_small$term_name]
  
  hc_small <- stats::hclust(stats::as.dist(1- jaccard_small),
                            method = "ward.D2")
  
  return(list(full = hc_full, small = hc_small, gse_res = gse_res))
  
}

plot_histogram <- function(hc, gse_res, ylim = c(2, -2), title = NULL){
  # Convert hclust to dendrogram
  dend <- as.dendrogram(hc)
  
  
  ddata <- dendro_data(dend, type = "rectangle")
  
  dot_data <- merge(ddata$labels, gse_res, by.x = "label", by.y = "term_name",
                    all.x = TRUE, all.y = FALSE)
  
  
  p <- ggplot(segment(ddata)) + 
    ggplot2::geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
    ggplot2::geom_text(data = ddata$labels, ggplot2::aes(x = x, y = y - 0.06,
                                                         label = label),
                       hjust = 0, angle = 0, size = 3) +
    ggplot2::geom_point(data = dot_data, ggplot2::aes(x = x, y = y,
                                                      color = -log(p_value),
                                                      size = intersection_size)) +
    ggplot2::scale_color_gradient(low = "blue", high = "red") +
    ggplot2::coord_flip() + 
    ggplot2::scale_y_reverse(expand = c(0.2, 0)) +
    ggplot2::ylim(ylim) +
    ggplot2::theme(axis.line = element_blank(),
                   axis.text = element_blank(),
                   axis.ticks = element_blank(),
                   axis.title = element_blank())
  
  if(!is.null(title)){
    p <- p +
      ggplot2::ggtitle(title) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  }
  
  return(p)
  
}


all_hist <- lapply(unique(all_kegg$comparison), function(rmats_test){
  top_kegg <- all_kegg  %>% 
    dplyr::filter(comparison == rmats_test) %>%
    dplyr::filter(representation > 2) %>%
    dplyr::mutate(query = rmats_test)
  
  res <- make_histogram(query = rmats_test, gse_res = top_kegg, p_val_cutoff = 0.01)
  
  title <- str_split(rmats_test, pattern = "_")[[1]]
  title <- paste0(title[[1]], " vs ", title[[2]], " KEGG pathways")
  hist <- plot_histogram(hc = res$full, ylim = c(2, -1),
                         title = title,
                         gse_res = res$gse_res) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank()
    )
  
  kegg_pathways <- nrow(top_kegg)
  
  pdf(file.path(save_dir, paste0(rmats_test, "_kegg_clustering.pdf")),
      width = 12, heigh = kegg_pathways/4)
  print(hist)
  dev.off()
  
  return(hist)
  
})

