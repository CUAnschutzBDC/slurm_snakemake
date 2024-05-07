# Functions --------------------------------------------------------------------
plot_pca <- function(vsd, group_by, color_palette = NULL){
  if(is.null(color_palette)){
    color_palette <- RColorBrewer::brewer.pal(length(levels(vsd[[group_by]])), "Set1")
    names(color_palette) <- levels(vsd[[group_by]])
  }
  pcaData <- plotPCA(vsd, intgroup = group_by, returnData = T)
  colnames(pcaData)[3] <- "group_by"
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  pca_plot <- ggplot(data = pcaData,
                     mapping = aes(x = PC1,
                                   y = PC2,
                                   color = group_by)) +
    theme_classic() +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    labs(color = group_by) +
    #coord_fixed() +
    scale_color_manual(values = color_palette)
  return(pca_plot)
}

# Make gene plots
plot_genes <- function(gene_id, gene_id_list, deseq_obj,
                       intgroup, plot_ggplot = TRUE,
                       color = NULL, return_data = TRUE,
                       print = TRUE, save_path = NULL,
                       normalized = FALSE, vst = NULL){
  if(normalized & is.null(vst)){
    stop(paste0("If wanted normalized values, you must provide a vst object ",
                "where you performed limma::removeBatchEffect"))
  }
  
  # Locate the index of the gene of interest
  index <- grep(paste0("^", gene_id, "$"), gene_id_list)
  if(length(index) == 0){
    print(paste0(gene_id, " not in deseq object"))
    return(NULL)
  } else if(length(index == 1)) {
    counts_plot <- make_plots(index = index, deseq_obj = deseq_obj,
                              intgroup = intgroup, plot_ggplot = plot_ggplot,
                              color = color, return_data = return_data,
                              print = print, save_path = save_path,
                              normalized = normalized, vst = vst)
    return(counts_plot)
  } else {
    # This is in case a gene has two ensembl values
    plot_list <- lapply(index, function(x) make_plots(index = x,
                                                      deseq_obj = deseq_obj,
                                                      intgroup = intgroup,
                                                      plot_ggplot = plot_ggplot,
                                                      color = color, 
                                                      return_data = return_data,
                                                      print = print,
                                                      save_path = save_path,
                                                      normalized = normalized,
                                                      vst = vst))
    return(plot_list)
  }
}
make_plots <- function(index, deseq_obj, intgroup, plot_ggplot,
                       color, return_data, print, save_path,
                       normalized = FALSE, vst = NULL){
  gene_name <- rownames(deseq_obj)[index]
  if(normalized){
    mat <- assay(vst)
    counts_plot <- t(mat[gene_name, , drop = FALSE]) %>% data.frame()
    colnames(counts_plot) <- "count"
    
    counts_plot$Group <- deseq_obj[[intgroup]]
    
  } else {
    counts_plot <- DESeq2::plotCounts(deseq_obj, gene = gene_name,
                                      intgroup = intgroup,
                                      returnData = TRUE) 
    colnames(counts_plot) <- c("count", "Group")
  }
  
  if(plot_ggplot){
    if(is.null(color)){
      color <- RColorBrewer::brewer.pal(length(levels(deseq_obj[[intgroup]])), "Set1")
      names(color) <- levels(deseq_obj[[intgroup]])
    }
    
    ggplot_counts_plot <- ggplot2::ggplot(counts_plot,
                                          ggplot2::aes(x = Group,
                                                       y = count)) +
      ggplot2::geom_point(ggplot2::aes(color = Group), size = 3) +
      ggplot2::scale_color_manual(values = color, name = "Group") +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                         vjust = 1,
                                                         hjust=1)) +
      ggplot2::ggtitle(gene_name)
    if(print){
      print(ggplot_counts_plot)
    }
    if (!is.null(save_path)){
      file_name <- paste0("counts_plot_", gene_name, ".pdf")
      file_path <- file.path(save_path, file_name)
      ggplot2::ggsave(filename = file_path, plot = ggplot_counts_plot)
    }
    if(return_data){
      return(ggplot_counts_plot)
    }
  } else {
    return(counts_plot)
  }
}



make_corrected_heatamp <- function(vsd, gene_list, meta_ave, meta_col,
                                   plot_meta_col = FALSE,
                                   max_val = 2.5,
                                   min_val = -2.5,
                                   cluster_rows = TRUE,
                                   cluster_cols = FALSE,
                                   coloring = NULL,
                                   plot_rownames = FALSE, ...){
  
  # Pull corrected values from vsd
  heatmap_df <- assay(vsd)[rownames(assay(vsd)) %in% gene_list,]
  cluster_order <- levels(meta_ave[[meta_col]])
  heatmap_scale <- t(scale(t(heatmap_df), scale = TRUE))
  blueYellow <- c("#352A86", "#343DAE", "#0262E0", "#1389D2", 
                  "#2DB7A3", "#A5BE6A", "#F8BA43", "#F6DA23", "#F8FA0D")
  
  meta_ave <- meta_ave[order(match(meta_ave[[meta_col]], 
                                   cluster_order)), , drop = FALSE]
  if(!all(rownames(meta_ave) %in% colnames(heatmap_scale))){
    heatmap_scale <- heatmap_scale[, colnames(heatmap_scale) %in%
                                     rownames(meta_ave)]
  }
  
  if (!identical(colnames(heatmap_scale), rownames(meta_ave))) {
    heatmap_scale <- heatmap_scale[, rownames(meta_ave)]
  }
  
  if (!plot_meta_col) {
    meta_ave[[meta_col]] <- NULL
  }
  heatmap_scale <- ifelse(heatmap_scale > max_val, max_val, 
                          heatmap_scale)
  heatmap_scale <- ifelse(heatmap_scale < min_val, min_val, 
                          heatmap_scale)
  heatmap <- pheatmap::pheatmap(heatmap_scale, cluster_rows = cluster_rows, 
                                cluster_cols = cluster_cols, 
                                show_rownames = plot_rownames, 
                                show_colnames = FALSE, 
                                annotation_col = meta_ave, 
                                annotation_colors = coloring, 
                                color = blueYellow, border_color = NA, 
                                clustering_method = "complete",
                                silent = TRUE, ...)
  
  return(heatmap)
  
}
