#' Plots the enrichment values of selected TFs.
#' @param enrich_motifs motif_enrich output dataframe from findEnrichMotifPair
#' @param tfs TFs for which to plot enrichment
#' @param metric Value to plot; either pval or logFC
#' @return a ggplot2 plot object
#' @examples
#' \dontrun{
#' motifs_df <- results$motif_enrich
#' plotEnrichment(enrich_motifs = motifs_df, tfs = c("SOX2", "NANOG", "MYC", "POU5F1"),
#' metric = "pval")
#' }
#' @importFrom dplyr %>%
#' @importFrom ggplot2 ggplot aes_string geom_tile xlab ylab theme theme_bw 
#' element_blank element_text scale_fill_gradient
#' @export
#' 
plotEnrichment = function(enrich_motifs, 
                          tfs
                          ){
  
  # filter for these selected TFs
  enrich_motifs_filtered <- enrich_motifs %>% 
    dplyr::filter(TF_name %in% tfs) %>% 
    dplyr::distinct(TF_name, .keep_all = TRUE)
  
  
  # replace very low values so that it is easy to visualize
  enrich_motifs_filtered[enrich_motifs_filtered < 1e-100] <- 1e-100

  enrich_motifs_filtered$log10FC = log10(enrich_motifs_filtered$fold_enrich)
  
  enrich_motifs_filtered = dplyr::arrange(enrich_motifs_filtered, log10FC)
  enrich_motifs_filtered$TF_name <- factor(enrich_motifs_filtered$TF_name, 
                                           levels = enrich_motifs_filtered$TF_name)
  
  # choose color palette
  col_pal <- RColorBrewer::brewer.pal(9,"Blues")
  # enriched motifs heatmap using ggplot2
  ggplot(data = enrich_motifs_filtered, aes_string(x = NA, y = "TF_name", 
                                                  fill = "log10FC"))  + 
    geom_tile() + scale_fill_gradient(low = col_pal[1], high = col_pal[9]) + 
    ylab("") + xlab("") + theme_bw() + 
    theme(plot.background = element_blank()
          ,panel.grid.major = element_blank()
          ,panel.grid.minor = element_blank()
          ,text = element_text(size=16), axis.text.x = element_blank()
    )
  
}

#' Plots the enrichment values of selected TFs.
#' @param enrich_pair motif_pair_enrich output dataframe from findEnrichMotifPair
#' @param tfs TFs for which to plot enrichment
#' @param metric Value to plot; either pval or logFC
#' @return a ggplot2 plot object
#' @examples
#' \dontrun{
#' motif_pairs_df <- results$motif_pair_enrich
#' plotEnrichment(enrich_motifs = motif_pairs_df, tfs = c("SOX2", "NANOG", "MYC", "POU5F1"),
#' metric = "pval")
#' }
#' @importFrom dplyr %>%
#' @importFrom ggplot2 ggplot aes_string geom_tile xlab ylab theme theme_bw 
#' element_blank element_text scale_fill_gradient
#' @export
#' 
plotEnrichPair = function(enrich_pairs, 
                          tfs
                          ){

  # select top 10 binding partners and remove redundant motifs
  enrich_pairs_filtered <- enrich_pairs %>% 
    dplyr::filter(TF_name_1 %in% tfs) %>% 
    dplyr::group_by(TF_name_1, TF_name_2) %>% 
    dplyr::distinct(TF_name_2, .keep_all = TRUE) %>% 
    dplyr::group_by(TF_name_1) %>% dplyr::top_n(-10, pval_adj)
  
  enrich_pairs_filtered$log10FC = log10(enrich_pairs_filtered$fold_enrich)
 
  order = aggregate(enrich_pairs_filtered$log10FC, list(enrich_pairs_filtered$TF_name_1), FUN=mean) %>% 
          dplyr::arrange(-x)
  order_var = order$Group.1
  
  enrich_pairs_filtered_mat <- enrich_pairs_filtered %>%
    dplyr::select(TF_name_1, TF_name_2, log10FC) %>%
    tidyr::spread(TF_name_2, log10FC) %>% 
    tibble::column_to_rownames(var = "TF_name_1") %>%
    as.matrix() 
  
  #idx <- match(enrich_pairs_filtered$TF_name, rownames(enrich_pairs_filtered_mat))
  #enrich_pairs_filtered_mat <- na.omit(enrich_pairs_filtered_mat[idx, ])
  
  enrich_pairs_filtered_mat[is.na(enrich_pairs_filtered_mat)] <- 0
  col_pal <- c("white", RColorBrewer::brewer.pal(9,"Reds"))
  heatmap = pheatmap::pheatmap(enrich_pairs_filtered_mat, cluster_rows=F, cluster_cols=T, silent = T)
  cluster_order = heatmap$tree_col$order
  enrich_pairs_filtered_mat <- enrich_pairs_filtered_mat[order_var, cluster_order]
  
  if (length(tfs) > 1){
    enrich_pairs_melt = suppressWarnings(reshape2::melt(enrich_pairs_filtered_mat))
  } else {
    enrich_pairs_melt = data.frame(tf1 = rep(tfs, length(enrich_pairs_filtered_mat)),
                                   tf2 = names(enrich_pairs_filtered_mat),
                                   value = as.numeric(enrich_pairs_filtered_mat))
  }
  
  colnames(enrich_pairs_melt) = c("TF_name_1", y = "TF_name_2", "log10FC")
  ggplot(data = enrich_pairs_melt,
               aes_string(x = "TF_name_2", y = "TF_name_1", fill = "log10FC")) +
    geom_tile() + scale_fill_gradient(low = col_pal[1], high = col_pal[9]) +
    ylab("") + xlab("") + theme_bw() +
    theme(plot.background = element_blank()
          ,panel.grid.major = element_blank()
          ,panel.grid.minor = element_blank()
          ,text = element_text(size=16), axis.text.x = element_text(angle=90, vjust=0.5)
    )
  
  #return(gg)
}

#' Plots network of TF enrichment interactions.
#' @param enrich_pair motif_pair_enrich output dataframe from findEnrichMotifPair
#' @param tfs TF for which to plot enrichment network
#' @param color_TF Color of input TF
#' @param color_bind_TF Color of binding partner TFs
#' @return a ggplot2 plot object
#' @examples
#' \dontrun{
#' motif_pairs_df <- results$motif_pair_enrich
#' plotNetwork(enrich_pairs, TF_name = "SP8")
#' }
#' @importFrom dplyr %>%
#' @export
#' 
plotNetwork <- function(enrich_pairs, 
                        TF_name, 
                        color_TF = "#70d9e0", 
                        color_bind_TF = "#e841da"){
  
  
  enrich_pairs_filtered <- enrich_pairs %>% 
    dplyr::filter(TF_name_1 == TF_name) %>% 
    dplyr::group_by(TF_name_1, TF_name_2) %>% 
    dplyr::distinct(TF_name_2, .keep_all = TRUE) %>% 
    dplyr::group_by(TF_name_1) %>% dplyr::top_n(-10, pval_adj)
  
  edges <- enrich_pairs_filtered %>% 
    dplyr::filter(TF_name_1 == TF_name) %>% 
    dplyr::mutate(sig = -log(pval_adj)) %>% 
    dplyr::mutate(sig = sig*8/max(sig)) %>% as.data.frame() %>% 
    dplyr::select(TF_name_1, TF_name_2, sig)
  
  n_edges = length(unique(c(edges$TF_name_1, edges$TF_name_2)))
  edge_df_n <- n_edges - 1
  nodes <- data.frame(
    name=unique(c(edges$TF_name_1, edges$TF_name_2)),
    role=c(rep("TF", 1), rep("partner", edge_df_n))
  )
  # Turn it into igraph object
  network <- igraph::graph_from_data_frame(d=edges, vertices=nodes, 
                                           directed=FALSE) 
  
  col <- c(color_TF, color_bind_TF)
  # Create a vector of color
  my_color <- col[as.numeric(as.factor(igraph::V(network)$role))]
  # plotting the network
  sig = igraph::E(network)$sig
  GGally::ggnet2(network, node.size = 12, node.color=my_color, edge.color = "grey",
         label = T, edge.size=log(sig), label.size = 4, layout.exp = 0.3)
    
}



