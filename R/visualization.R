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
                          tfs,
                          metric = "pval"){
  
  if (metric == "pval"){
    metric = "pval_adj"
  }
  
  # filter for these selected TFs
  enrich_motifs_filtered <- enrich_motifs %>% 
    dplyr::filter(TF_name %in% tfs) %>% 
    dplyr::distinct(TF_name, .keep_all = TRUE)
  
  
  # replace very low values so that it is easy to visualize
  enrich_motifs_filtered[enrich_motifs_filtered < 1e-100] <- 1e-100
  
  enrich_motifs_filtered$logp = -log(enrich_motifs_filtered$pval_adj)
  enrich_motifs_filtered$logFC = log(enrich_motifs_filtered$fold_enrich)
  
  if (metric == "pval_adj"){
    hm_value = "logp"
  } else {hm_value = "logFC"}
  
  enrich_motifs_filtered$x <- NA
  # choose color palette
  col_pal <- RColorBrewer::brewer.pal(9,"Blues")
  # enriched motifs heatmap using ggplot2
  gg <- ggplot(data = enrich_motifs_filtered, aes_string(x = NA, y = "TF_name", 
                                                  fill = hm_value))  + 
    geom_tile() + scale_fill_gradient(low = col_pal[1], high = col_pal[9]) + 
    ylab("") + xlab("") + theme_bw() + 
    theme(plot.background = element_blank()
          ,panel.grid.major = element_blank()
          ,panel.grid.minor = element_blank()
          ,text = element_text(size=16), axis.text.x = element_blank()
    )
  return(gg)
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
                          tfs,
                          metric = "pval"){
  
  if (metric == "pval"){
    metric = "pval_adj"
  }
  
  # select top 10 binding partners and remove redundant motifs
  enrich_pairs_filtered <- enrich_pairs %>% 
    dplyr::filter(TF_name_1 %in% tfs) %>% 
    group_by(TF_name_1, TF_name_2) %>% 
    dplyr::distinct(TF_name_2, .keep_all = TRUE) %>% 
    dplyr::group_by(TF_name_1) %>% top_n(-10, pval_adj)
  
  
  # choose color palette
  col_pal <- RColorBrewer::brewer.pal(9,"Blues")
  # enriched motif pairs heaatmap using ggplot2
  enrich_pairs_filtered$logp = -log(enrich_pairs_filtered$pval_adj)
  enrich_pairs_filtered$logFC = log(enrich_pairs_filtered$fold_enrich)
  
  if (metric == "pval_adj"){
    hm_value = "logp"
  } else {hm_value = "logFC"}
  
  gg <- ggplot(data = enrich_pairs_filtered, 
               aes_string(x = "TF_name_2", y = "TF_name_1", fill = hm_value)) + 
    geom_tile() + scale_fill_gradient(low = col_pal[1], high = col_pal[9]) + 
    ylab("") + xlab("") + theme_bw() + 
    theme(plot.background = element_blank()
          ,panel.grid.major = element_blank()
          ,panel.grid.minor = element_blank()
          ,text = element_text(size=16), axis.text.x = element_text(angle=90, vjust=0.5)
    )
  
  return(gg)
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
  enrich_pairs = results$motif_pair_enrich
  TF_name = "BACH2"
  enrich_pairs_filtered <- enrich_pairs %>% 
    dplyr::filter(TF_name_1 == TF_name) %>% 
    group_by(TF_name_1, TF_name_2) %>% 
    dplyr::distinct(TF_name_2, .keep_all = TRUE) %>% 
    dplyr::group_by(TF_name_1) %>% top_n(-10, pval_adj)
  
  edges <- enrich_pairs_filtered %>% 
    dplyr::filter(TF_name_1 == TF_name) %>% 
    dplyr::mutate(sig = -log(pval_adj)) %>% 
    dplyr::mutate(sig = sig*8/max(sig)) %>% as.data.frame() %>% 
    dplyr::select(TF_name_1, TF_name_2, sig)
  
  n_edges = length(unique(c(edges$TF_name_1, edges$TF_name_2)))
  nodes <- data.frame(
    name=unique(c(edges$TF_name_1, edges$TF_name_2)),
    role=c(rep("TF",1),rep("partner", n_edges - 1))
  )
  # Turn it into igraph object
  network <- igraph::graph_from_data_frame(d=edges, vertices=nodes, 
                                           directed=FALSE) 
  
  col <- c(color_TF, color_bind_TF)
  # Create a vector of color
  my_color <- col[as.numeric(as.factor(igraph::V(network)$role))]
  # plotting the network
  plot(network, vertex.color=my_color, vertex.shape = c("rectangle"), 
       vertex.size = 40, vertex.size2 = 30, vertex.label.cex=0.8, 
       vertex.label.color="black", edge.width=igraph::E(network)$sig, 
       edge.color="black")
}

install.packages('extrafont')
library(extrafont)
font_addpackage(pkg = "serif")
fonts_installed <- fonts()
df_font <- fonttable()
loadfonts()
font_paths()






