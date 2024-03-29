#' Finds enriched motif pairs in given a
#' genomic regions.
#'
#' @param target_data A data frame containing peaks coordinates with chr, start,
#' and end as first three columns
#' @param background_data Same as "target_data"
#' @param genome_ver Genome version, it should be either "hg19" or "hg38"
#' @param scramble_data A logical value (TRUE/FALSE) whether background data
#' need to be generated by scrambling target_data, by default FALSE
#' @param Pvalue_computation A character representing P value computations are
#' either based on binomial or hyper geometric distributions, it should be one
#' of c("binom", "hyper"), default = "binom"
#' @param Pvalue_threshold A numeric value for filtering the enriched motifs and
#' motif pairs
#' @param Pvalue_adjust_method A character representing method for adjustig the
#' P values for multple comparisons, it should be one or more of
#' c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
#' default = "BH"
#' @return A list with two data frames containing enriched motifs and motif
#' pairs with P-values, fold-enrichment and annotations
#' @examples
#' \dontrun{
#' findEnrichMotifPair_Jolma2015(target_data, background_data = NULL,
#' genome_ver = "hg38", scramble_data = TRUE)
#' }
#' @import stats
#' @importFrom motifmatchr matchMotifs
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom rlang .data
#' @importFrom SummarizedExperiment assays
#' @importFrom dplyr inner_join
#' @importFrom dplyr %>%
#' @importFrom assertthat assert_that
#' @export
#'
#'
findEnrichMotifPair_Jolma2015 = function(target_data, background_data = NULL,
                               genome_ver = "hg38", scramble_data = TRUE,
                               Pvalue_computation = "hyper",
                               Pvalue_threshold = 0.05,
                               Pvalue_adjust_method = "BH"){
    # checking input data and parameters
    assertthat::assert_that(is.data.frame(target_data))
    assertthat::assert_that(is.character(genome_ver))
    assertthat::assert_that(is.logical(scramble_data))
    assertthat::assert_that(is.character(Pvalue_computation))
    assertthat::assert_that(is.numeric(Pvalue_threshold))
    assertthat::assert_that(is.character(Pvalue_adjust_method))

    # loading the right genome based on user choice
    if(genome_ver == "hg19"){
        genome_data <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
    } else {genome_data <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38}

    # checking if background data is provided by user
    if(is.null(background_data) & !scramble_data){
        bg_gr <- generateBackgroundSeqs(input_data = target_data,
                                        genome_ver = genome_ver)
        message("background data is generated by matched GC, length and chromosome distribution")
    } else {message("background data is provided by the user or generated by scrambling the input sequences")}
    
    if(is.null(background_data) & scramble_data){
        bg_gr <- generateScrambleSeqs(input_data = target_data,
                                      genome_ver = genome_ver)
        message("background data is generated by scrambling the input sequences")
    }
    
    # convert data into GRanges objects
    if(!is.null(background_data) & !scramble_data){
        bg_gr <- background_data %>% GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE, ignore.strand = TRUE)
    }
    
    
    message("Check point 1 passed")
    tg_gr <- target_data %>% GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE, ignore.strand = TRUE)
    message("Check point 2 passed")
    # extract relevant motifs PWM and theri annotation based on motif database
    filter_col <- function(df1, col_name, val) {
        col_name <- dplyr::enquo(col_name) # captures the environment in which the function was called
        df1 %>% dplyr::filter((!!col_name) %in% val)
    }

    filter_col2 <- function(df1, col_name, val) {
        col_name <- dplyr::enquo(col_name) # captures the environment in which the function was called
        df1 %>% dplyr::filter((!!col_name) == val)
    }
    motif_anno <- Jolma_2015_motifs_pairs_anno
    #motif_anno <- base::subset(motif_anno_data, motif_source == motif_database)
    message("Check point 3 passed")
    #idx <- base::match(motif_anno[,1], names(Jolma_2015_motifs_pairs))
    message("Check point 4 passed")
    motifs_PWMs <- Jolma_2015_motifs_pairs
    message("Check point 5 passed")

    # find motif locations
    motif_ix_scores_tg <- motifmatchr::matchMotifs(motifs_PWMs, tg_gr, out = "scores", genome = genome_data)
    motif_ix_scores_bg <- motifmatchr::matchMotifs(motifs_PWMs, bg_gr, out = "scores", genome = genome_data)
    message("Check point 6 passed")

    f1 <- SummarizedExperiment::assays(motif_ix_scores_tg)$motifMatches
    f2 <- SummarizedExperiment::assays(motif_ix_scores_bg)$motifMatches
    # save(f1, f2, file = "tt.rda")
    # compute P-value and Fold enrichment
    motif_pval <- function(tg2 = f1, bg2 = f2, Pval_dist = Pvalue_computation){
        tg_mf_c1 <- data.frame(motif_pair_name = colnames(tg2), tg_motif_count = Matrix::colSums(tg2), row.names = NULL, stringsAsFactors = FALSE)
        bg_mf_c1 <- data.frame(motif_pair_name = colnames(bg2), bg_motif_count = Matrix::colSums(bg2), row.names = NULL, stringsAsFactors = FALSE)

        tg_bg_mf_c <- dplyr::inner_join(tg_mf_c1, bg_mf_c1)
        tg_bg_mf_c$fold_enrich <- ((tg_bg_mf_c$tg_motif_count + 0)/(dim(tg2)[1] + 0))/((tg_bg_mf_c$bg_motif_count + 0)/(dim(bg2)[1] + 0))

        # computes P-values based on binomial distribution
        compute_pval_binom <- function(data1 = data1, tot_tg = dim(tg2)[1], tot_bg = dim(bg2)[1]){
            return(binom.test(x = data1[1], n = tot_tg, p = data1[2]/tot_bg, alternative = "greater")$p.value)
        }

        # computes P-values based on hypergeometric distribution or fisher's exact test
        compute_pval_hyper <- function(data1 = data1, tot_tg = dim(tg2)[1], tot_bg = dim(bg2)[1]){
            return(phyper(q = data1[1]-1, m = (data1[1] + data1[2]), n = (tot_tg + tot_bg) - (data1[1] + data1[2]), k = tot_tg, lower.tail = FALSE, log.p = FALSE))
        }

        if(Pvalue_computation == "binom"){
            tg_bg_mf_c$pval <- apply(X = tg_bg_mf_c[,2:3], MARGIN = 1, FUN = compute_pval_binom, tot_tg = dim(tg2)[1], tot_bg = dim(bg2)[1])
        } else (tg_bg_mf_c$pval <- apply(X = tg_bg_mf_c[,2:3], MARGIN = 1, FUN = compute_pval_hyper, tot_tg = dim(tg2)[1], tot_bg = dim(bg2)[1]))
        return(tg_bg_mf_c)
    }

    motif_enrich <- motif_pval(tg2 = f1, bg2 = f2, Pval_dist = Pvalue_computation)
    motif_enrich$pval_adj <- p.adjust(motif_enrich$pval, method = Pvalue_adjust_method)
    motif_enrich <- dplyr::inner_join(motif_enrich, motif_anno[, 1:4])
    motif_enrich <- motif_enrich %>% dplyr::filter(.data$pval_adj < Pvalue_threshold) %>% dplyr::arrange(pval_adj)
    message("Check point 8 passed")

    # get rid off the inf values if any
    motif_res <- motif_enrich[, c(1,7,8,9,4,5,6)] %>% dplyr::filter_if(~is.numeric(.), dplyr::all_vars(!is.infinite(.)))

    return(motif_res)

}




