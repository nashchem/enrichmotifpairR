#' Finds and removes duplicate motif pairs from the output of
#' "findEnrichMotifPair" or "findEnrichMotifPairAll"
#'
#' @param data A data frame, the output of "findEnrichMotifPair" or "findEnrichMotifPairAll"
#' @return A data frame same as input containing unique (non-redundant) motif pairs with
#' P-values, fold-enrichment and annotations
#' @examples
#' \dontrun{
#' removeDuplicateTFPairs(data)
#' }

#' @import stats
#' @importFrom rlang .data
#' @importFrom dplyr inner_join
#' @importFrom dplyr %>%
#' @importFrom assertthat assert_that
#'

removeDuplicateTFPairs <- function(data = data){
    data <- data %>% dplyr::arrange(pval_adj) %>%
        dplyr::mutate(key = paste0(pmin(TF_name_1, TF_name_2), 
                                   pmax(TF_name_1, TF_name_2), sep = "")) %>% 
        dplyr::distinct(key, .keep_all = TRUE) %>% 
        dplyr::filter(TF_name_1 != TF_name_2) %>% 
        dplyr::select(-key)
    return(data)
}
