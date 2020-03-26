#' Generates scrambled sequences from genomic regions
#'
#' @param input_data A data frame containing peaks coordinates with chr, start,
#' and end as first three columns
#' @param genome_ver Genome version, it should be either "hg19" or "hg38"
#' @return A DNAStringSet instance with scrambled seqs
#' @examples
#' \dontrun{
#' generateScrambleSeqs(input_data, genome_ver = "hg38")
#' }
#' @export
#' @importFrom S4Vectors endoapply

generateScrambleSeqs <- function(input_data, genome_ver = "hg38"){
    if(genome_ver == "hg19"){
        genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
    } else {genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38}
    input_gr <- input_data %>% GenomicRanges::makeGRangesFromDataFrame(ignore.strand=TRUE, keep.extra.columns=TRUE)
    input_seqs = Biostrings::getSeq(genome, input_gr)
    scr_seqs <- endoapply(input_seqs, sample)
    return(scr_seqs)
}
