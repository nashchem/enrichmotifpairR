#' Finds GC and length matched random genomic sequences
#'
#' @param input_data A data frame containing peaks coordinates with chr, start,
#' and end as first three columns
#' @param genome_ver Genome version, it should be either "hg19" or "hg38"
#' @param xfold An integer, number of times the backgrouns seqs should be
#' generated
#' @param batchsize An integer, number of random sequences tested in each trial
#' @param nMaxTrials An integer, number of maximum trials
#' @return A GRanges object with coordinates
#' @examples
#' \dontrun{
#' generateBackgroundSeqs(input_data, genome_ver = "hg38",
#' xfold = 1, batchsize = 10000, nMaxTrials = 50)
#' }
#' @export

generateBackgroundSeqs <- function (input_data, genome_ver = "hg38", xfold = 1,
                                    batchsize = 30000, nMaxTrials = 50) {
    # all borrowed from gkmSVM package and modified
    # checking input data and parameters
    assertthat::assert_that(is.data.frame(input_data))
    assertthat::assert_that(is.character(genome_ver))
    assertthat::assert_that(assertthat::is.number(xfold))
    assertthat::assert_that(assertthat::is.number(batchsize))
    assertthat::assert_that(assertthat::is.number(nMaxTrials))
    GC_match_tol = 0.03
    length_match_tol = 0.03
    if(genome_ver == "hg19"){
        genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
    } else {genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38}
    seqnams = GenomeInfoDb::seqnames(genome)
    chrlens = GenomeInfoDb::seqlengths(genome)
    chrpos = cumsum(as.numeric(chrlens))
    pmax = max(chrpos)
    chrpos = c(chrpos, 1e+12)
    chrpos0 = c(0, chrpos)
    ichrA = as.character(names(chrlens))
    getichrpos = function(ipos) {
        j = order(ipos)
        ipos = sort(ipos)
        ci = 1
        res = rep(NA, length(ipos))
        for (i in 1:length(ipos)) {
            while (ipos[i] > chrpos[ci]) {
                ci = ci + 1
            }
            res[j[i]] = ci
        }
        return(res)
    }
    generateRandomGenSeqs = function(seqlens) {
        rpos = sample(pmax, length(seqlens), replace = TRUE)
        ichr1 = getichrpos(rpos)
        ichr2 = getichrpos(rpos + seqlens)
        jj = which(ichr1 != ichr2)
        while (length(jj) > 0) {
            rpos[jj] = sample(pmax, length(jj), replace = TRUE)
            ichr1 = getichrpos(rpos)
            ichr2 = getichrpos(rpos + seqlens)
            jj = which(ichr1 != ichr2)
        }
        chr = ichrA[ichr1]
        start = rpos - chrpos0[ichr1]
        names <- chr
        ranges <- IRanges::IRanges(start = start, width = seqlens)
        strand <- BiocGenerics::strand(sample(c("+", "-"),
                                              length(names), replace = TRUE))
        gr <- GenomicRanges::GRanges(seqnames = names, ranges = ranges,
                                     strand = strand)
    }
    inbed <- input_data
    inBed <- inbed %>%
        GenomicRanges::makeGRangesFromDataFrame(ignore.strand=TRUE,
                                                keep.extra.columns=TRUE)
    jj = which(is.na(match(as.character(inbed$seqnames),
                           as.character(seqnams))))
    if (length(jj) > 0) {
        cat(paste("ERROR: Chromosome name not recognized for",
                  length(jj), "sequences.\n"))
        cat(unique(as.character(inbed$seqnames[jj])))
        return(NULL)
    }
    jj = which(inbed$end >
                   GenomeInfoDb::seqlengths(genome)[as.character(inbed$seqnames)])
    if (length(jj) > 0) {
        cat("ERROR: Region outside chromosome. (Check the genome version) \n")
        print(inbed[jj, ])
        return(NULL)
    }
    gcContent <- function(seqs) {
        alf <- Biostrings::alphabetFrequency(seqs, as.prob = TRUE)
        gc = rowSums(alf[, c("G", "C"), drop = FALSE])
    }
    matchSeqs = function(gc1, gc2, len1, len2,
                         gc_th = 0.02, len_th = 0.02) {
        len_th = len_th * len1
        i1 = order(gc1)
        i2 = order(gc2)
        gc1 = gc1[i1]
        gc2 = gc2[i2]
        len1 = len1[i1]
        len2 = len2[i2]
        #rpt1 = rpt1[i1]
        #rpt2 = rpt2[i2]
        gc2 = c(gc2, 1e+10)
        len_th = len_th[i1]
        m2 = 1
        N = length(i1)
        N2 = length(i2)
        mtc1 = rep(NA, N)
        mtc2 = rep(0, length(i2))
        for (i in 1:N) {
            gc1i = gc1[i]
            len1i = len1[i]
            #rpt1i = rpt1[i]
            len_thi = len_th[i]
            while (gc1i - gc2[m2] > gc_th) {
                m2 = m2 + 1
            }
            if (m2 <= N2) {
                m2b = m2
                while (gc2[m2b] - gc1i <= gc_th) {
                    if ((mtc2[m2b] == 0) & (abs(len1i - len2[m2b]) <=
                                            len_thi)) {
                        mtc2[m2b] = i
                        mtc1[i] = m2b
                        if (m2b == m2) {
                            m2 = m2 + 1
                        }
                        break
                    }
                    m2b = m2b + 1
                }
            }
            else {
                break
            }
        }
        mtc1 = i2[mtc1]
        res = rep(NA, N)
        res[i1] = mtc1
        return(res)
    }

    inbed = GenomicRanges::as.data.frame(inBed)
    cat("importing sequences", "from", GenomeInfoDb::bsgenomeName(genome),
        "\n")
    inSeqs = Biostrings::getSeq(genome, inBed)
    seqlens = inbed$width
    cat(" calculating GC distributions\n")
    inGC = gcContent(inSeqs)
    nout = round(nrow(inbed) * xfold)
    outbed = matrix(ncol = ncol(inbed), nrow = nout)
    outSeq = rep(inSeqs, length = nout)
    colnames(outbed) = colnames(inbed)
    unmatched = 1:length(outSeq)
    desGC = rep(inGC, length = nout)
    desLens = rep(seqlens, length = nout)
    for (iter in 1:nMaxTrials) {
        if (length(unmatched) > 0) {
            cat(" Trial", iter, "out of", nMaxTrials, "\n")
            rndBed = generateRandomGenSeqs(rep(desLens[unmatched],
                                               length.out = batchsize))
            rndbed = GenomicRanges::as.data.frame(rndBed)
            cat(" importing sequences\n")
            rndSeqs = Biostrings::getSeq(genome, rndBed)
            rndGC = gcContent(rndSeqs)
            cat(" matching sequences\n")
            mtc = matchSeqs(desGC[unmatched], rndGC, desLens[unmatched],
                            BiocGenerics::width(rndBed), gc_th = GC_match_tol,
                            len_th = length_match_tol
            )
            jj = which(!is.na(mtc))
            if (length(jj) > 0) {
                outbed[unmatched[jj], 1:5] = as.matrix(rndbed[mtc[jj],
                                                              ])
                outSeq[unmatched[jj], ] = rndSeqs[mtc[jj],
                                                  ]
                unmatched = unmatched[-jj]
            }
            cat(nrow(outbed) - length(unmatched), "sequences found so far, ",
                length(unmatched), " remaining.\n")
        }
    }
    if (length(unmatched) > 0) {
        outbed = outbed[-unmatched, ]
        outSeq = outSeq[-unmatched, ]
    }
    outbed = gsub(" ", "", outbed)
    out_file <- as.data.frame(outbed[, 1:3], stringsAsFactors = FALSE) %>%
        GenomicRanges::makeGRangesFromDataFrame(ignore.strand=TRUE,
                                                keep.extra.columns=TRUE)
    return(out_file)
}

