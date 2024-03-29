% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{plotEnrichment}
\alias{plotEnrichment}
\title{Plots the enrichment values of selected TFs.}
\usage{
plotEnrichment(enrich_motifs, tfs)
}
\arguments{
\item{enrich_motifs}{motif_enrich output dataframe from findEnrichMotifPair}

\item{tfs}{TFs for which to plot enrichment}

\item{metric}{Value to plot; either pval or logFC}
}
\value{
a ggplot2 plot object
}
\description{
Plots the enrichment values of selected TFs.
}
\examples{
\dontrun{
motifs_df <- results$motif_enrich
plotEnrichment(enrich_motifs = motifs_df, tfs = c("SOX2", "NANOG", "MYC", "POU5F1"),
metric = "pval")
}
}
