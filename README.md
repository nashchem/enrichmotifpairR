README
================

## Introduction

The identification and differential enrichment of transcription factor
(TF) motifs in a given set of genomic regions relative to control
regions is a common task in regulatory genomics. Here, we present
enrichmotifpairR, an R package for identification of differentially
enriched TF motifs and their binding partner motifs, or enriched motif
pairs, in a given set of genomic regions relative to a control set of
genomic regions.

![enrichmotifpairR framework](supple_Figure1.png)

## Installation

First install the genome annotation packages.

``` r
# install devtools and biocmanager if necessary
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages(c("BiocManager", "devtools"))

# Install genome annotations either hg19 or hg38 depending on your input genome 

## Install hg19 
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")

## Install hg38
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")


```

Next install the `enrichmotifpairR` package

``` r
# installation
devtools::install_github("nashchem/enrichmotifpairR")

# load
library(enrichmotifpairR)
```

## Use case examples

See the `vignettes` directory for several use case examples.


## Shiny app
Web-based interactive shiny app is also available at [shiny app](https://hawkinslab.shinyapps.io/EnrichMotifPair/).


