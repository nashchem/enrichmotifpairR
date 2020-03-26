Introduction:
-------------

Add text here.

Installation:
-------------

Install the package from Github.

add installation command here.

Quick start:
------------

Load the package into R session.

     library(enrichmotifpairR)

The enrichmotifpairR package includes an example dataset called
peak\_data.

    data(peak_data)

    set.seed(109)
    data_results <- findEnrichMotifPair(peak_data[1:1000, ], background_data = NULL,
      genome_ver = "hg38", scramble_data = TRUE,
      motif_database = "JASPAR")
