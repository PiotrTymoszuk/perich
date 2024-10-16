# Development of the Gibbs samples for simulation of the expected event
# probabilities given that the variables are dependent

# tools and data --------

  ## R packages

  library(tidyverse)
  library(perich)
  library(bindata)

  library(future)

  library(microViz)

  ## data

  data("tcga_mutations")

  tcga_mutations <- tcga_mutations %>%
    column_to_rownames('sample_id')

# preparing objects for testing -------

  ## a numeric matrix with binary indices for all mutations

  background <- tcga_mutations %>%
    select(-clust_id) %>%
    as.matrix

  ## a group assignment vector

  g <- tcga_mutations$clust_id

  ## features of interest: genes found to be mutated in at least 10 samples

  mutation_counts <- background %>%
    colSums %>%
    sort(decreasing = TRUE)

  genes <- names(mutation_counts)[mutation_counts >= 10]

# development --------

  ## matrix of conditional probabilities

  background[, genes[1:100], drop = FALSE] %>%
    condprob()
