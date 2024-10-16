# Description of the example data sets.

# TCGA mutation data ------

#' Binary somatic mutation indices in TCGA BLCA urothelial cancers.
#'
#' @description
#' A data frame with binary indexes of mutations in over 17K genes in a set of
#' urothelial cancers assigned to three molecular clusters.
#'
#' @format
#' A data frame with 405 rows (cancer samples) and 17005 variables:
#' * __sample_id__ contains unique sample identifiers
#' * __clust_id__ codes for molecular cluster assignment
#' * the remaining columns are named after HUGO gene symbols and store binarized
#' indexes of mutation indices. 0 represents the wild-type variant, 1 represents
#' at least one somatic mutation in the given gene.
#'
#' @source TCGA BLCA mutation data set generated with the text file deposited at
#' [cBioportal](https://www.cbioportal.org/study/summary?id=blca_tcga_pan_can_atlas_2018).
#'
#' @docType data
#'
#' @name tcga_mutations
#'
#' @usage data(tcga_mutations)

  NULL

# IMvigor mutation data --------

#' Binary somatic mutation indices in TCGA BLCA urothelial cancers.
#'
#' @description
#' A data frame with binary indexes of mutations in over 300 genes in a set of
#' urothelial cancers assigned to three molecular clusters.
#'
#' @format
#' A data frame with 180 rows (cancer samples) and 394 variables:
#' * __sample_id__ contains unique sample identifiers
#' * __clust_id__ codes for molecular cluster assignment
#' * the remaining columns are named after HUGO gene symbols and store binarized
#' indexes of mutation indices. 0 represents the wild-type variant, 1 represents
#' at least one somatic mutation in the given gene.
#'
#' @source IMvigor cohort mutation data set authored by Si Yangming and others
#' (DOI: 10.1038/nature25501) and available as CC 3.0 licensed R package
#' [IMvigor210CoreBiologies](https://github.com/SiYangming/IMvigor210CoreBiologies).
#'
#' @docType data
#'
#' @name imvigor_mutations
#'
#' @usage data(imvigor_mutations)

  NULL
