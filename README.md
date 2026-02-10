# perich
_Permutation Enrichment Testing for Binary Data_

## A bigger picture

Many statistical tools like Fisher's exact and $\chi^2$ test or logistic regression used for 
comparison of binary information between some analysis groups assume independent and identical distributions (IID). 
This assumption may be severely violated in many biological data sets such as genetic alterations in cancer samples. 
In this particular example, it not hard to imagine that mutations, amplifications and deletions may arise at overall rates 
that differ deeply between cancers.
In particular, in cancers with genetically unstable phenotype, e.g. due to inactivation of DNA repair mechanisms or 
cell cycle checkpoints, many genetic alterations may simply arise by chance. 
Tools of `perich` package try to overcome the IID assumption by weighted permutation testing, i.e. comparison of the actual 
rates of a binary event (e.g. presence of a mutation) with N random reshuffles of the binary event vector. 
By default, the permutation weights or probability of an event for an observation equal to the overall event rate 
for all available binary variables. 
In our practical example of genetic alterations, the weights correspond to the fraction of mutated genes in the cancer samples.

## Installation

You may easily fetch the package via `devtools`: 

```r

devtools::install_github('PiotrTymoszuk/perich')

```

## Terms of use

The package is available under a [GPL-3 license](https://github.com/PiotrTymoszuk/perich/blob/main/LICENSE).

## Contact

The package maintainer is [Piotr Tymoszuk](mailto:piotr.s.tymoszuk@gmail.com).

## Acknowledgements

Many thanks to authors, maintainers and contributors of [Rcpp](https://www.rcpp.org/) and [tidyverse evironment](https://www.tidyverse.org/)

## Basic usage

As noted above, weighted permutation testing may be especially useful at analysis of binary genetic data, 
e.g. comparison of mutation frequency between analysis groups or enrichment in of somatic mutations in 
an analysis group that is unlikely attributed to random differences in overall mutation rates between cancer samples. 
`perich` package comes with two publicly available somatic mutation data sets of urothelial cancers, 
the TCGA BLCA and IMvigor cohort. 
Here, we will test for non-accidental enrichment of somatic mutations in [three molecular clusters](https://github.com/PiotrTymoszuk/BLCA-cluster-paper). 
The mutation data are 0/1-coded, where 0 denotes a wild-type (WT) variant of a gene and 1 codes for presence of 
at least one somatic mutation. 
Rows represent cancer samples, genes are stored in columns:

```r
  ## R packages

  library(tidyverse)
  library(perich)

  library(future)

  library(microViz)

  ## data

  data("tcga_mutations")
  data("imvigor_mutations")

  data <- list(tcga = tcga_mutations,
               imvigor = imvigor_mutations) %>%
    map(column_to_rownames, 'sample_id')

```
```r
> data$tcga[1:10, 1:10]
                clust_id A1BG A1CF A2M A2ML1 A3GALT2 A4GALT A4GNT AAAS AACS
TCGA-2F-A9KO-01       #1    0    0   0     0       0      0     0    0    0
TCGA-2F-A9KP-01       #3    0    0   0     0       0      0     0    0    0
TCGA-2F-A9KQ-01       #3    0    0   0     0       0      0     0    0    0
TCGA-2F-A9KR-01       #3    0    0   0     0       0      0     0    0    0
TCGA-2F-A9KT-01       #2    0    0   0     0       0      0     0    0    0
TCGA-2F-A9KW-01       #1    0    0   0     0       0      0     0    0    0
TCGA-4Z-AA7M-01       #3    0    0   0     0       0      0     0    0    0
TCGA-4Z-AA7N-01       #1    0    0   0     0       0      0     0    0    0
TCGA-4Z-AA7O-01       #3    0    0   0     0       0      0     0    0    0
TCGA-4Z-AA7Q-01       #2    0    0   0     0       0      0     0    0    0
```

For the analysis, a background object is required. 
It is recommended to provide a matrix of all available binary variables as a background object. 
In our case, the background objects will be simply 0/1 mutations for all genes measured in cancer cohorts. 
Additionally, a group-defining factor is required, in our case a vector with molecular subset labels. 
We are going to assess enrichment for genes with mutations in at least 2.5% of cancer samples in each cohort.

```r
  ## numeric matrices with binary indices for all mutations

  background <- data %>%
    map(select, -clust_id) %>%
    map(as.matrix)

  ## group assignment vectors

  g <- data %>%
    map(~.x$clust_id)

  ## features of interest: genes found to be mutated in at least 2.5% of samples

  mutation_counts <- background %>%
    map(colMeans) %>%
    map(sort, decreasing = TRUE)

  genes <- mutation_counts %>%
    map(~.x[.x >= 0.025]) %>%
    map(names)

```

To test for enrichment in the molecular clusters, we are calling `perTest()` function. 
In a fist step, we are checking for rates of mutations in the tumor suppressor gene _TP53_. 
The function takes a binary vector with the binary variable of interest, 
the group-defining factor `f`, a `background` object as obligatory arguments.
The user is welcome to adjust numbers of algorithm iterations, type of confidence intervals or 
testing alternative (both by default: both enrichment and depletion in the clusters will be investigated):

```r

  perTest(data$tcga$TP53,
          f = g$tcga,
          background = background$tcga,
          laplace = 1,
          n_iter = 1000,
          ci_type = 'perc',
          alternative = 'both',
          as_data_frame = TRUE)
```
```r
        strata        ES lower_ci  upper_ci p_value n_observed n_expected n_total    chisq df  cramer_v global_p_value
strata1     #1 1.0150430 0.880597 1.1800000   0.501         58     57.462     119 7.854945  2 0.0984757     0.01969339
strata2     #2 1.2939231 1.119403 1.5306122   0.001         74     57.312     117 7.854945  2 0.0984757     0.01969339
strata3     #3 0.8029207 0.718750 0.9078947   0.002         68     85.226     169 7.854945  2 0.0984757     0.01969339

```
The effect size of enrichment in the clusters is measured by an enrichment score `ES`, 
which is an Laplace smoother-corrected ratio of the observed to the expected event counts obtained by weighted permutations. 
The `ES` is returned with 95% confidence intervals and a p value for enrichment for the clusters. 
The `chisq`, `df`, and `global_p_value` of the output refer to a slightly modified Pearson's $\chi^2$ test 
for differences in mutation frequency between the clusters. 
Cramer's V statistic measures the effect size of such global differences in event rates. 

`perTest()` function accepts also a binary matrix as the first argument. 
In our example, if a matrix or a data frame with 0/1-coded mutations is provided, each of the columns is treated as 
a separate variable. By setting `compress = TRUE`, the output will be coerced to a single data frame. 
`adj_method` allows for specifying a multiple testing adjustment method (Benjamini-Hochberg in our case): 

```r
perTest(data$imvigor[, genes$imvigor[1:10]],
          f = g$imvigor,
          background = background$imvigor,
          laplace = 1,
          n_iter = 1000,
          ci_type = 'bca',
          alternative = 'both',
          as_data_frame = TRUE,
          compress = TRUE,
          adj_method = 'BH') %>%
head
```
```r
  variable strata        ES  lower_ci upper_ci p_value n_observed n_expected n_total    chisq df   cramer_v
1   ARID1A     #1 0.5714481 0.4117647 0.875000   0.014          6     11.787      57 5.382812  2 0.12227942
2   ARID1A     #2 1.6894412 1.2500000 3.750000   0.019         14      8.457      44 5.382812  2 0.12227942
3   ARID1A     #3 1.0417782 0.7727273 1.545455   0.544         16     15.756      79 5.382812  2 0.12227942
4    CCND1     #1 0.7711331 0.5454545 1.500000   0.197          5      7.292      57 0.780918  2 0.04657485
5    CCND1     #2 1.4501887 0.8888889 4.000000   0.222          7      5.173      44 0.780918  2 0.04657485
6    CCND1     #3 1.0942245 0.7333333 1.833333   0.502         10      9.535      79 0.780918  2 0.04657485
  global_p_value p_adjusted global_p_adjusted
1     0.06778557 0.08142857         0.2259519
2     0.06778557 0.08142857         0.2259519
3     0.06778557 0.56275862         0.2259519
4     0.67674619 0.39176471         0.8459327
5     0.67674619 0.39176471         0.8459327
6     0.67674619 0.55777778         0.8459327

```
Finally, by calling a parallel backend via `future()`, we make the function run in parallel, 
which enables a 'genome-wide' analysis: 

```r

 enrichment <- list()

  plan('multisession')

  for(i in names(data)) {

    enrichment[[i]] <- perTest(data[[i]][, genes[[i]]],
                               f = g[[i]],
                               background = background[[i]],
                               laplace = 1,
                               n_iter = 1000,
                               ci_type = 'bca',
                               alternative = 'both',
                               as_data_frame = TRUE,
                               compress = TRUE,
                               adj_method = 'BH',
                               .parallel = TRUE)

  }

  plan('sequential')

```
```r
  ## significant enrichment/depletion may be defined by raw p < 0.01.
  ## We're extracting significantly enriched,
  ## and depleted mutations in each of the molecular subsets

  enrichment <- enrichment %>%
    map(mutate,
        log_ES = log(ES),
        status = ifelse(p_value >= 0.01, 'ns',
                        ifelse(log_ES > 0, 'enriched',
                               ifelse(log_ES < 0, 'depleted', 'ns'))),
        status = factor(status, c('enriched', 'depleted', 'ns')),
        strata = factor(strata)) %>%
    map(as_tibble)

  ## significant enrichment for the molecular subsets

  significant_results <- enrichment %>%
    map(filter, status %in% c('enriched', 'depleted')) %>%
    map(~split(.x, .x$strata))

  significant_genes <- significant_results %>%
    map(map, ~.x$variable)

```
```r
> significant_genes$tcga
$`#1`
 [1] "ASAP1"     "CDH16"     "FGFR3"     "GPR98"     "HIST1H2BK" "KIAA0556"  "LRRN1"     "MEGF8"     "PHLDB2"   
[10] "TCF20"     "ZNF791"   

$`#2`
 [1] "AHNAK"    "ATP2B4"   "CDHR2"    "CDKN2A"   "CMYA5"    "FAM179B"  "FAM83G"   "FANCM"    "FGFR3"    "INADL"   
[11] "IPO9"     "NDST4"    "NFE2L2"   "PALB2"    "PAXBP1"   "TNPO3"    "TOR1AIP1" "TP53"     "TRAPPC10" "ZNF462"  

$`#3`
 [1] "ACTN4"   "ASAP1"   "COL16A1" "FGFR3"   "FLG2"    "FREM1"   "MYO16"   "NBPF1"   "NFE2L2"  "PHLDB2"  "PLA2R1" 
[12] "PLEKHH2" "RB1"     "SUPT20H" "SYNE2"   "TIAM1"   "TP53"    "UBE3A"   "ZNF423" 
```

Subsequently, the significantly enriched genes may be visualized in oncoplots e.g. with function `plot_bionco()` from [`microViz` package](https://github.com/PiotrTymoszuk/microViz): 

```r
  ## oncoplots for the genes found to be significantly enriched in the clusters
  ## working with safely(), bacause in the IMvigor cohorts there may be no
  ## significant enrichment in one of the clusters

  cluster_enriched <- significant_results %>%
    map(map, filter, status == 'enriched') %>%
    map(map, ~.x$variable)

  onco_plots <- list()

  for(i in names(cluster_enriched)) {

    onco_plots[[i]] <-
      list(variables = cluster_enriched[[i]],
           plot_title = paste0(toupper(i),
                               ', subset ',
                               names(cluster_enriched[[i]]),
                               ' signature mutations')) %>%
      pmap(safely(plot_bionco),
           data = data[[i]],
           split_fct = 'clust_id',
           x_lab = 'cancer sample',
           y_lab = NULL,
           hide_x_axis_text = TRUE,
           y_text_face = 'italic') %>%
      map(~.x$result)

  }

```
![image](https://github.com/user-attachments/assets/3f98a981-6e37-4d8a-aa32-ef3d9fd6087d)
