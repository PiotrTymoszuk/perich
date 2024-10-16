# perich
Permutation Enrichment Testing for Binary Data

## A bigger picture

Many statistical tools like Fisher's exact and $\chi^2$ test or logistic regression used for comparison of binary information between some analysis groups assume independent and identical distributions (IID). This assumption may be severely violated in many biological data sets such as genetic alterations in cancer samples. In this particular example, it not hard to imagine that mutations, amplifications and deletions may arise at overall rated that differ deeply between cancers. In particular, in cancers with genetically unstable phenotype, e.g. due to inactivation of DNA repair mechanisms or cell cycle checkpoints, many genetic alterations may simply arise by chance. 
Tools of `perich` package try to overcome the IID assumption by weighted permutation testing, i.e. comparision of the actual rates of a binary event (e.g. presence of a mutation) with N radom reshuffles of the binary event vector. By default, the permutation weights or probability of an event for an observation equal to the overall event rate for all available binary variables. In our practical example of genetic alterations, the weights correspond to the fraction of mutated genes in the cancer samples.

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

