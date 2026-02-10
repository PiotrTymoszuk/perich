# Basic usage of the package

# tools and data --------

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

# preparing objects for testing -------

  ## numeric matrices with binary indices for all mutations,
  ## they'll serve as backgrounds for estimating expacted rates of
  ## alterations occuring at random

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

# permutation testing for single vectors of binary events ------

  set.seed(12345)

  ## testing for enrichment with input vectors with NAs:
  ## any NA values are silently removed

  tricky_vec <- data$tcga$FGFR3
  tricky_vec[10] <- NA
  tricky_vec[199] <- NA

  tricky_f <- g$tcga
  tricky_f[11] <- NA
  tricky_f[200] <- NA

  perTest(tricky_vec,
          f = tricky_f,
          background = background$tcga,
          laplace = 1,
          n_iter = 1000,
          ci_type = 'perc',
          alternative = 'both')

  ## differences in mutation rates of the key tumor suppressor TP53
  ## output as a data frame

  perTest(data$tcga$TP53,
          f = g$tcga,
          background = background$tcga,
          laplace = 1,
          n_iter = 1000,
          ci_type = 'perc',
          alternative = 'both',
          as_data_frame = TRUE)

  perTest(data$imvigor$TP53,
          f = g$imvigor,
          background = background$imvigor,
          laplace = 1,
          n_iter = 1000,
          ci_type = 'perc',
          alternative = 'both',
          as_data_frame = TRUE)

  ## background provided as a numeric vector: zeros are not allowed!

  perTest(data$tcga$FGFR3,
          f = g$tcga,
          background = rowMeans(background$tcga) + 1e6,
          laplace = 1,
          n_iter = 1000,
          ci_type = 'perc',
          alternative = 'both',
          as_data_frame = TRUE)

  perTest(data$imvigor$FGFR3,
          f = g$imvigor,
          background = rowMeans(background$imvigor) + 1e6,
          laplace = 1,
          n_iter = 1000,
          ci_type = 'perc',
          alternative = 'both',
          as_data_frame = TRUE)

# permutation testing for a matrix of events --------

  ## serial runs for the top 10 mutated genes,
  ## output as a single data frame

  perTest(data$tcga[, genes$tcga[1:10]],
          f = g$tcga,
          background = background$tcga,
          laplace = 1,
          n_iter = 1000,
          ci_type = 'bca',
          alternative = 'both',
          as_data_frame = TRUE,
          compress = TRUE,
          adj_method = 'BH')

  perTest(data$imvigor[, genes$imvigor[1:10]],
          f = g$imvigor,
          background = background$imvigor,
          laplace = 1,
          n_iter = 1000,
          ci_type = 'bca',
          alternative = 'both',
          as_data_frame = TRUE,
          compress = TRUE,
          adj_method = 'BH')

  ## parallel run for all genes of interest

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

# Significant enrichment or depletion --------

  ## significant enrichment/depletion is defined by raw p < 0.01.
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

# visualization of the significant results -------

  ## volcano plots of enrichment scores and raw p values
  ## for the significant genes

  volcano_plots <- list()

  for(i in names(significant_results)) {

    volcano_plots[[i]] <-
      list(data = significant_results[[i]],
           plot_title = paste0('Subset ',
                               names(significant_results[[i]]),
                               ', ', toupper(i))) %>%
      pmap(plot_volcano,
           regulation_variable = 'log_ES',
           p_variable = 'p_value',
           label_variable = 'variable',
           signif_level = 0.05,
           regulation_level = 0,
           top_regulated = 100,
           label_type = 'text',
           txt_face = 'italic',
           x_lab = 'log Enrichment Score')

  }

  ## classification of variables by cluster and heat map visualization

  hm_variables <- significant_genes %>%
    map(unlist, use.names = FALSE)

  hm_data <-
    map2(data,
         hm_variables,
         ~select(.x, clust_id, all_of(.y)))

  hm_plots <-
    list(data = hm_data,
         variables = hm_variables,
         plot_title = paste(c('TCGA,', 'IMVIGOR,'),
                            'significant mutations')) %>%
    pmap(heat_map,
         split_fct = 'clust_id',
         normalize = FALSE) %>%
    map2(., c(2, 1), ~.x +
          scale_fill_gradient(low = 'gray60',
                              high = 'orangered2',
                              breaks = c(0, 1),
                              labels = c('WT', 'mutated'),
                              name = '') +
          guides(y = guide_axis(n.dodge = .y),
                 fill = guide_legend()) +
          theme(axis.text.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.y = element_text(face = 'italic')))

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

  ## stack plots for the genes significantly enriched in the clusters

  stack_plots <- list()

  for(i in names(cluster_enriched)) {

    stack_plots[[i]] <-
      list(variables = cluster_enriched[[i]],
           plot_title = paste0(toupper(i),
                               ', subset ',
                               names(cluster_enriched[[i]]),
                               ' signature mutations')) %>%
      pmap(safely(plot_bistack),
           data = data[[i]],
           split_fct = 'clust_id',
           scale = 'percent',
           x_lab = '% of samples with mutatations',
           y_lab = NULL) %>%
      map(~.x$result)

  }

# END ---------
