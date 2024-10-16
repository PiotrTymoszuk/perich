# Permutation testing functions and methods

#' Permutation testing for differences in frequency of a binary event.
#'
#' @description
#' Permutation testing for an enrichment or depletion of a binary event
#' (e.g. presence of a somatic mutation) in analysis groups.
#'
#' @details
#'
#' __General principle__
#'
#' By principle, the permutation test compares the observed rate of a binary
#' event in an analysis group with N random permutations of the event vector
#' weighted by overall event rate in the observations.
#' In an illustrative example, the test may be used to assess enrichment or
#' depletion of a given somatic mutation in analysis group as compared with
#' random permutations which are weighted by the general mutation rates in
#' biological samples. In this particular context, the enrichment estimates
#' (Enrichment Scores, which resemble odds ratios) and p values returned
#' by the test are expected to be more reliable as compared with classical
#' enrichment assessment methods such as Fisher's exact test or logistic
#' regression, which assume that the mutation rates are independently and
#' identically distributed in the analysis groups.
#'
#' __Enrichment or depletion of events in an analysis group__
#'
#' The \eqn{H_0} hypothesis of the test is that the event frequency in an
#' analysis group is attributed solely to the random variability and differences
#' in overall event rate between single samples.
#'
#' Lets code the vector of 0/1-coded events for \eqn{n} observations as
#' \eqn{X = [x_1, x_2, ..., x_n]}.
#' The observations are assigned to \eqn{p} analysis groups defined by the
#' vector \eqn{F = [f_1, f_2, ..., f_n]}.
#' For the observations, weights \eqn{WT = [wt_1, wt_2, ..., wt_n]} are computed,
#' where each weight corresponds to the column mean of the background matrix,
#' i.e. the overall probability of the event.
#' Next, \eqn{N_{iter}} permutations \eqn{P = [p_1, p_2, ..., p_n]} of \eqn{X},
#' i.e. random reshuffles, are drawn where the event probabilities for the
#' observations are weighted by \eqn{WT} vector.
#' Those permutations are used to estimate the expected frequency of the event
#' in each of the analysis group.
#'
#' To test for enrichment in an analysis group, we compare the observed sum of
#' events in the \eqn{f^{th}} group:
#'
#' \deqn{O_f = \sum X_f}
#'
#' with the expected sum of events in the \eqn{f^{th}} group for iterations
#' \eqn{i, ..., N_{iter}}:
#'
#'  \deqn{E_{f, i} = \sum P_{f, i}}
#'
#' If \eqn{O_f > E_{f, i}} in a particular iteration, this provides a case
#' for \eqn{H_0} rejection, which we code \eqn{H_{f, i} = 1}. Otherwise, we
#' code \eqn{H_{f, i} = 0}. The p-value for significant enrichment in the
#' \eqn{f^{th}} group is simply computed by summing the \eqn{H_0} rejection
#' cases over iterations and dividing it by the total iteration number:
#'
#'  \deqn{p_f = \frac{(\sum_{i}^{N_{iter}} H_{f, i}) + 1}{N_{iter}}}
#'
#' The enrichment score in the \eqn{f^{th}} group is calculated as a ratio of
#' the observed and expected counts of events averaged over the iterations:
#'
#' \deqn{ES_f = \frac{\sum_{i}^{N_{iter}} \frac{O_f + l}{E_{f, i} + l}}{N_[iter]}}
#'
#' Note that in the formula above, the expected and observed counts are
#' corrected by addition of a constant \eqn{l} Laplace smoother, which is there
#' to prevent division by 0 and expansion of \eqn{ES_f} for small values of
#' \eqn{E_{f, i}}. Confidence intervals of \eqn{ES_f} are derived either
#' from percentiles or computed with a more robust
#' [BCA algorithm](https://blogs.sas.com/content/iml/2017/07/12/bootstrap-bca-interval.html).
#'
#' __Global difference between the observed and expected event rates__
#'
#' To assess the global difference between the observed rate of the event and
#' the expected event rate attributable solely to varying overall event rates in
#' the observations, a modified Pearson's \eqn{\chi^2} test is used.
#' The Yates-corrected test statistic \eqn{\chi^2} is computed as follows:
#'
#' \deqn{\chi^2 = \sum_{f = 1}^{k} \frac{(\mid O_f - E_f \mid - 0.5)^2}{E_f}}
#'
#' for \eqn{f = 1, ..., k} analysis groups, \eqn{O_f} observed frequency of the
#' event in the f^th^ group, and \eqn{E_f} expected frequency of the event in
#' the f^th^ group obtained by weighted permutation.
#' Effect size of the global difference is measured by Cramer's V:
#'
#' \deqn{V = \sqrt{\frac{\chi^2/N}{min(k - 1, r - 1)}}},
#'
#' where \eqn{N} denotes the total number of observation, \eqn{k} stands for
#' the analysis group number, and \eqn{r} is the minimum size of the analysis
#' groups.
#'
#' __Background object__
#'
#' This object provided as `background` argument stores event rates for all
#' measured variables (e.g. all somatic mutations returned in an sequencing
#' experiment) and will be used to estimate observation-specific weights
#' applied to the algorithm permutations.
#' It may be provided in three forms:
#'
#' * as a numeric vector, whose length corresponds to the number of observations.
#' This option may be used when the user intends to used pre-computed
#' observation weights such as total mutation burden. Zero and negative values
#' are not allowed.
#'
#' * as a binary numeric matrix, i.e. matrix that contains 0 and 1 values only.
#' Features are represented by the columns, observations are represented by
#' the rows
#'
#' * as a binary data frame or a data frame that stores two level factors.
#' Analogically, features are represented by the columns and observations are
#' coded by the rows
#'
#' If no background is
#' provided, the X matrix or a data frame will be used instead; for just few
#' variables of interest, this is certainly not an optimal solution
#'
#' @return a data frame or a matrix, or a named list of data frames or matrices
#' with the following columns:
#'
#' * __variable__: names of the variables
#'
#' * __strata__: the analysis group, named after levels of factor `f`, that
#' defines the analysis groups
#'
#' * __ES__: enrichment score
#'
#' * __lower_ci__ and __upper_ci__: lower and upper confidence interval of
#' the enrichment score
#'
#' * __p_value__: p value for enrichment or depletion in the analysis groups
#'
#' * __n_observed__: observed numbers of events in the analysis groups
#'
#' * __n_expected__: expected numbers of events in the analysis groups averaged
#' by all iterations
#'
#' * __n_total__: total numbers of complete observations in the analysis groups
#'
#' * __chisq__: value of \eqn{\chi^2} statistic for global differences between
#' the observed and expected event rates
#'
#' * __df__: degrees of freedom of the \eqn{\chi^2} statistic. Essentially
#' the number of analysis group minus one
#'
#' * __cramer_v__: effect size of the global difference in event rates
#'
#' * __global_p_value__: p value of the global difference in event rates
#'
#' If `as_data_frame = FALSE`, a matrix or a list of matrices is returned, with
#' strata names in row names, as well as p and global p values corrected for
#' multiple testing.
#'
#' @param x a binary numeric vector, a two-level factor, a binary numeric
#' matrix, or a data frame of numeric vectors or two-level factors. If a matrix
#' or a data frame is provided, the columns will define the variables and
#' observations are represented by rows.
#' @param f a factor that defines the analysis groups
#' @param background a numeric vector of length equal to `f`'s length, a binary
#' numeric matrix, a data frame of numeric vectors or two-level factors.
#' See Details.
#' @param alternative a string that specifies if the test should check for
#' enrichment ('greater'), depletion ('less') or both ('both', the default
#' option).
#' @param n_iter number of algorithm iterations.
#' @param laplace Laplace smoother used for calculation of the enrichment score.
#' @param ci_type type of the confidence intervals: bias-corrected and
#' accelerated ('bca', the default option), or percentile confidence intervals
#' ('percentile').
#' @param ci_level alpha level for calculation of the confidence intervals of
#' the enrichment scores.
#' @param as_data_frame logical, should the output of the function be coerced
#' to a data frame or a list of data frames? If `TRUE`, this may make the
#' function  slower and memory-consuming.
#' @param compress logical, should the output list be compressed to a single
#' data frame? If `TRUE`, this may make the function  slower and
#' memory-consuming. Ignored if `as_data_frame = FALSE`.
#' @param adj_method method of p value adjustment as defined by
#' \code{\link[stats]{p.adjust}}. By default no multiple testing adjustment is
#' applied to the result.
#' Ignored if `compress = FALSE` or `as_data_frame = FALSE`.
#' @param .parallel logical, should the function be run in parallel? This option
#' works only if a suitable parallel backend is declared via
#' \code{\link[future]{plan}}.
#' @param .n_chunks number of chunks for the parallel run. If not provided,
#' the number will be inferred from the number of available cores.
#' @param .paropts options for the parallel run provided as a
#' \code{\link[furrr]{furrr_options}} object.
#' @param ... additional arguments passed to the methods.
#'
#' @export

  perTest <- function(x, ...) UseMethod('perTest')

#' @rdname perTest
#' @export

  perTest.default <- function(x,
                              f,
                              background,
                              alternative = c('both', 'greater', 'less'),
                              n_iter = 1000,
                              laplace = 1,
                              ci_type = c('bca', 'percentile'),
                              ci_level = 0.95,
                              as_data_frame = FALSE, ...) {

    ## input control for the x vector --------

    stopifnot(is.atomic(x))

    err_txt <- "'x' has to be a binary numeric vector or a two-level factor"

    if(!is.numeric(x)) {

      if(!is.factor(x)) {

        stop(err_txt, call. = FALSE)

      }

    }

    if(is.numeric(x)) {

      x <- as.integer(x)

    } else {

      x <- as.integer(x) - 1

    }

    if(any(na.omit(x) < 0)) stop(err_txt, call. = FALSE)
    if(any(na.omit(x) > 1)) stop(err_txt, call. = FALSE)

    ## input control for the f factor --------

    if(!is.factor(f)) stop("'f' has to be a factor", call. = FALSE)

    f_levels <- levels(f);

    if(length(f_levels) == 1) {

      stop("'f' has to have at least two levels", call. = FALSE)

    }

    if(length(f) != length(x)) {

      stop("Lengths of 'f' and 'x' must be equal", call. = FALSE)

    }

    ## input control for the background object ---------

    err_txt <- paste("'background' has to be a numeric vector",
                     "a binary numeric matrix, or",
                     "a binary numeric data frame or a data frame of two-level",
                     "factors")

    if(!is.matrix(background) & !is.data.frame(background)) {

      if(!is.numeric(background)) {

        stop(err_txt, call. = FALSE)

      }

    }

    if(is.matrix(background) | is.data.frame(background)) {

      if(nrow(background) != length(f)) {

        stop("The row number in 'background' must be equal to the length of 'f'",
             call. = FALSE)

      }

      if(is.data.frame(background)) background <- process_df(background)

      if(!is.numeric(background)) stop(err_txt, call. = FALSE)

      if(any(!na.omit(background) %in% c(0, 1))) stop(err_txt, call. = FALSE)

    } else {

      if(length(background) != length(f)) {

        stop("The length of 'background' must be equal to the length of 'f'",
             call. = FALSE)

      }

    }

    ## input control for the remaining arguments --------

    alternative <- match.arg(alternative[1],
                             c('both', 'greater', 'less'))

    n_iter = as.integer(n_iter[1])

    stopifnot(n_iter > 0)

    stopifnot(is.numeric(laplace))

    laplace <- laplace[1]

    if(laplace < 0) stop("'laplace' must be > 0", call. = FALSE)

    ci_type <- match.arg(ci_type[1], c('bca', 'percentile'))

    stopifnot(is.numeric(ci_level))

    ci_levels <- ci_level[1]

    if(ci_level < 0 | ci_level > 1) {

      stop("'ci_level' must be within [0, 1] range", call. = FALSE)

    }

    ## permutation testing ----------

    if(is.matrix(background)) {

      test_res <- testVector(x,
                             f = f,
                             background = background,
                             weights = NULL,
                             laplace = laplace,
                             n_iter = n_iter,
                             ci_type = ci_type,
                             alternative = ci_type)

    } else {

      test_res <- testVector(x,
                             f = f,
                             background = NULL,
                             weights = background,
                             laplace = laplace,
                             n_iter = n_iter,
                             ci_type = ci_type,
                             alternative = ci_type)

    }

    ## formatting the output -------

    if(!as_data_frame) {

      rownames(test_res) <- f_levels

      return(test_res)

    }

    col_names <- colnames(test_res)

    test_res <- cbind(c(strata = f_levels),
                      as.data.frame(test_res))

    colnames(test_res) <- c('strata', col_names)

    test_res

  }

#' @rdname perTest
#' @export

  perTest.matrix <- function(x,
                             f,
                             background = NULL,
                             alternative = c('both', 'greater', 'less'),
                             n_iter = 1000,
                             laplace = 1,
                             ci_type = c('bca', 'percentile'),
                             ci_level = 0.95,
                             as_data_frame = FALSE,
                             compress = FALSE,
                             adj_method = 'none',
                             .parallel = TRUE,
                             .n_chunks = NULL,
                             .paropts = furrr_options(seed = TRUE), ...) {

    ## entry control for the x matrix -------

    err_txt <- paste("'x' has to be a binary numeric matrix, a binary numeric",
                     "data frame, or a data frame of two-level vectors.")

    if(!is.matrix(x)) stop(err_txt, call. = FALSE)

    if(!is.numeric(x)) stop(err_txt, call. = FALSE)

    if(any(!na.omit(x) %in% c(0, 1))) stop(err_txt, call. = FALSE)

    if(is.null(colnames(x))) {

      colnames(x) <- paste0('feature_', 1:ncol(x))

    }

    ## entry control for the f factor --------

    if(!is.factor(f)) stop("'f' has to be a factor", call. = FALSE)

    f_levels <- levels(f);

    if(length(f_levels) == 1) {

      stop("'f' has to have at least two levels", call. = FALSE)

    }

    if(length(f) != nrow(x)) {

      stop("Lengths of 'f' and number of rows in 'x' must be equal",
           call. = FALSE)

    }

    ## input control for the background object ---------

    err_txt <- paste("'background' has to be a numeric vector",
                     "a binary numeric matrix, or",
                     "a binary numeric data frame or a data frame of two-level",
                     "factors")

    if(!is.matrix(background) & !is.data.frame(background)) {

      if(!is.numeric(background)) {

        stop(err_txt, call. = FALSE)

      }

    }

    if(is.matrix(background) | is.data.frame(background)) {

      if(nrow(background) != length(f)) {

        stop("The row number in 'background' must be equal to the length of 'f'",
             call. = FALSE)

      }

      if(is.data.frame(background)) background <- process_df(background)

      if(!is.numeric(background)) stop(err_txt, call. = FALSE)

      if(any(!na.omit(background) %in% c(0, 1))) stop(err_txt, call. = FALSE)

    } else {

      if(length(background) != length(f)) {

        stop("The length of 'background' must be equal to the length of 'f'",
             call. = FALSE)

      }

    }

    ## input control for the remaining arguments --------

    alternative <- match.arg(alternative[1],
                             c('both', 'greater', 'less'))

    n_iter = as.integer(n_iter[1])

    stopifnot(n_iter > 0)

    stopifnot(is.numeric(laplace))

    laplace <- laplace[1]

    if(laplace < 0) stop("'laplace' must be > 0", call. = FALSE)

    ci_type <- match.arg(ci_type[1], c('bca', 'percentile'))

    stopifnot(is.numeric(ci_level))

    ci_levels <- ci_level[1]

    if(ci_level < 0 | ci_level > 1) {

      stop("'ci_level' must be within [0, 1] range", call. = FALSE)

    }

    ## testing ---------

    if(!.parallel | inherits(plan(), 'sequential')) {

      if(is.matrix(background)) {

        test_res <- testMatrix(x,
                               f = f,
                               background = background,
                               weights = NULL,
                               laplace = laplace,
                               n_iter = n_iter,
                               ci_type = ci_type,
                               alternative = alternative)

      } else {

        test_res <- testMatrix(x,
                               f = f,
                               background = NULL,
                               weights = background,
                               laplace = laplace,
                               n_iter = n_iter,
                               ci_type = ci_type,
                               alternative = alternative)

      }

    } else {

      if(is.null(.n_chunks)) .n_chunks = availableCores() - 1

      index_lst <- split_vector(1:ncol(x), .n_chunks)

      if(is.matrix(background)) {

        test_res <-
          future_map(index_lst,
                     ~testMatrix(x[, .x, drop = FALSE],
                                 f = f,
                                 background = background,
                                 weights = NULL,
                                 laplace = laplace,
                                 n_iter = n_iter,
                                 ci_type = ci_type,
                                 alternative = alternative),
                     .options = .paropts)

      } else {

        test_res <-
          future_map(index_lst,
                     ~testMatrix(x[, .x, drop = FALSE],
                                 f = f,
                                 background = NULL,
                                 weights = weights,
                                 laplace = laplace,
                                 n_iter = n_iter,
                                 ci_type = ci_type,
                                 alternative = alternative),
                     .options = .paropts)

      }

      test_res <- do.call('c', test_res)

      names(test_res) <- colnames(x)

    }

    ## formatting the output: a list of matrices ---------

    if(!as_data_frame) {

      for(i in seq_along(test_res)) {

        rownames(test_res[[i]]) <- f_levels

      }

      return(test_res)

    }

    ## formatting the output: a list of data frames or a single data frame -------

    col_names <- c('variable', 'strata', colnames(test_res[[1]]))

    test_res <- map2(test_res, colnames(x),
                     ~cbind(cbind(variable = rep(.y, nrow(.x)),
                                  strata = f_levels),
                            as.data.frame(.x)))

    test_res <- map(test_res,
                    set_colnames_,
                    col_names)

    if(!compress) return(test_res)

    test_res <- do.call('rbind', test_res)

    if(adj_method == 'none') {

      test_res[['p_adjusted']] <- test_res[['p_value']]
      test_res[['global_p_adjusted']] <- test_res[['global_p_value']]

      return(test_res)

    }

    ## adjustment of the category p values of differences

    test_res[['p_adjusted']] <-
      p.adjust(test_res[['p_value']], adj_method)

    ## adjustment for the global p values: the adjustment is required only
    ## for the single p values per variable

    global_p_tbl <- map_dfr(split(test_res[c('variable', 'global_p_value')],
                                  test_res[['variable']]),
                            ~.x[1, , drop = FALSE])

    global_p_tbl[['global_p_adjusted']] <-
      p.adjust(global_p_tbl[['global_p_value']], adj_method)

    test_res <- merge(test_res,
                      global_p_tbl[c('variable', 'global_p_adjusted')],
                      by = 'variable')

    test_res

  }

#' @rdname perTest
#' @export

  perTest.data.frame <- function(x,
                                 f,
                                 background = NULL,
                                 alternative = c('both', 'greater', 'less'),
                                 n_iter = 1000,
                                 laplace = 1,
                                 ci_type = c('bca', 'percentile'),
                                 ci_level = 0.95,
                                 as_data_frame = FALSE,
                                 compress = FALSE,
                                 adj_method = 'none',
                                 .parallel = TRUE,
                                 .n_chunks = NULL,
                                 .paropts = furrr_options(seed = TRUE), ...) {

    ## the input control and calculations are done by the matrix method

    x <- process_df(x)

    perTest.matrix(x,
                   f,
                   background = background,
                   alternative = alternative,
                   n_iter = n_iter,
                   laplace = laplace,
                   ci_type = ci_type,
                   ci_level = ci_level,
                   as_data_frame = as_data_frame,
                   compress = compress,
                   adj_method = adj_method,
                   .parallel = .parallel,
                   .n_chunks = .n_chunks,
                   .paropts = .paropts)

  }

# END ---------
