# Non-exported utilities

# Vector splitting -----

#' Split a vector into a list of its fragments.
#'
#' @description
#' The function takes a vector and returns a list of its chunks of possibly
#' similar sizes.
#'
#' @return a list of vectors
#'
#' @param x a vector.
#' @param n number of chunks.

  split_vector <- function(x, n) {

    chunk_size <- ceiling(length(x)/n)

    split_list <- split(x, ceiling(seq_along(x)/chunk_size))

    split_list

  }

# Columns names -------

#' Set column names.
#'
#' @description Sets column names in a matrix or a data frame.
#'
#' @param x a data frame or a matrix.
#' @param nm a character vector with column names.
#'
#' @return a matrix or a data frame.

  set_colnames_ <- function(x, nm) {

    stopifnot(length(nm) == ncol(x))

    colnames(x) <- nm

    x

  }

# Pre-processing of a data frame ---------

#' Pre-processing of a data frame for permutation testing.
#'
#' @description
#' Throws errors if any any of the columns is in an unsupported
#' non-binary format. Converts two-level factors to 0/1 integers.
#'
#' @param x a data frame.
#'
#' @return a 0/1 matrix.

  process_df <- function(x) {

    ## checking for non-numeric and non-factor columns

    stopifnot(is.data.frame(x))

    class_check <-
      map_lgl(x, function(x) is.numeric(x) | is.factor(x))

    if(any(!class_check)) {

      stop("'x' stores non-numeric or non-factor data.", call. = FALSE)

    }

    ## conversion of factors to numeric vectors

    factor_check <- map_lgl(x, is.factor)

    if(sum(factor_check) > 0) {

      x[factor_check] <- map_dfc(x[factor_check],
                                 ~as.numeric(.x) - 1)

    }

    x <- as.matrix(x)

    if(any(!na.omit(x) %in% c(0, 1))) {

      stop("'x' stores non-binary numerics or multi-level factors.",
           call. = FALSE)

    }

    x

  }

# END -------
