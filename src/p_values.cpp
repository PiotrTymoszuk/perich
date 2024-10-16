/*** Calculation of sample- and feature-wise probabilities of an event (feature == 1) */

#include <Rcpp.h>
#include <Rmath.h>

#include "utils.h"
#include "p_values.h"

using namespace Rcpp;

// feature-wise probabilities of the event: if no event was found,
// a default, very low probability is inserted instead of p = 0.0

// [[Rcpp::export]]

NumericVector pFeature (NumericMatrix x, double min_p = 1e-6) {

  double n_samples = x.nrow();
  int n_features = x.ncol();

  NumericVector colVec(n_samples);
  NumericVector pValues(n_features);

  for(int i = 0; i < n_features; ++i) {

    colVec = x(_, i);

    pValues[i] = sum(na_omit(colVec))/n_samples;

    if(pValues[i] == 0.0) pValues[i] = min_p;

  }

  if(checkMtxNames(x)[1]) pValues.names() = colnames(x);

  return pValues;

}

// sample-wise probabilities of the event

// [[Rcpp::export]]

NumericVector pSample (NumericMatrix x, double min_p = 1e-6) {

  double n_samples = x.nrow();
  int n_features = x.ncol();

  NumericVector rowVec(n_features);
  NumericVector pValues(n_samples);

  for(int i = 0; i < n_samples; ++i) {

    rowVec = x(i, _);

    pValues[i] = sum(na_omit(rowVec))/n_features;

    if(pValues[i] == 0.0) pValues[i] = min_p;

  }

  if(checkMtxNames(x)[0]) pValues.names() = rownames(x);

  return pValues;

}

// END
