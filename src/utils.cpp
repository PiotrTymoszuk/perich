/*** Helper functions for transformation of vectors and matrices */

#include <Rcpp.h>
#include <Rmath.h>

#include "utils.h"

using namespace Rcpp;

// checking if a matrix has column and row names

LogicalVector checkMtxNames(const NumericMatrix &x) {

  List s = x.attr("dimnames");  // could be nil or list

  LogicalVector res(2, false);

  if(s.length() == 0) return res;

  if(!Rf_isNull(s[0])) res[0] = true;
  if(!Rf_isNull(s[1])) res[1] = true;

  return res;

}

// checking if a numeric vector has names

bool checkVecNames(const NumericVector &x) {

  CharacterVector xNames = x.names();

  return xNames.length() != 0;

}

// splitting of a numeric vector by categories defined by an integer vector

List splitNumVec(NumericVector x, IntegerVector f) {

  // splitting

  IntegerVector f_cat = unique(f).sort();
  int fLen = f.length();
  int catLen = f_cat.length();

  List res(catLen);

  if(catLen == 1) {

    res[0] = x;
    res.names() = f_cat[0];

    return res;

  }

  if(x.length() != fLen) stop("Improper length of f");

  for(int k = 0; k < catLen; ++k) {

    LogicalVector indexes(fLen, false);

    for(int i = 0; i < fLen; ++i) {

      if(f[i] == f_cat[k]) indexes[i] = true;

    }

    // omitting the NA-storing positions from the output

    LogicalVector na_check = !is_na(f);

    indexes = na_check & indexes;

    res[k] = x[indexes];

  }

  res.names() = f_cat;

  return res;

}

// computation of quantiles and confidence intervals

// [[Rcpp::export]]

NumericVector Quantile(NumericVector x, NumericVector probs) {

  // calculation of quantiles
  // many thanks to https://github.com/RcppCore/Rcpp/issues/967

  const size_t n = x.size(), np = probs.size();

  if (n == 0) return x;
  if (np == 0) return probs;

  NumericVector index = (n - 1.0) * probs, y = x.sort(), x_hi(np), qs(np);
  NumericVector lo = floor(index), hi = ceiling(index);

  for (size_t i = 0; i < np; ++i) {

    qs[i] = y[lo[i]];
    x_hi[i] = y[hi[i]];

    if ((index[i] > lo[i]) && (x_hi[i] != qs[i])) {

      double h;
      h = index[i] - lo[i];
      qs[i] = (1.- h) * qs[i] + h * x_hi[i];

    }

  }

  return qs;

}

// [[Rcpp::export]]

NumericVector perCI(NumericVector theta, double conf_level = 0.95) {

  // computes percentile confidence intervals

  NumericVector ci_probs{(1 - conf_level)/2, (1 + conf_level)/2};

  return Quantile(theta, ci_probs);

}

// [[Rcpp::export]]

NumericVector BCA(NumericVector theta, double conf_level = 0.95) {

  // computes BCA confidence intervals based on the R code
  // provided by the coxed package
  // https://rdrr.io/cran/coxed/src/R/bca.R

  double low;
  double high;

  low = (1 - conf_level)/2;
  high = 1 - low;

  int sims = theta.size();

  NumericVector low_theta;

  low_theta = ifelse(theta < mean(theta), 1.0, 0.0);

  double z_inv;

  z_inv = sum(low_theta)/sims;

  double z;

  z = R::qnorm(z_inv, 0, 1, 1, 0);

  NumericVector U;

  U = (sims - 1) * (mean(theta) - theta);

  double top;
  double under;
  double a;

  top = sum(pow(U, 3));

  under = sum(pow(U, 2));

  under = 6 * std::pow(under, 1.5);

  a = top/under;

  double lower_inv;
  double upper_inv;

  double q_low;
  double q_high;

  q_low = R::qnorm(low, 0, 1, 1, 0);

  lower_inv = z + (z + q_low)/(1 - a * (z + q_low));

  lower_inv = R::pnorm(lower_inv, 0, 1, 1, 0);

  q_high = R::qnorm(high, 0, 1, 1, 0);

  upper_inv = z + (z + q_high)/(1 - a * (z + q_high));

  upper_inv = R::pnorm(upper_inv, 0, 1, 1, 0);

  NumericVector probs = {lower_inv, upper_inv};

  return Quantile(theta, probs);

}

// END
