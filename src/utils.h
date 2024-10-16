#ifndef __utils__
#define __utils__

#include <Rcpp.h>

using namespace Rcpp;

// checks if a matrix has row and columns names

LogicalVector checkMtxNames(const NumericMatrix &x);

// checks if a vector has names

bool checkVecNames(const NumericVector &x);

// splitting of a numeric vector by unique values of an integer vector

List splitNumVec(NumericVector x, IntegerVector f);

// quantiles and confidence intervals

NumericVector Quantile(NumericVector x, NumericVector probs);
NumericVector perCI(NumericVector theta, double conf_level);
NumericVector BCA(NumericVector theta, double conf_level);

#endif // __utils__
