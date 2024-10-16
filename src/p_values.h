#ifndef __p_values__
#define __p_values__

#include <Rcpp.h>

using namespace Rcpp;

// feature- and sample-wise probabilities of an event

NumericVector pFeature (NumericMatrix x, double min_p);
NumericVector pSample (NumericMatrix x, double min_p);

#endif // __p_values__
