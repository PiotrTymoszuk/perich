#ifndef __test__
#define __test__

#include <Rcpp.h>

using namespace Rcpp;

// testing: a function that executes a single algorithm iteration

NumericMatrix testIter(NumericVector obs,
                       IntegerVector f,
                       NumericVector wt,
                       double laplace);

// testing a helper function that executes n iterations of the algorithm

NumericMatrix testMultiIter(NumericVector obs,
                            IntegerVector f,
                            NumericVector wt,
                            String alternative,
                            int n_iter,
                            double laplace,
                            String ci_type,
                            double ci_level);

#endif // __utils__
