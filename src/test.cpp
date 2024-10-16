/*** Enrichment testing */

#include <Rcpp.h>
#include <Rmath.h>

#include "utils.h"
#include "p_values.h"
#include "test.h"

using namespace Rcpp;

// enrichment testing: a function executing a single algorithm iteration
// the function returns a five-row numeric matrix:
// 1) enrichment scores (actually Laplace-corrected odds ratios),
// 2) 0/1 indexes of meeting the H0 hypothesis for the 'greater' scenario
// (i.e. N observed > N expected; 1: H0 met, 0: H0 falsified),
// 3) 0/1 indexes for the 'less' scenario (N observed < N expected),
// 4) observed counts in the categories defined by f,
// 5) expected counts obtained via permutations

NumericMatrix testIter(NumericVector obs,
                       IntegerVector f,
                       NumericVector wt,
                       double laplace = 1) {

  // no special entry control: done by the upstream function
  //
  // creating a weighted random permutation (sampling without replacement)
  // of the obs vector;
  // the weights are sample-specific event probabilities.

  int n_samples = obs.size();

  NumericVector permVec = sample(obs, n_samples, false, wt);

  // splits of the vector with observed events and expected events;
  // comparing the frequencies in the groups

  List obsSplits = splitNumVec(obs, f);
  List permSplits = splitNumVec(permVec, f);

  // comparison of event sums in the splits of the observed event and the
  // permuted event vectors

  int n_groups = obsSplits.size();

  NumericVector obsSums(n_groups);
  NumericVector permSums(n_groups);

  NumericVector enrichScores(n_groups);
  NumericVector greaterHits(n_groups, 1.0);
  NumericVector lessHits(n_groups, 1.0);

  for(int i = 0; i < n_groups; ++i) {

    NumericVector obsSubset = obsSplits[i];
    NumericVector permSubset = permSplits[i];

    obsSums[i] = sum(na_omit(obsSubset));
    permSums[i] = sum(na_omit(permSubset));

    enrichScores[i] = (obsSums[i] + laplace)/(permSums[i] + laplace);

    if(obsSums[i] > permSums[i]) greaterHits[i] = 0.0;
    if(obsSums[i] < permSums[i]) lessHits[i] = 0.0;

  }

  NumericMatrix resMat(5, n_groups);

  resMat(0, _) = enrichScores;
  resMat(1, _) = greaterHits;
  resMat(2, _) = lessHits;
  resMat(3, _) = obsSums;
  resMat(4, _) = permSums;

  return resMat;

}

// enrichment testing: multi-iteration testing function.
// it returns a seven column matrix:
// 1) mean enrichment scores in categories of f,
// 2) lower and 3) upper confidence interval of the enrichment score,
// 4) p value for enrichment, depletion or both, in particular categories,
// 5) stores the observed numbers of events,
// 6) stores the expected numbers of events obtained by permutations,
// 7) numbers of samples in the categories,
// 8) Yates-corrected value of the chi-square test statistic of differences
// between the observed and expected events,
// 9) numbers of degrees of freedom,
// 10) Cramer's V effect size statistic for the global difference
// 11) p value of global difference between the observed and expected number of
// events.
// Categories of f are represented by rows

NumericMatrix testMultiIter(NumericVector obs,
                            IntegerVector f,
                            NumericVector wt,
                            String alternative = "both",
                            int n_iter = 1000,
                            double laplace = 1,
                            String ci_type = "bca",
                            double ci_level = 0.95) {

  // declaration of the output matrix

  int nCat = unique(na_omit(f)).length();

  NumericMatrix resMat(nCat, 11);
  CharacterVector colNames = {"ES", "lower_ci", "upper_ci",
                              "p_value",
                              "n_observed", "n_expected", "n_total",
                              "chisq", "df", "cramer_v", "global_p_value"};

  colnames(resMat) = colNames;

  // no entry control; it will be done by an upstream function
  // but we need to select complete cases

  LogicalVector complete_obs = is_na(obs);
  LogicalVector complete_f = is_na(f);
  LogicalVector complete_wt = is_na(wt);

  LogicalVector complete_cases = !(complete_f | complete_f | complete_wt);

  if(sum(complete_cases) < 2) {

    warning("Not enough observations");

    return resMat;

  }

  obs = obs[complete_cases];
  f = f[complete_cases];
  wt = wt[complete_cases];

  double n_samples = 1.0 * obs.size();

  // enrichment scores and H0 hypothesis hits in weighted permutations

  NumericMatrix esMat(n_iter, nCat);
  NumericMatrix greaterMat(n_iter, nCat);
  NumericMatrix lessMat(n_iter, nCat);
  NumericMatrix obsMat(n_iter, nCat);
  NumericMatrix permMat(n_iter, nCat);

  NumericMatrix iterRes(5, nCat);

  for(int i = 0; i < n_iter; ++i) {

    iterRes = testIter(obs, f, wt, laplace);

    esMat(i, _) = iterRes(0, _);
    greaterMat(i, _) = iterRes(1, _);
    lessMat(i, _) = iterRes(2, _);
    obsMat(i, _) = iterRes(3, _);
    permMat(i, _) = iterRes(4, _);

  }

  // containers for numbers of hits, p values,
  // enrichment scores and their 95% confidence intervals,
  // observed and expected N numbers obtained by permutation,
  // as well as a double container for chi-square statistic

  NumericVector pVals(nCat);

  NumericVector hits(n_iter);

  NumericVector (*ciFun)(NumericVector, double);

  if(ci_type == "bca") {

    ciFun = BCA;

  } else {

    ciFun = perCI;

  }

  NumericVector esMeans(nCat);
  NumericMatrix esCI(2, nCat);

  NumericVector esScores(n_iter);

  NumericVector obsMeans(nCat);
  NumericVector permMeans(nCat);

  NumericVector observed(n_iter);
  NumericVector expected(n_iter);

  double df = 1.0 * (nCat - 1);
  double chisq = 0.0;

  // calculation of the statistics

  for(int i = 0; i < nCat; ++i) {

    // enrichment scores with confidence intervals

    esScores = esMat(_, i);

    esMeans[i] = sum(esScores)/n_iter;

    esCI(_, i) = ciFun(esScores, ci_level);

    // p values

    if(alternative == "greater") {

      hits = greaterMat(_, i);

    } else if(alternative == "less") {

      hits = lessMat(_, i);

    } else {

      if(esMeans[i] > 1) {

        hits = greaterMat(_, i);

      } else {

        hits = lessMat(_, i);

      }

    }

    pVals[i] = (sum(hits) + 1)/n_iter;

    if(pVals[i] > 1.0) pVals[i] = 1.0;

    // observed and expected N numbers of events

    observed = obsMat(_, i);
    expected = permMat(_, i);

    obsMeans[i] = sum(observed)/n_iter;
    permMeans[i] = sum(expected)/n_iter;

    // chi-square statistic

    chisq +=
      std::pow(std::abs(obsMeans[i] - permMeans[i]) - 0.5, 2.0)/permMeans[i];

  }

  // Cramer's V and p value for the global difference

  NumericVector categoryN = as<NumericVector>(table(f));

  double cramer_denom = min(categoryN) - 1.0;

  if(df < cramer_denom) cramer_denom = df;

  double cramer_v = std::sqrt(chisq/n_samples/cramer_denom);

  double global_p = R::pchisq(chisq, df, false, false);

  NumericVector chisqVec(nCat, chisq);
  NumericVector dfVec(nCat, df);
  NumericVector cramerVec(nCat, cramer_v);
  NumericVector pGlobalVec(nCat, global_p);

  // the output matrix

  resMat(_, 0) = esMeans;
  resMat(_, 1) = esCI(0, _);
  resMat(_, 2) = esCI(1, _);
  resMat(_, 3) = pVals;
  resMat(_, 4) = obsMeans;
  resMat(_, 5) = permMeans;
  resMat(_, 6) = categoryN;
  resMat(_, 7) = chisqVec;
  resMat(_, 8) = dfVec;
  resMat(_, 9) = cramerVec;
  resMat(_, 10) = pGlobalVec;

  return resMat;

}

// enrichment testing for a single vector of observed items and the splitting
// vector f. This function will be exposed to R.
// It returns a four column matrix: 1) mean enrichment scores in categories
// of f, 2) lower and 3) upper confidence interval of the enrichment score, and
// 4) enrichment p value. Categories of f are represented by rows

//[[Rcpp::export]]

NumericMatrix testVector(NumericVector x,
                         IntegerVector f,
                         Nullable<NumericMatrix> background = R_NilValue,
                         Nullable<NumericVector> weights = R_NilValue,
                         String alternative = "both",
                         int n_iter = 1000,
                         double laplace = 1,
                         String ci_type = "bca",
                         double ci_level = 0.95) {

  // entry control

  int n_samples = x.length();

  if(f.length() != n_samples) stop("Incompatible lengths of 'x' and 'f'");

  NumericVector wt;

  if(background.isNotNull()) {

    NumericMatrix bg(background);

    wt = pSample(bg, 1e-6);

  } else if(weights.isNotNull()) {

    wt = weights;

  } else {

   stop("Either 'background' or 'weights' need to be provided.");

  }

  if(wt.length() != n_samples) {

    stop("Length of 'x' must be equal to the number of rows in 'background'");

  }

  int nCat = unique(f).length();

  if(nCat == 1) stop("'f' has to have at least two categories");

  // testing

  NumericMatrix testMat = testMultiIter(x, f, wt, alternative,
                                        n_iter, laplace,
                                        ci_type, ci_level);

  return testMat;

}

// enrichment testing for a binary matrix of events; the function returns
// a list of matrices as described for the vector-variant

//[[Rcpp::export]]

List testMatrix(NumericMatrix x,
                IntegerVector f,
                Nullable<NumericMatrix> background = R_NilValue,
                Nullable<NumericVector> weights = R_NilValue,
                String alternative = "both",
                int n_iter = 1000,
                double laplace = 1,
                String ci_type = "bca",
                double ci_level = 0.95) {

  // entry control

  int n_samples = x.nrow();
  int n_features = x.ncol();

  if(f.length() != n_samples) stop("Incompatible lengths of 'x' and 'f'");

  NumericVector wt;

  if(background.isNotNull()) {

    NumericMatrix bg(background);

    wt = pSample(bg, 1e-6);

  } else if(weights.isNotNull()) {

    wt = weights;

  } else {

    stop("Either 'background' or 'weights' need to be provided.");

  }

  if(wt.length() != n_samples) {

    stop("Length of 'x' must be equal to the number of rows in 'background'");

  }

  int nCat = unique(f).length();

  if(nCat == 1) stop("'f' has to have at least two categories");

  // testing

  List resList(n_features);

  for(int i = 0; i < n_features; ++i) {

    resList[i] = testMultiIter(x(_, i),
                               f,
                               wt,
                               alternative,
                               n_iter,
                               laplace,
                               ci_type,
                               ci_level);

  }

  if(checkMtxNames(x)[1]) resList.names() = colnames(x);

  return resList;

}

// END

