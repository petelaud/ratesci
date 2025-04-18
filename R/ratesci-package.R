#' @keywords internal
"_PACKAGE"
#' ratesci: A package for computing confidence intervals for various
#' comparisons of binomial and Poisson rates.
#'
#' @section ratesci functions:
#'   * scoreci(): for score-based confidence intervals
#'   * scasci(): wrapper function to compute SCAS interval
#'   * tdasci(): wrapper function to compute TDAS random effects stratified
#'    interval
#'   * moverci(): for the MOVER method
#'   * moverbci(): wrapper function to compute MOVER-B interval
#'   * jeffreysci(): wrapper function to compute Jeffreys interval for a
#'   single rate
#'   * scaspci(): non-iterative SCAS method for a single rate
#'   * rateci(): wrapper function for SCAS, Jeffreys or 'exact' methods
#'   for a single rate
#'   * pairbinci(): for paired binomial data (includes asymptotic score and
#'   MOVER options)
#'   * clusterpci(): for estimation of binomial proportions based on
#'   clustered data
#' @name ratesci-package
#' @references
#' Laud PJ. Equal-tailed confidence intervals for comparison of
#' rates. Pharmaceutical Statistics 2017; 16:334-348.
#'
#' Laud PJ. Corrigendum: Equal-tailed confidence intervals for comparison of
#' rates. Pharmaceutical Statistics 2018; 17:290-293.
#'
#' Tang Y. Score confidence intervals and sample sizes for stratified
#' comparisons of binomial proportions. Statistics in Medicine 2020;
#' 39:3427–3457.
#'
#' Tang Y. Comments on “Equal-tailed confidence intervals for
#' comparison of rates”. Pharmaceutical Statistics 2021;20:1288-1292.
#'
#' Laud PJ. Author's reply to the letter to the editor by Yongqiang Tang:
#' Comments on “Equal-tailed confidence intervals for
#' comparison of rates”. Pharmaceutical Statistics 2021; 20:1293-1297
#'
#' Miettinen OS, Nurminen M. Comparative analysis of two rates. Statistics in
#' Medicine 1985; 4:213-226.
#'
#' Gart JJ. Analysis of the common odds ratio: corrections for bias and
#' skewness. Bulletin of the International Statistical Institute 1985,
#' 45th session, book 1, 175-176.
#'
#' Gart JJ, Nam JM. Approximate interval estimation of the ratio of binomial
#' parameters: A review and corrections for skewness. Biometrics 1988;
#' 44(2):323-338.
#'
#' Gart JJ, Nam JM. Approximate interval estimation of the difference in
#' binomial parameters: correction for skewness and extension to multiple
#' tables. Biometrics 1990; 46(3):637-643.
#'
#' Farrington CP, Manning G. Test statistics and sample size formulae
#' for comparative binomial trials with null hypothesis of non-zero risk
#' difference or non-unity relative risk. Statistics in Medicine 1990;
#' 9(12):1447–1454.
#'
#' Newcombe RG. Interval estimation for the difference between independent
#' proportions: comparison of eleven methods. Statistics in Medicine 1998;
#' 17(8):873-890.
#'
#' Donner A, Zou G. Closed-form confidence intervals for functions of the
#' normal mean and standard deviation. Statistical Methods in Medical Research
#' 2012; 21(4):347-359.
#'
#' Tango T. Equivalence test and confidence interval for the difference
#' in proportions for the paired-sample design.
#' Statistics in Medicine 1998; 17:891-908
#'
#' Tang N-S, Tang M-L, Chan ISF. On tests of equivalence via non-unity
#' relative risk for matched-pair design.
#' Statistics in Medicine 2003; 22:1217-1233
#'
#' Laud PJ. Improved confidence intervals and tests for paired binomial
#' proportions. (2025, Under review)
#'
#' Saha K, Miller D and Wang S. A comparison of some approximate confidence
#' intervals for a single proportion for clustered binary outcome data.
#' Int J Biostat 2016; 12:1–18.

## usethis namespace: start
## usethis namespace: end
NULL
