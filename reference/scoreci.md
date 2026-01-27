# Score confidence intervals and tests for a single binomial or Poisson rate, or for comparisons of independent rates, with or without stratification.

Score-based confidence intervals for the rate (or risk) difference
("RD") or ratio ("RR") for independent binomial or Poisson rates, or for
odds ratio ("OR", binomial only). Including options for variance bias
correction (from Miettinen & Nurminen), skewness correction ("GNbc"
method from Laud & Dane, developed from Gart & Nam, and generalised as
"SCAS" in Laud 2017) and continuity adjustment (for strictly
conservative coverage).

Also includes score intervals for a single binomial proportion or
Poisson rate ("p"). These are based on the Wilson score interval, and
when corrected for skewness, coverage is almost identical to the mid-p
method, or to Clopper-Pearson when also continuity-adjusted.

Hypothesis tests for association or non-inferiority are provided using
the same score, to ensure consistency between test and CI. This function
is vectorised in x1, x2, n1, and n2. Vector inputs may also be combined
into a single stratified analysis (e.g. meta-analysis), either using
fixed effects, or the more general random effects "TDAS" method, which
incorporates stratum variability using a t-distribution score (inspired
by Hartung-Knapp-Sidik-Jonkman). For fixed-effects analysis of
stratified datasets, with weighting = "MH" for RD or RR, or weighting =
"INV" for OR, omitting the skewness correction produces the CMH test,
together with a coherent confidence interval for the required contrast.
Alternatively, weighting = "INV" for any contrast gives intervals
consistent with the efficient score test.

## Usage

``` r
scoreci(
  x1,
  n1,
  x2 = 0,
  n2 = 0,
  distrib = "bin",
  contrast = "RD",
  level = 0.95,
  skew = TRUE,
  simpleskew = FALSE,
  or_bias = TRUE,
  ORbias = NULL,
  rr_tang = NULL,
  RRtang = NULL,
  bcf = ifelse(contrast != "p", TRUE, FALSE),
  cc = FALSE,
  theta0 = NULL,
  precis = 6,
  plot = FALSE,
  plotmax = 100,
  hetplot = FALSE,
  xlim = NULL,
  ylim = NULL,
  stratified = FALSE,
  weighting = NULL,
  mn_tol = 1e-08,
  MNtol = NULL,
  wt = NULL,
  sda = NULL,
  fda = NULL,
  dropzeros = FALSE,
  random = FALSE,
  prediction = FALSE,
  warn = TRUE,
  ...
)
```

## Arguments

- x1, x2:

  Numeric vectors of numbers of events in group 1 & group 2
  respectively.

- n1, n2:

  Numeric vectors of sample sizes (for binomial rates) or exposure times
  (for Poisson rates) in each group.

- distrib:

  Character string indicating distribution assumed for the input data:  
  "bin" = binomial (default),  
  "poi" = Poisson.

- contrast:

  Character string indicating the contrast of interest:  
  "RD" = rate difference (default);  
  "RR" = rate ratio;  
  "OR" = odds ratio;  
  "p" gives an interval for the single proportion or rate `x1/n1`.

- level:

  Number specifying confidence level (between 0 and 1, default 0.95).

- skew:

  Logical (default TRUE) indicating whether to apply skewness correction
  (for the SCAS or Gart-Nam method) or not (for the Miettinen-Nurminen
  method).

- simpleskew:

  Logical (default FALSE) indicating whether to use the "simplified"
  skewness correction instead of the quadratic solution. See Laud 2021
  for details.  
  NOTE: this version of the score is only suitable for obtaining
  confidence limits, not p-values.

- or_bias:

  Logical (default is TRUE for `contrast = "OR"`, otherwise NULL)
  indicating whether to apply additional bias correction for OR derived
  from Gart 1985. (Laud 2018). Only applies if contrast is "OR".

- ORbias:

  (deprecated: argument renamed to or_bias.)

- rr_tang:

  Logical indicating whether to use Tang's score for RR: Stheta =
  (p1hat - p2hat \* theta) / p2d (see Tang 2020). Default TRUE for
  `stratified = TRUE`, with weighting = "IVS" or "INV". Forced to FALSE
  for `stratified = TRUE` with other weightings. Has no effect when
  `stratified = FALSE`, as p2d terms cancel out. Experimental for
  `distrib = "poi"`.

- RRtang:

  (deprecated: argument renamed to rr_tang.)

- bcf:

  Logical (default TRUE) indicating whether to apply 'N-1' variance
  correction in the score denominator. Applicable to `distrib = "bin"`
  only.  
  NOTE: `bcf = FALSE` option is really only included for legacy
  validation against previous published methods (i.e. Gart & Nam, Mee,
  or standard Chi-squared test) and for `contrast = "p"`.

- cc:

  Number or logical (default FALSE) specifying (amount of) continuity
  adjustment. Numeric value between 0 and 0.5 is taken as the gamma
  parameter in Laud 2017, Appendix S2 (`cc = TRUE` translates to 0.5 for
  'conventional' Yates adjustment).  
  IMPORTANT NOTES:

  1.  This adjustment (conventionally but controversially termed
      'continuity correction') is aimed at approximating strictly
      conservative coverage, NOT for dealing with zero cell counts. Such
      'sparse data adjustments' are not needed in the score method,
      except to deal with double-zero cells for stratified RD (&
      double-100% cells for binomial RD & RR) with IVS/INV weights.

  2.  The continuity adjustments provided here have not been fully
      tested for stratified methods, but are found to match the
      continuity-adjusted version of the Mantel-Haenszel test, when
      `cc = 0.5` for any of the binomial contrasts. Flexibility is
      included for a less conservative adjustment, such as `cc = 0.25`
      suggested in Laud 2017 (see Appendix S3.4), or
      `cc = 3/16 = 0.1875` in Mehrotra & Railkar (2000).

- theta0:

  Number to be used in a one-sided significance test (e.g.
  non-inferiority margin). 1-sided p-value will be \<0.025 iff 2-sided
  95\\ excludes theta0. (If `bcf = FALSE` and `skew = FALSE` this gives
  a Farrington-Manning test.)  
  By default, a two-sided test for association against theta0 = 0 (for
  RD) or 1 (for RR/OR) is also output:

  - If `bcf = FALSE` and `skew = FALSE` this is the same as K. Pearson's
    Chi-squared test in the single stratum case.

  - `bcf = TRUE` gives E. Pearson's 'N-1' Chi-squared test for a single
    stratum, (Recommended by Campbell 2007:
    https://doi.org/10.1002/sim.2832) and (with default weighting and
    `random = FALSE`) the CMH test for stratified tables.

  - Default `bcf = TRUE` and \`skew = TRUE produces a skewness-corrected
    version of the 'N-1' Chi-squared test or CMH. This correction will
    only change the p-value if group sizes are unequal.

- precis:

  Number (default 6) specifying precision (i.e. number of decimal
  places) to be used in optimisation subroutine for the confidence
  interval.

- plot:

  Logical (default FALSE) indicating whether to output plot of the score
  function

- plotmax:

  Numeric value indicating maximum value to be displayed on x-axis of
  plots (useful for ratio contrasts which can be infinite).

- hetplot:

  Logical (default FALSE) indicating whether to output plots for
  evaluating heterogeneity of stratified datasets.

- xlim:

  pair of values indicating range of values to be plotted.

- ylim:

  pair of values indicating range of values to be plotted.

- stratified:

  Logical (default FALSE) indicating whether to combine vector inputs
  into a single stratified analysis.  
  IMPORTANT NOTE: The mechanism for stratified calculations is enabled
  for contrast = "p", but the performance of the resulting intervals has
  not been fully evaluated.

- weighting:

  String indicating which weighting method to use if stratified =
  "TRUE":  
  "IVS" = Inverse Variance of Score (see Laud 2017 for details);  
  "INV" = Inverse Variance (bcf omitted, default for contrast = "OR"
  giving CMH test);  
  "MH" = Mantel-Haenszel (n1j \* n2j) / (n1j + n2j) (default for
  contrast = "RD" or "RR" giving CMH test); (= sample size for contrast
  = "p");  
  "MN" = Miettinen-Nurminen weights. (similar to MH for contrast = "RD"
  or "RR", similar to INV for contrast = "OR");  
  "Tang" = (n1j \* n2j) / (n1j + n2j) / (1 - pj) from Tang 2020, for an
  optimal test of RD if RRs are constant across strata. (Included only
  for validation purposes. In general, such a test would more logically
  use contrast = "RR" with weighting = "INV") For CI consistent with a
  CMH test, select `skew = FALSE`, `random = FALSE`, and use default MH
  weighting for RD/RR and INV for OR.  
  `Weighting = "MN"` also matches the CMH test.  
  For the Radhakrishna optimal (most powerful) test, select INV
  weighting.  
  Note: Alternative user-specified weighting may also be applied, via
  the 'wt' argument.

- mn_tol:

  Numeric value indicating convergence tolerance to be used in iteration
  with weighting = "MN".

- MNtol:

  (deprecated: argument renamed to mn_tol)

- wt:

  Numeric vector containing (optional) user-specified weights.  
  Overrides `weighting` if non-empty.

- sda:

  Sparse data adjustment to avoid zero variance when `x1 + x2 = 0`: Only
  applied when `stratified = TRUE`. Default 0.5 for RD with IVS/INV
  weights. Not required for RR/OR, default is to remove double-zero
  strata instead.

- fda:

  Full data adjustment to avoid zero variance when x1 + x2 = n1 + n2:
  Only applied when `stratified = TRUE`. Default 0.5 for RD & RR with
  IVS/INV weights. Not required for OR, default is to remove affected
  strata.

- dropzeros:

  Logical (default FALSE) indicating whether to drop uninformative
  strata for RR/OR (i.e. strata with `x1 + x2 = 0`), even when the
  choice of weights would allow them to be retained for a fixed effects
  analysis. Has no effect on estimates, just the heterogeneity test.

- random:

  Logical (default FALSE) indicating whether to perform random effects
  meta-analysis for stratified data, using the t-distribution (TDAS)
  method for stratified data (defined in Laud 2017).  
  NOTE: If `random = TRUE`, then `skew = TRUE` only affects the
  per-stratum estimates.

- prediction:

  Logical (default FALSE) indicating whether to produce a prediction
  interval (work in progress).

- warn:

  Logical (default TRUE) giving the option to suppress warnings.

- ...:

  Other arguments.

## Value

A list containing the following components:

- estimates:

  a matrix containing estimates of the requested contrast and its
  confidence interval, and the estimated rates in each group: (p1hat,
  p2hat) are (r1, r0) from Miettinen-Nurminen, or (r1\*, r0\*) when
  stratified; (p1mle, p2mle) are (R1, R0), or (R1\*, R0\*) when
  stratified, evaluated at the MLE for the contrast parameter,
  incorporating any specified skewness/bias corrections.

- pval:

  a matrix containing details of the corresponding 2-sided significance
  test against the null hypothesis that `p_1 = p_2`, and one-sided
  significance tests against the null hypothesis that theta \>= or \<=
  theta0.

- call:

  details of the function call.

If `stratified = TRUE`, the following outputs are added:

- Qtest:

  a vector of values describing and testing heterogeneity, including a
  score-based version of a Q statistic and p-value, I^2 and tau^2 to
  quantify heterogeneity, and a test for qualitative interaction
  analogous to the Gail and Simon test.

- weighting:

  a string indicating the selected weighting method.

- stratdata:

  a matrix containing stratum estimates and weights.

## References

Laud PJ. Equal-tailed confidence intervals for comparison of rates.
Pharmaceutical Statistics 2017; 16:334-348.

Laud PJ. Corrigendum: Equal-tailed confidence intervals for comparison
of rates. Pharmaceutical Statistics 2018; 17:290-293.

Laud PJ, Dane A. Confidence intervals for the difference between
independent binomial proportions: comparison using a graphical approach
and moving averages. Pharmaceutical Statistics 2014; 13(5):294-308.

Miettinen OS, Nurminen M. Comparative analysis of two rates. Statistics
in Medicine 1985; 4:213-226.

Farrington CP, Manning G. Test statistics and sample size formulae for
comparative binomial trials with null hypothesis of non-zero risk
difference or non-unity relative risk. Statistics in Medicine 1990;
9(12):1447-1454.

Gart JJ. Analysis of the common odds ratio: corrections for bias and
skewness. Bulletin of the International Statistical Institute 1985, 45th
session, book 1, 175-176.

Gart JJ, Nam Jm. Approximate interval estimation of the ratio of
binomial parameters: a review and corrections for skewness. Biometrics
1988; 44(2):323-338.

Gart JJ, Nam Jm. Approximate interval estimation of the difference in
binomial parameters: correction for skewness and extension to multiple
tables. Biometrics 1990; 46(3):637-643.

Tang Y. Score confidence intervals and sample sizes for stratified
comparisons of binomial proportions. Statistics in Medicine 2020;
39:3427-3457.

## Author

Pete Laud, <p.j.laud@sheffield.ac.uk>

## Examples

``` r
# Binomial RD, SCAS method:
scoreci(
  x1 = c(12, 19, 5), n1 = c(16, 29, 56),
  x2 = c(1, 22, 0), n2 = c(16, 30, 29)
)
#> $estimates
#>            lower         est     upper level x1 n1 x2 n2      p1hat     p2hat
#> [1,]  0.38567597  0.68162256 0.8778839  0.95 12 16  1 16 0.75000000 0.0625000
#> [2,] -0.31189939 -0.07795525 0.1600794  0.95 19 29 22 30 0.65517241 0.7333333
#> [3,] -0.01861751  0.09168753 0.1867180  0.95  5 56  0 29 0.08928571 0.0000000
#>           p1mle      p2mle
#> [1,] 0.74553163 0.06390906
#> [2,] 0.65528437 0.73323962
#> [3,] 0.09168753 0.00000000
#> 
#> $pval
#>           chisq   pval2sided theta0  scorenull pval_left   pval_right
#> [1,] 15.1862348 9.741092e-05      0  3.8969520 0.9999513 4.870546e-05
#> [2,]  0.4181622 5.178555e-01      0 -0.6466546 0.2589278 7.410722e-01
#> [3,]  3.0248621 8.199729e-02      0  1.7392131 0.9590014 4.099865e-02
#> 
#> $call
#>  distrib contrast    level      bcf     skew       cc 
#>    "bin"     "RD"   "0.95"   "TRUE"   "TRUE"  "FALSE" 
#> 

# Binomial RD, MN method:
scoreci(
  x1 = c(12, 19, 5), n1 = c(16, 29, 56),
  x2 = c(1, 22, 0), n2 = c(16, 30, 29), skew = FALSE
)
#> $estimates
#>            lower         est    upper level x1 n1 x2 n2      p1hat     p2hat
#> [1,]  0.37497827  0.68749997 0.862899  0.95 12 16  1 16 0.75000000 0.0625000
#> [2,] -0.30859360 -0.07816094 0.158172  0.95 19 29 22 30 0.65517241 0.7333333
#> [3,] -0.03259659  0.08928570 0.193331  0.95  5 56  0 29 0.08928571 0.0000000
#>          p1mle      p2mle
#> [1,] 0.7500000 0.06250001
#> [2,] 0.6551724 0.73333334
#> [3,] 0.0892857 0.00000000
#> 
#> $pval
#>           chisq   pval2sided theta0  scorenull pval_left   pval_right
#> [1,] 15.1862348 9.741092e-05      0  3.8969520 0.9999513 4.870546e-05
#> [2,]  0.4177055 5.180842e-01      0 -0.6463014 0.2590421 7.409579e-01
#> [3,]  2.7187500 9.917565e-02      0  1.6488632 0.9504122 4.958783e-02
#> 
#> $call
#>  distrib contrast    level      bcf     skew       cc 
#>    "bin"     "RD"   "0.95"   "TRUE"  "FALSE"  "FALSE" 
#> 

# Poisson RR, SCAS method:
scoreci(x1 = 5, n1 = 56, x2 = 0, n2 = 29, distrib = "poi", contrast = "RR")
#> $estimates
#>          lower      est upper level x1 n1 x2 n2      p1hat p2hat      p1mle
#> [1,] 0.7264486 16.05357   Inf  0.95  5 56  0 29 0.08928571     0 0.08649554
#>            p2mle
#> [1,] 0.005387931
#> 
#> $pval
#>         chisq pval2sided theta0 scorenull pval_left pval_right
#> [1,] 2.904393 0.08833848      1  1.704228 0.9558308 0.04416924
#> 
#> $call
#>  distrib contrast    level      bcf     skew  rr_tang       cc 
#>    "poi"     "RR"   "0.95"   "TRUE"   "TRUE"   "TRUE"  "FALSE" 
#> 

# Poisson RR, MN method:
scoreci(
  x1 = 5, n1 = 56, x2 = 0, n2 = 29, distrib = "poi",
  contrast = "RR", skew = FALSE
)
#> $estimates
#>          lower est upper level x1 n1 x2 n2      p1hat p2hat p1mle p2mle
#> [1,] 0.6740371 Inf   Inf  0.95  5 56  0 29 0.08928571     0   NaN 1e-08
#> 
#> $pval
#>         chisq pval2sided theta0 scorenull pval_left pval_right
#> [1,] 2.589286  0.1075888      1  1.609126 0.9462056 0.05379442
#> 
#> $call
#>  distrib contrast    level      bcf     skew  rr_tang       cc 
#>    "poi"     "RR"   "0.95"   "TRUE"  "FALSE"   "TRUE"  "FALSE" 
#> 

# Binomial rate, SCAS method:
scoreci(x1 = c(5, 0), n1 = c(56, 29), contrast = "p")
#> $estimates
#>           lower         est      upper level x1 n1      p1hat       p1mle
#> [1,] 0.03396264 0.091715962 0.18585265  0.95  5 56 0.08928571 0.091715962
#> [2,] 0.00000000 0.005681813 0.09170711  0.95  0 29 0.00000000 0.005681813
#> 
#> $pval
#>         chisq   pval2sided theta0 scorenull    pval_left pval_right
#> [1,] 37.78571 7.895787e-10    0.5 -6.147009 3.947893e-10          1
#> [2,] 29.00000 7.237830e-08    0.5 -5.385165 3.618915e-08          1
#> 
#> $call
#>  distrib contrast    level      bcf     skew       cc 
#>    "bin"      "p"   "0.95"  "FALSE"   "TRUE"  "FALSE" 
#> 

# Binomial rate, Wilson score method:
scoreci(x1 = c(5, 0), n1 = c(56, 29), contrast = "p", skew = FALSE)
#> $estimates
#>           lower       est     upper level x1 n1      p1hat     p1mle
#> [1,] 0.03874213 0.0892857 0.1925600  0.95  5 56 0.08928571 0.0892857
#> [2,] 0.00000000 0.0000000 0.1169698  0.95  0 29 0.00000000 0.0000000
#> 
#> $pval
#>         chisq   pval2sided theta0 scorenull    pval_left pval_right
#> [1,] 37.78571 7.895787e-10    0.5 -6.147009 3.947893e-10          1
#> [2,] 29.00000 7.237830e-08    0.5 -5.385165 3.618915e-08          1
#> 
#> $call
#>  distrib contrast    level      bcf     skew       cc 
#>    "bin"      "p"   "0.95"  "FALSE"  "FALSE"  "FALSE" 
#> 

# Poisson rate, SCAS method:
scoreci(x1 = c(5, 0), n1 = c(56, 29), distrib = "poi", contrast = "p")
#> $estimates
#>           lower         est      upper level x1 n1      p1hat       p1mle
#> [1,] 0.03314556 0.092261900 0.19710988  0.95  5 56 0.08928571 0.092261900
#> [2,] 0.00000000 0.005747108 0.09705604  0.95  0 29 0.00000000 0.005747108
#> 
#> $pval
#>         chisq   pval2sided theta0 scorenull    pval_left pval_right
#> [1,] 26.52974 2.595125e-07    0.5 -5.150703 1.297562e-07  0.9999999
#> [2,] 22.58937 2.005914e-06    0.5 -4.752828 1.002957e-06  0.9999990
#> 
#> $call
#>  distrib contrast    level      bcf     skew       cc 
#>    "poi"      "p"   "0.95"  "FALSE"   "TRUE"  "FALSE" 
#> 

# Stratified example, using data from Hartung & Knapp:
scoreci(
  x1 = c(15, 12, 29, 42, 14, 44, 14, 29, 10, 17, 38, 19, 21),
  x2 = c(9, 1, 18, 31, 6, 17, 7, 23, 3, 6, 12, 22, 19),
  n1 = c(16, 16, 34, 56, 22, 54, 17, 58, 14, 26, 44, 29, 38),
  n2 = c(16, 16, 34, 56, 22, 55, 15, 58, 15, 27, 45, 30, 38),
  stratified = TRUE
)
#> $estimates
#>          lower       est     upper level     p1hat     p2hat     p1mle
#> [1,] 0.2455228 0.3087549 0.3703304  0.95 0.7168521 0.4079934 0.7143654
#>          p2mle
#> [1,] 0.4056105
#> 
#> $pval
#>         chisq  pval2sided theta0 scorenull pval_left pval_right
#> [1,] 84.63601 3.58673e-20      0  9.199783         1          0
#> 
#> $Qtest
#>            Q         Q_df     pval_het           I2         tau2           Qc 
#> 4.434094e+01 1.200000e+01 1.335827e-05 7.293697e+01 3.249909e-02 4.381431e-01 
#> pval_qualhet 
#> 9.874622e-01 
#> 
#> $weighting
#> [1] "MH"
#> 
#> $stratdata
#>       x1j n1j x2j n2j    p1hatj    p2hatj  wt_fixed wtpct_fixed wtpct_rand
#>  [1,]  15  16   9  16 0.9375000 0.5625000  8.000000    3.761235   3.761235
#>  [2,]  12  16   1  16 0.7500000 0.0625000  8.000000    3.761235   3.761235
#>  [3,]  29  34  18  34 0.8529412 0.5294118 17.000000    7.992625   7.992625
#>  [4,]  42  56  31  56 0.7500000 0.5535714 28.000000   13.164324  13.164324
#>  [5,]  14  22   6  22 0.6363636 0.2727273 11.000000    5.171699   5.171699
#>  [6,]  44  54  17  55 0.8148148 0.3090909 27.247706   12.810630  12.810630
#>  [7,]  14  17   7  15 0.8235294 0.4666667  7.968750    3.746543   3.746543
#>  [8,]  29  58  23  58 0.5000000 0.3965517 29.000000   13.634478  13.634478
#>  [9,]  10  14   3  15 0.7142857 0.2000000  7.241379    3.404567   3.404567
#> [10,]  17  26   6  27 0.6538462 0.2222222 13.245283    6.227328   6.227328
#> [11,]  38  44  12  45 0.8636364 0.2666667 22.247191   10.459615  10.459615
#> [12,]  19  29  22  30 0.6551724 0.7333333 14.745763    6.932786   6.932786
#> [13,]  21  38  19  38 0.5526316 0.5000000 19.000000    8.932934   8.932934
#>           theta_j     lower_j   upper_j         V_j    Stheta_j         Q_j
#>  [1,]  0.37432668  0.07947651 0.6366192 0.019934516  0.06624508  0.22014132
#>  [2,]  0.68162256  0.38567597 0.8778839 0.026998683  0.37874508  5.31314199
#>  [3,]  0.32258016  0.10699961 0.5210504 0.011266832  0.01477449  0.01937418
#>  [4,]  0.19597310  0.01948127 0.3650952 0.007443573 -0.11232635  1.69504736
#>  [5,]  0.36101505  0.06819779 0.6120265 0.020831237  0.05488144  0.14458924
#>  [6,]  0.50425580  0.33270276 0.6528154 0.008184236  0.19696898  4.74042804
#>  [7,]  0.35432470  0.02152622 0.6408332 0.026682150  0.04810782  0.08673824
#>  [8,]  0.10316080 -0.07844922 0.2802846 0.007795273 -0.20530665  5.40722789
#>  [9,]  0.50878373  0.15475377 0.7736135 0.032012523  0.20553079  1.31957442
#> [10,]  0.42909852  0.16911775 0.6487974 0.017068011  0.12286901  0.88450806
#> [11,]  0.59488156  0.41394266 0.7429163 0.009980466  0.28821478  8.32303392
#> [12,] -0.07795525 -0.31189939 0.1600794 0.013943229 -0.38691584 10.73667113
#> [13,]  0.05240601 -0.17214090 0.2734169 0.012035530 -0.25612334  5.45045942
#> 
#> $call
#>  distrib contrast    level      bcf     skew       cc   random 
#>    "bin"     "RD"   "0.95"   "TRUE"   "TRUE"  "FALSE"  "FALSE" 
#> 

# "Random effects" TDAS example, using data from Hartung & Knapp:
scoreci(
  x1 = c(15, 12, 29, 42, 14, 44, 14, 29, 10, 17, 38, 19, 21),
  x2 = c(9, 1, 18, 31, 6, 17, 7, 23, 3, 6, 12, 22, 19),
  n1 = c(16, 16, 34, 56, 22, 54, 17, 58, 14, 26, 44, 29, 38),
  n2 = c(16, 16, 34, 56, 22, 55, 15, 58, 15, 27, 45, 30, 38),
  stratified = TRUE, random = TRUE
)
#> $estimates
#>          lower       est     upper level     p1hat     p2hat    p1mle     p2mle
#> [1,] 0.1746892 0.3088587 0.4430283  0.95 0.7168521 0.4079934 0.714412 0.4055533
#> 
#> $pval
#>         chisq   pval2sided theta0 scorenull pval_left   pval_right
#> [1,] 25.15659 0.0003013239      0  5.015634 0.9998493 0.0001506619
#> 
#> $Qtest
#>            Q         Q_df     pval_het           I2         tau2           Qc 
#> 4.434094e+01 1.200000e+01 1.335827e-05 7.293697e+01 3.249909e-02 4.381431e-01 
#> pval_qualhet 
#> 9.874622e-01 
#> 
#> $weighting
#> [1] "MH"
#> 
#> $stratdata
#>       x1j n1j x2j n2j    p1hatj    p2hatj  wt_fixed wtpct_fixed wtpct_rand
#>  [1,]  15  16   9  16 0.9375000 0.5625000  8.000000    3.761235   3.761235
#>  [2,]  12  16   1  16 0.7500000 0.0625000  8.000000    3.761235   3.761235
#>  [3,]  29  34  18  34 0.8529412 0.5294118 17.000000    7.992625   7.992625
#>  [4,]  42  56  31  56 0.7500000 0.5535714 28.000000   13.164324  13.164324
#>  [5,]  14  22   6  22 0.6363636 0.2727273 11.000000    5.171699   5.171699
#>  [6,]  44  54  17  55 0.8148148 0.3090909 27.247706   12.810630  12.810630
#>  [7,]  14  17   7  15 0.8235294 0.4666667  7.968750    3.746543   3.746543
#>  [8,]  29  58  23  58 0.5000000 0.3965517 29.000000   13.634478  13.634478
#>  [9,]  10  14   3  15 0.7142857 0.2000000  7.241379    3.404567   3.404567
#> [10,]  17  26   6  27 0.6538462 0.2222222 13.245283    6.227328   6.227328
#> [11,]  38  44  12  45 0.8636364 0.2666667 22.247191   10.459615  10.459615
#> [12,]  19  29  22  30 0.6551724 0.7333333 14.745763    6.932786   6.932786
#> [13,]  21  38  19  38 0.5526316 0.5000000 19.000000    8.932934   8.932934
#>           theta_j     lower_j   upper_j         V_j    Stheta_j         Q_j
#>  [1,]  0.37432668  0.07947651 0.6366192 0.019934516  0.06624508  0.22014132
#>  [2,]  0.68162256  0.38567597 0.8778839 0.026998683  0.37874508  5.31314199
#>  [3,]  0.32258016  0.10699961 0.5210504 0.011266832  0.01477449  0.01937418
#>  [4,]  0.19597310  0.01948127 0.3650952 0.007443573 -0.11232635  1.69504736
#>  [5,]  0.36101505  0.06819779 0.6120265 0.020831237  0.05488144  0.14458924
#>  [6,]  0.50425580  0.33270276 0.6528154 0.008184236  0.19696898  4.74042804
#>  [7,]  0.35432470  0.02152622 0.6408332 0.026682150  0.04810782  0.08673824
#>  [8,]  0.10316080 -0.07844922 0.2802846 0.007795273 -0.20530665  5.40722789
#>  [9,]  0.50878373  0.15475377 0.7736135 0.032012523  0.20553079  1.31957442
#> [10,]  0.42909852  0.16911775 0.6487974 0.017068011  0.12286901  0.88450806
#> [11,]  0.59488156  0.41394266 0.7429163 0.009980466  0.28821478  8.32303392
#> [12,] -0.07795525 -0.31189939 0.1600794 0.013943229 -0.38691584 10.73667113
#> [13,]  0.05240601 -0.17214090 0.2734169 0.012035530 -0.25612334  5.45045942
#> 
#> $call
#>  distrib contrast    level      bcf     skew       cc   random 
#>    "bin"     "RD"   "0.95"   "TRUE"   "TRUE"  "FALSE"   "TRUE" 
#> 

# Stratified example, with extremely rare instance of non-calculable skewness
# correction seen on plot of score function:
if (FALSE) { # interactive()
scoreci(
  x1 = c(1, 16), n1 = c(20, 40), x2 = c(0, 139), n2 = c(80, 160),
  contrast = "RD", skew = TRUE, simpleskew = FALSE,
  distrib = "bin", stratified = TRUE, plot = TRUE, weighting = "IVS"
)
}
```
