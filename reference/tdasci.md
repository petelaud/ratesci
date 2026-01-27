# t-distribution asymptotic score ("TDAS") confidence intervals for random effects stratified comparisons of independent binomial or Poisson rates.

Wrapper function for the TDAS method. Score-based stratified confidence
intervals for the rate (or risk) difference ("RD") or ratio ("RR") for
independent binomial or Poisson rates, or for odds ratio ("OR", binomial
only), or for prevalence or incidence rate ("p"). This function combines
vector inputs into a single stratified random effects analysis (e.g.
meta-analysis), incorporating any stratum variability into the
confidence interval.

## Usage

``` r
tdasci(
  x1,
  n1,
  x2 = NULL,
  n2 = NULL,
  distrib = "bin",
  contrast = "RD",
  level = 0.95,
  cc = FALSE,
  theta0 = NULL,
  precis = 6,
  plot = FALSE,
  hetplot = FALSE,
  plotmax = 100,
  xlim = NULL,
  ylim = NULL,
  weighting = NULL,
  mn_tol = 1e-08,
  MNtol = NULL,
  wt = NULL,
  skew = TRUE,
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
  95\\ excludes theta0. By default, a two-sided test against theta0 = 0
  (for RD) or 1 (for RR/OR) is also output.

- precis:

  Number (default 6) specifying precision (i.e. number of decimal
  places) to be used in optimisation subroutine for the confidence
  interval.

- plot:

  Logical (default FALSE) indicating whether to output plot of the score
  function

- hetplot:

  Logical (default FALSE) indicating whether to output plots for
  evaluating heterogeneity of stratified datasets.

- plotmax:

  Numeric value indicating maximum value to be displayed on x-axis of
  plots (useful for ratio contrasts which can be infinite).

- xlim:

  pair of values indicating range of values to be plotted.

- ylim:

  pair of values indicating range of values to be plotted.

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

- skew:

  Logical (default TRUE) indicating whether to apply skewness correction
  (for the SCAS method recommended in Laud 2017) or not (for the
  Miettinen-Nurminen method) to the per-stratum estimates provided in
  the output. Has no effect on the TDAS interval itself.

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

  a matrix containing estimates of the rates in each group and of the
  requested contrast, with its confidence interval

- pval:

  a matrix containing details of the corresponding 2-sided significance
  test against the null hypothesis that p_1 = p_2, and one-sided
  significance tests against the null hypothesis that theta \>= or \<=
  theta0

- Qtest:

  a vector of values describing and testing heterogeneity

- weighting:

  a string indicating the selected weighting method

- stratdata:

  a matrix containing stratum estimates and weights

- call:

  details of the function call

## References

Laud PJ. Equal-tailed confidence intervals for comparison of rates.
Pharmaceutical Statistics 2017; 16:334-348.

Laud PJ. Corrigendum: Equal-tailed confidence intervals for comparison
of rates. Pharmaceutical Statistics 2018; 17:290-293.

## Author

Pete Laud, <p.j.laud@sheffield.ac.uk>
