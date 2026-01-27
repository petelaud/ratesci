# Skewness-corrected asymptotic score ("SCAS") confidence intervals for single binomial or Poisson rate using closed-form calculations.

Closed-form function for computing confidence intervals for a single
rate. Note: For associated hypothesis tests, use
[`scoreci()`](https://petelaud.github.io/ratesci/reference/scoreci.md)
with `contrast = "p"`. This function is vectorised in x, n.

## Usage

``` r
scaspci(
  x,
  n,
  distrib = "bin",
  level = 0.95,
  bcf = FALSE,
  bign = n,
  xihat = 1,
  cc = FALSE,
  ...
)
```

## Arguments

- x:

  Numeric vector of number of events.

- n:

  Numeric vector of sample sizes (for binomial rates) or exposure times
  (for Poisson rates).

- distrib:

  Character string indicating distribution assumed for the input data:  
  "bin" = binomial (default);  
  "poi" = Poisson.

- level:

  Number specifying confidence level (between 0 and 1, default 0.95).

- bcf:

  Logical (default TRUE) indicating whether to apply bias correction in
  the score denominator. Applicable to `distrib = "bin"` only.

- bign:

  Sample size N to be used in the calculation of bcf, if different
  from n. (Used by transformed SCASp method for paired conditional OR in
  [`pairbinci()`](https://petelaud.github.io/ratesci/reference/pairbinci.md).)

- xihat:

  Number specifying estimated variance inflation factor for a skewness
  corrected version of the Saha Wilson Score interval for clustered
  binomial proportions. Need to calculate using BMS and WMS as per
  Saha 2016. Used by
  [`clusterpci()`](https://petelaud.github.io/ratesci/reference/clusterpci.md)
  function for data entered per cluster.

- cc:

  Number or logical (default FALSE) specifying (amount of) continuity
  adjustment. Numeric value is taken as the gamma parameter in Laud
  2017, Appendix S2 (default 0.5 for 'conventional' adjustment if cc =
  TRUE).

- ...:

  Other arguments.

## Value

A list containing the following components:

- estimates:

  a matrix containing estimated rate(s), the SCAS confidence interval,
  and the input values x and n.

- call:

  details of the function call.

## References

Laud PJ. Equal-tailed confidence intervals for comparison of rates.
Pharmaceutical Statistics 2017; 16:334-348. (Appendix A.4)

## Author

Pete Laud, <p.j.laud@sheffield.ac.uk>
