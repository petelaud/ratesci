# Selected confidence intervals for the single binomial or Poisson rate.

Confidence intervals for the single binomial or Poisson rate. Including
SCAS and Jeffreys intervals, with or without continuity adjustment, and
'exact' Clopper-Pearson/Garwood or mid-p intervals, and another version
of the exact or mid-p interval derived from Beta distributions (from
p.115 of Brown et al.), with an equivalent using Gamma distributions for
a Poisson rate. Note that these closed-form calculations exactly match
the iterative calculations for the exact interval (when `cc = TRUE`),
but not for the mid-p interval (`cc = FALSE`) This function is
vectorised in x, n.

## Usage

``` r
rateci(
  x,
  n,
  distrib = "bin",
  level = 0.95,
  std_est = TRUE,
  cc = FALSE,
  precis = 8
)
```

## Arguments

- x:

  Numeric vector of number of events.

- n:

  Numeric vector of sample size (for binomial rate) or exposure times
  (for Poisson rate).

- distrib:

  Character string indicating distribution assumed for the input data:
  "bin" = binomial (default), "poi" = Poisson.

- level:

  Number specifying confidence level (between 0 and 1, default 0.95).

- std_est:

  logical, specifying if the crude point estimate for the proportion
  value x/n should be returned (TRUE, default) or the method-specific
  alternative point estimate consistent with a 0% confidence interval
  (FALSE).

- cc:

  Number or logical (default FALSE) specifying continuity adjustment.

- precis:

  Number (default 8) specifying precision (i.e. number of decimal
  places) to be used in root-finding subroutine for the exact confidence
  interval. (Note all other methods use closed-form calculations so are
  not affected.)

## Value

A list containing, for each method, a matrix containing lower and upper
confidence limits and point estimate of p for each value of x and n.
Methods shown depend on the cc parameter, which specifies whether the
continuity adjustment is applied to the SCAS and Jeffreys methods. The
corresponding 'exact' method is Clopper-Pearson/Garwood if cc = TRUE and
mid-p if cc = FALSE. The last list item contains details of the function
call.

## References

Laud PJ. Equal-tailed confidence intervals for comparison of rates.
Pharmaceutical Statistics 2017; 16:334-348. (Appendix A.4)

Brown LD, Cai TT and DasGupta A. Interval estimation for a binomial
proportion. Statistical Science 2001; 16(2):101-133.

Garwood F. Fiducial limits for the Poisson distribution. Biometrika
1936; 28(3-4):437, doi:10.1093/biomet/28.3-4.437.

## Author

Pete Laud, <p.j.laud@sheffield.ac.uk>
