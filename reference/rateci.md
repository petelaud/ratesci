# Selected confidence intervals for the single binomial or Poisson rate.

Confidence intervals for the single binomial or Poisson rate. This
convenience wrapper function produces a selection of alternative
methods. The first three are recommended for achieving 1-sided and
2-sided coverage probability close to the nominal levels (see Laud 2017
and Laud 2018):

- SCAS (skewness-corrected asymptotic score)

- Jeffreys

- mid-p (two versions, using exact calculation or approximation via
  Beta/Gamma distribution, see p.115 of Brown et al.)) The following
  more approximate methods are included for users wishing to use a more
  established or commonly used method:

- Wilson score

- Agresti-Coull

- Wald (strongly advise this is not used for any purpose but included
  for reference)

All methods can be made more conservative with a 'continuity
adjustment', which may either be specified as TRUE, or an intermediate
'compromise' value between 0 and 0.5 may be selected. When `cc` is
`TRUE` or `0.5`, the mid-p method becomes the Clopper-Pearson interval
(or Garwood for Poisson rates). Note that Brown et al's Beta formulation
perfectly matches the exact interval when `cc` is TRUE (i.e. for
Clopper-Pearson) but not when `cc` is `FALSE` (for mid-p) All methods
except Agresti-Coull have equivalent formulae for the Poisson
distribution: Garwood for Clopper-Pearson, Rao score for Wilson score.
Jeffreys has a Poisson equivalent using the Gamma distribution. e.g. See
Brown et al. 2003, Swift 2009 and Laud 2017. The formulation for the
approximate mid-p interval using Gamma distribution for a Poisson rate
has been deduced by the package author from the corresponding formulae
from Brown et al., and has not (to the best of my knowledge) been
published.

This function is vectorised in x, n.

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
