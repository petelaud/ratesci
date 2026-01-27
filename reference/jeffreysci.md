# Jeffreys and other approximate Bayesian confidence intervals for a single binomial or Poisson rate.

Generalised approximate Bayesian confidence intervals based on a Beta
(for binomial rates) or Gamma (for Poisson rates) conjugate priors.
Encompassing the Jeffreys method (with Beta(0.5, 0.5) or Gamma(0.5)
respectively), as well as any user-specified prior distribution.
Clopper-Pearson method (as quantiles of a Beta distribution as described
in Brown et al. 2001) also included by way of a "continuity adjustment"
parameter.

## Usage

``` r
jeffreysci(
  x,
  n,
  ai = 0.5,
  bi = 0.5,
  cc = 0,
  level = 0.95,
  distrib = "bin",
  adj = TRUE,
  ...
)
```

## Arguments

- x:

  Numeric vector of number of events.

- n:

  Numeric vector of sample sizes (for binomial rates) or exposure times
  (for Poisson rates).

- ai, bi:

  Numbers defining the Beta prior distribution (default \`ai = bi = 0.5“
  for Jeffreys interval). Gamma prior for Poisson rates requires only
  ai.

- cc:

  Number or logical specifying (amount of) "continuity adjustment". cc =
  0 (default) gives Jeffreys interval, `cc = 0.5` gives the
  Clopper-Pearson interval (or Garwood for Poisson). A value between 0
  and 0.5 allows a compromise between proximate and conservative
  coverage.

- level:

  Number specifying confidence level (between 0 and 1, default 0.95).

- distrib:

  Character string indicating distribution assumed for the input data:  
  "bin" = binomial (default);  
  "poi" = Poisson.

- adj:

  Logical (default TRUE) indicating whether to apply the boundary
  adjustment recommended on p108 of Brown et al. (set to FALSE if
  informative priors are used).

- ...:

  Other arguments.

## Value

A list containing the following components:

- estimates:

  a matrix containing estimated rate(s), and corresponding approximate
  Bayesian confidence interval, and the input values x and n.

- call:

  details of the function call.

## References

Laud PJ. Equal-tailed confidence intervals for comparison of rates.
Pharmaceutical Statistics 2017; 16:334-348.

Brown LD, Cai TT, DasGupta A. Interval estimation for a binomial
proportion. Statistical Science 2001; 16(2):101-133

## Author

Pete Laud, <p.j.laud@sheffield.ac.uk>

## Examples

``` r
# Jeffreys method:
jeffreysci(x = 5, n = 56)
#> $estimates
#>           lower        est     upper x  n
#> [1,] 0.03489147 0.09177968 0.1846509 5 56
#> 
#> $call
#> distrib   level      cc     adj      ai      bi 
#>   "bin"  "0.95"     "0"  "TRUE"   "0.5"   "0.5" 
#> 
```
