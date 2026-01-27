# Method of Variance Estimates Recovery ("MOVER") confidence intervals for comparisons of independent binomial or Poisson rates.

Confidence intervals applying the MOVER method ("Method of Variance
Estimates Recovery", developed from the Newcombe method for binomial RD)
across different contrasts (RD, RR, OR) and distributions (binomial,
Poisson) using equal-tailed Jeffreys intervals instead of the Wilson
score method for the event rates. Also allows more general Beta and
Gamma priors for an approximate Bayesian confidence interval
incorporating prior beliefs about the group event rates. This function
is vectorised in x1, x2, n1, and n2.

## Usage

``` r
moverci(
  x1,
  n1,
  x2 = NULL,
  n2 = NULL,
  distrib = "bin",
  contrast = "RD",
  level = 0.95,
  a1 = 0.5,
  b1 = 0.5,
  a2 = 0.5,
  b2 = 0.5,
  type = "jeff",
  adj = FALSE,
  cc = FALSE,
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
  "bin" = binomial (default);  
  "poi" = Poisson.

- contrast:

  Character string indicating the contrast of interest:  
  "RD" = rate difference (default);  
  "RR" = rate ratio;  
  "OR" = odds ratio;  
  "p" gives an interval for the single proportion `x1/n1`.

- level:

  Number specifying confidence level (between 0 and 1, default 0.95).

- a1, b1, a2, b2:

  Numbers defining the Beta(ai,bi) prior distributions for each group
  (default ai = bi = 0.5 for Jeffreys method). Gamma priors for Poisson
  rates require only a1, a2.

- type:

  Character string indicating the method used for the intervals for the
  individual group rates.  
  "jeff" = Jeffreys equal-tailed intervals (default);  
  "exact" = Clopper-Pearson/Garwood exact intervals (note this does NOT
  result in a strictly conservative interval for the contrast, except
  for contrast = "p". The scoreci function with `cc = TRUE` is
  recommended as a superior approximation of 'exact' methods);  
  "midp" = mid-p intervals;  
  "SCAS" = SCAS non-iterative intervals;  
  "wilson" = Wilson score intervals (as per Newcombe 1998). (Rao score
  is used for `distrib = "poi"`)  
  NB: "wilson" option is included only for legacy validation against
  previous published method by Newcombe. It is not recommended, as
  `type = "jeff"` or other equal-tailed options achieve much better
  coverage properties.

- adj:

  Logical (default FALSE) indicating whether to apply the boundary
  adjustment for Jeffreys intervals recommended on p108 of Brown et al.
  (`type = "jeff"` only: set to FALSE if using informative priors.)

- cc:

  Number or logical specifying (amount of) continuity adjustment
  (default FALSE). Numeric value is taken as the gamma parameter in Laud
  2017, Appendix S2 (default 0.5 if `cc = TRUE`). Forced equal to 0.5 if
  `type = "exact"`.

- ...:

  Additional arguments.

## Value

A list containing the following components:

- estimates:

  a matrix containing estimates of the rates in each group and of the
  requested contrast, with its confidence interval.

- call:

  details of the function call.

## References

Laud PJ. Equal-tailed confidence intervals for comparison of rates.
Pharmaceutical Statistics 2017; 16:334-348.

Newcombe RG. Interval estimation for the difference between independent
proportions: comparison of eleven methods. Statistics in Medicine 1998;
17(8):873-890.

Donner A, Zou G. Closed-form confidence intervals for functions of the
normal mean and standard deviation. Statistical Methods in Medical
Research 2012; 21(4):347-359.

Fagerland MW, Newcombe RG. Confidence intervals for odds ratio and
relative risk based on the inverse hyperbolic sine transformation.
Statistics in Medicine 2013; 32(16):2823-2836.

Li HQ, Tang ML, Wong WK. Confidence intervals for ratio of two Poisson
rates using the method of variance estimates recovery. Computational
Statistics 2014; 29(3-4):869-889.

## Author

Pete Laud, <p.j.laud@sheffield.ac.uk>

## Examples

``` r
# Binomial RD, MOVER-J method:
moverci(x1 = 5, n1 = 56, x2 = 0, n2 = 29)
#> $estimates
#>             lower        est     upper level x1 n1 x2 n2      p1hat       p2hat
#> [1,] -0.009734968 0.08403347 0.1772258  0.95  5 56  0 29 0.09177968 0.007746206
#> 
#> $call
#>  distrib contrast    level     type      adj       cc       a1       b1 
#>    "bin"     "RD"   "0.95"   "jeff"  "FALSE"  "FALSE"    "0.5"    "0.5" 
#>       a2       b2 
#>    "0.5"    "0.5" 
#> 

# Binomial RD, Newcombe method:
moverci(x1 = 5, n1 = 56, x2 = 0, n2 = 29, type = "wilson")
#> $estimates
#>            lower        est   upper level x1 n1 x2 n2      p1hat p2hat
#> [1,] -0.03813715 0.08928571 0.19256  0.95  5 56  0 29 0.08928571     0
#> 
#> $call
#>  distrib contrast    level     type      adj       cc       a1       b1 
#>    "bin"     "RD"   "0.95" "wilson"  "FALSE"  "FALSE"    "0.5"    "0.5" 
#>       a2       b2 
#>    "0.5"    "0.5" 
#> 
```
