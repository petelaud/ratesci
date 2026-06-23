# Selected confidence intervals for the single binomial or Poisson rate.

Confidence intervals for the single binomial or Poisson rate. This
convenience wrapper function produces a selection of alternative
methods. The first three are recommended for achieving 1-sided and
2-sided coverage probability close to the nominal levels (see Laud 2017
and Laud 2018):

- SCAS (skewness-corrected asymptotic score)

- Jeffreys

- midp (two versions, using exact calculation or approximation via
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
`TRUE` or `0.5`, the mid-p method becomes the 'exact' Clopper-Pearson
interval (or Garwood for Poisson rates), and the output also includes
the slightly less conservative Blaker interval. Hence the full list of
conservative methods produced (when `cc` is `TRUE`) is:

- SCAS_cc

- Jeffreys_cc

- Clopper-Pearson 'exact' (two identical versions, using exact
  calculation or approximation via Beta/Gamma distribution))

- Wilson_cc

- Wald_cc (strongly advise this is not used for any purpose but included
  for reference)

- Blaker 'exact'

Note that Brown et al's Beta formulation perfectly matches the exact
interval when `cc` is TRUE (i.e. for Clopper-Pearson) but not when `cc`
is `FALSE` (for mid-p) All methods except Agresti-Coull have equivalent
formulae for the Poisson distribution: Garwood for Clopper-Pearson, Rao
score for Wilson score. Jeffreys has a Poisson equivalent using the
Gamma distribution. e.g. See Brown et al. 2003, Swift 2009 and Laud
2017. The formulation for the approximate mid-p interval using Gamma
distribution for a Poisson rate has been deduced by the package author
from the corresponding formulae from Brown et al., and has not (to the
best of my knowledge) been published.

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
continuity adjustment is applied. The corresponding 'exact' method is
Clopper-Pearson/Garwood if cc = TRUE and mid-p if cc = FALSE. An
additional output object 'estimates' is provided for a side-by-side
comparison of all methods. These are also grouped depending on the cc
argument (if cc = TRUE then the continuity-adjusted and exact strictly
conservative methods are included) The last list item contains details
of the function call.

## References

Laud PJ. Equal-tailed confidence intervals for comparison of rates.
Pharmaceutical Statistics 2017; 16:334-348. (Appendix A.4)

Brown LD, Cai TT and DasGupta A. Interval estimation for a binomial
proportion. Statistical Science 2001; 16(2):101-133.

Garwood F. Fiducial limits for the Poisson distribution. Biometrika
1936; 28(3-4):437, doi:10.1093/biomet/28.3-4.437.

## Author

Pete Laud, <p.j.laud@sheffield.ac.uk>

## Examples

``` r
# Selected example datasets from Newcombe 1998
rateci(
  x = c(15, 1), n = c(148, 29),
  precis = 4
)
#> $scas
#>            lower        est     upper  x   n
#> [1,] 0.060212460 0.10135135 0.1579114 15 148
#> [2,] 0.001991554 0.03448276 0.1548909  1  29
#> 
#> $jeff
#>            lower        est     upper  x   n
#> [1,] 0.060448517 0.10135135 0.1576431 15 148
#> [2,] 0.003746174 0.03448276 0.1500777  1  29
#> 
#> $midp
#>            lower        est     upper  x   n
#> [1,] 0.060062408 0.10135135 0.1581230 15 148
#> [2,] 0.001724243 0.03448276 0.1585388  1  29
#> 
#> $midp_beta
#>            lower        est     upper  x   n
#> [1,] 0.060456514 0.10135135 0.1576351 15 148
#> [2,] 0.004668305 0.03448276 0.1485446  1  29
#> 
#> $wilson
#>            lower        est     upper  x   n
#> [1,] 0.062386400 0.10135135 0.1604872 15 148
#> [2,] 0.006113214 0.03448276 0.1717552  1  29
#> 
#> $wald
#>                                        x   n
#> [1,]  0.05273001 0.10135135 0.1499727 15 148
#> [2,] -0.03192673 0.03448276 0.1008922  1  29
#> 
#> $estimates
#> , , 1
#> 
#>                lower    est  upper  x   n
#> SCAS          0.0602 0.1014 0.1579 15 148
#> Jeffreys      0.0604 0.1014 0.1576 15 148
#> midp          0.0601 0.1014 0.1581 15 148
#> midp(beta)    0.0605 0.1014 0.1576 15 148
#> Wilson        0.0624 0.1014 0.1605 15 148
#> Wald          0.0527 0.1014 0.1500 15 148
#> Agresti-Coull 0.0614 0.1014 0.1615 15 148
#> 
#> , , 2
#> 
#>                 lower    est  upper x  n
#> SCAS           0.0020 0.0345 0.1549 1 29
#> Jeffreys       0.0037 0.0345 0.1501 1 29
#> midp           0.0017 0.0345 0.1585 1 29
#> midp(beta)     0.0047 0.0345 0.1485 1 29
#> Wilson         0.0061 0.0345 0.1718 1 29
#> Wald          -0.0319 0.0345 0.1009 1 29
#> Agresti-Coull -0.0084 0.0345 0.1863 1 29
#> 
#> 
#> $call
#> distrib   level      cc std_est 
#>   "bin"  "0.95" "FALSE"  "TRUE" 
#> 
# With conventional continuity adjustment
rateci(
  x = c(15, 1), n = c(148, 29),
  precis = 4, cc = TRUE
)
#> $scas_cc
#>             lower        est     upper  x   n
#> [1,] 5.760455e-02 0.10135135 0.1619137 15 148
#> [2,] 6.185396e-06 0.03448276 0.1816531  1  29
#> 
#> $jeff_cc
#>             lower        est     upper  x   n
#> [1,] 0.0578440101 0.10135135 0.1616505 15 148
#> [2,] 0.0008726469 0.03448276 0.1776443  1  29
#> 
#> $cp
#>            lower        est     upper  x   n
#> [1,] 0.057842255 0.10135135 0.1616516 15 148
#> [2,] 0.000869751 0.03448276 0.1776466  1  29
#> 
#> $cp_beta
#>             lower        est     upper  x   n
#> [1,] 0.0578440101 0.10135135 0.1616505 15 148
#> [2,] 0.0008726469 0.03448276 0.1776443  1  29
#> 
#> $blaker
#>            lower        est     upper  x   n
#> [1,] 0.057944010 0.10135135 0.1601505 15 148
#> [2,] 0.001722647 0.03448276 0.1660443  1  29
#> 
#> $estimates
#> , , 1
#> 
#>                        lower    est  upper  x   n
#> SCAS_cc               0.0576 0.1014 0.1619 15 148
#> Jeffreys_cc           0.0578 0.1014 0.1617 15 148
#> Clopper-Pearson       0.0578 0.1014 0.1617 15 148
#> Clopper-Pearson(beta) 0.0578 0.1014 0.1617 15 148
#> Wilson_cc             0.0598 0.1014 0.1644 15 148
#> Wald_cc               0.0494 0.1014 0.1534 15 148
#> Blaker                0.0579 0.1014 0.1602 15 148
#> 
#> , , 2
#> 
#>                         lower    est  upper x  n
#> SCAS_cc                0.0000 0.0345 0.1817 1 29
#> Jeffreys_cc            0.0009 0.0345 0.1776 1 29
#> Clopper-Pearson        0.0009 0.0345 0.1776 1 29
#> Clopper-Pearson(beta)  0.0009 0.0345 0.1776 1 29
#> Wilson_cc              0.0018 0.0345 0.1963 1 29
#> Wald_cc               -0.0492 0.0345 0.1181 1 29
#> Blaker                 0.0017 0.0345 0.1660 1 29
#> 
#> 
#> $call
#> distrib   level      cc std_est 
#>   "bin"  "0.95"   "0.5"  "TRUE" 
#> 
# With intermediate continuity adjustment
rateci(
  x = c(15, 1), n = c(148, 29),
  precis = 4, cc = 0.25
)
#> $scas_cc
#>             lower        est     upper  x   n
#> [1,] 0.0589064809 0.10135135 0.1599145 15 148
#> [2,] 0.0006046622 0.03448276 0.1685137  1  29
#> 
#> $jeff_cc
#>            lower        est     upper  x   n
#> [1,] 0.059144221 0.10135135 0.1596488 15 148
#> [2,] 0.002052028 0.03448276 0.1641487  1  29
#> 
#> $beta_cc
#>            lower        est     upper  x   n
#> [1,] 0.059150262 0.10135135 0.1596428 15 148
#> [2,] 0.002770476 0.03448276 0.1630944  1  29
#> 
#> $estimates
#> , , 1
#> 
#>                      lower    est  upper  x   n
#> SCAS_cc(0.25)       0.0589 0.1014 0.1599 15 148
#> Jeffreys_cc(0.25)   0.0591 0.1014 0.1596 15 148
#> midp_cc(0.25)       0.0589 0.1014 0.1600 15 148
#> midp(beta)_cc(0.25) 0.0592 0.1014 0.1596 15 148
#> Wilson_cc(0.25)     0.0611 0.1014 0.1625 15 148
#> Wald_cc(0.25)       0.0510 0.1014 0.1517 15 148
#> 
#> , , 2
#> 
#>                       lower    est  upper x  n
#> SCAS_cc(0.25)        0.0006 0.0345 0.1685 1 29
#> Jeffreys_cc(0.25)    0.0021 0.0345 0.1641 1 29
#> midp_cc(0.25)        0.0012 0.0345 0.1694 1 29
#> midp(beta)_cc(0.25)  0.0028 0.0345 0.1631 1 29
#> Wilson_cc(0.25)      0.0038 0.0345 0.1841 1 29
#> Wald_cc(0.25)       -0.0405 0.0345 0.1095 1 29
#> 
#> 
#> $call
#> distrib   level      cc std_est 
#>   "bin"  "0.95"  "0.25"  "TRUE" 
#> 
```
