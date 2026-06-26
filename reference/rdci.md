# Confidence intervals for rate difference (RD) with independent binomial or Poisson rates.

Confidence intervals for comparisons of two binomial or Poisson rates.
This convenience wrapper function produces a selection of the methods
below as appropriate for the selected distribution (binomial or Poisson)
for the rate difference (RD) contrast, with or without continuity
adjustment.

- SCAS (skewness-corrected asymptotic score)

- Miettinen-Nurminen, Mee, Koopman, Gart-Nam Asymptotic Score methods

- MOVER-W (Method of Variance Estimates Recovery based on Wilson
  intervals, aka Newcombe Hybrid Score or 'square-and-add')

- MOVER-J (based on Jeffreys intervals)

- Agresti-Caffo (binomial RD only)

- Approximate normal (Wald) method (strongly advise this is not used for
  any purpose but included for reference)

## Usage

``` r
rdci(
  x1,
  n1,
  x2,
  n2,
  distrib = "bin",
  level = 0.95,
  std_est = TRUE,
  cc = FALSE,
  precis = 8
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

- level:

  Number specifying confidence level (between 0 and 1, default 0.95).

- std_est:

  logical, specifying if the crude point estimate for the contrast value
  should be returned (TRUE, default) or the method-specific alternative
  point estimate consistent with a 0% confidence interval (FALSE).

- cc:

  Number or logical (default FALSE) specifying (amount of) continuity
  adjustment. Numeric value between 0 and 0.5 is taken as the gamma
  parameter in Laud 2017, Appendix S2 (`cc = TRUE` translates to 0.5 for
  'conventional' Yates adjustment).  

- precis:

  Number (default 8) specifying precision (i.e. number of decimal
  places) to be used in root-finding subroutine for the score confidence
  intervals. (Note other methods use closed-form calculations so are not
  affected.)

## Value

A list containing the following components:

- estimates:

  an array containing the confidence interval for RD using various
  methods. The methods shown depends on the cc argument (if cc = TRUE
  then the continuity-adjusted methods are given).

- call:

  details of the function call.

## References

Laud PJ. Equal-tailed confidence intervals for comparison of rates.
Pharmaceutical Statistics 2017; 16:334-348.

Newcombe RG. Interval estimation for the difference between independent
proportions: comparison of eleven methods. Statistics in Medicine 1998;
17(8):873-890.

## Author

Pete Laud, <pete@sheffstat.co.uk>

## Examples

``` r
# Selected example datasets from Newcombe 1998 and Fagerland et al. 2011
# (note Fagerland et al. have the Mee method labelled as Miettinen-Nurminen)
rdci(
  x1 = c(5, 7), n1 = c(56, 34),
  x2 = c(0, 1), n2 = c(29, 34),
  precis = 4
)
#> $estimates
#> , , 5/56 vs 0/29
#> 
#>                      lower    est  upper
#> SCAS               -0.0186 0.0893 0.1867
#> Gart-Nam           -0.0172 0.0893 0.1859
#> Miettinen-Nurminen -0.0326 0.0893 0.1933
#> Mee                -0.0313 0.0893 0.1926
#> MOVER-W            -0.0381 0.0893 0.1926
#> MOVER-J            -0.0083 0.0893 0.1847
#> Wald                0.0146 0.0893 0.1640
#> Agresti-Caffo      -0.0289 0.0893 0.1712
#> 
#> , , 7/34 vs 1/34
#> 
#>                     lower    est  upper
#> SCAS               0.0261 0.1765 0.3424
#> Gart-Nam           0.0274 0.1765 0.3410
#> Miettinen-Nurminen 0.0270 0.1765 0.3453
#> Mee                0.0284 0.1765 0.3439
#> MOVER-W            0.0189 0.1765 0.3404
#> MOVER-J            0.0278 0.1765 0.3310
#> Wald               0.0292 0.1765 0.3238
#> Agresti-Caffo      0.0116 0.1765 0.3217
#> 
#> 
#> $call
#> distrib   level      cc 
#>   "bin"  "0.95" "FALSE" 
#> 
# With conventional continuity adjustment
rdci(
  x1 = c(5, 7), n1 = c(56, 34),
  x2 = c(0, 1), n2 = c(29, 34),
  precis = 4, cc = TRUE
)
#> $estimates
#> , , 5/56 vs 0/29
#> 
#>                         lower    est  upper
#> SCAS_cc               -0.0482 0.0893 0.2087
#> Gart-Nam_cc           -0.0468 0.0893 0.2079
#> Miettinen-Nurminen_cc -0.0613 0.0893 0.2147
#> Mee_cc                -0.0601 0.0893 0.2139
#> MOVER-W_cc            -0.0667 0.0893 0.2037
#> MOVER-J_cc            -0.0429 0.0893 0.1962
#> Wald_cc               -0.0116 0.0893 0.1901
#> Hauck-Anderson        -0.0033 0.0893 0.1819
#> 
#> , , 7/34 vs 1/34
#> 
#>                         lower    est  upper
#> SCAS_cc                0.0089 0.1765 0.3587
#> Gart-Nam_cc            0.0102 0.1765 0.3573
#> Miettinen-Nurminen_cc  0.0091 0.1765 0.3612
#> Mee_cc                 0.0105 0.1765 0.3598
#> MOVER-W_cc            -0.0040 0.1765 0.3568
#> MOVER-J_cc             0.0042 0.1765 0.3478
#> Wald_cc               -0.0002 0.1765 0.3532
#> Hauck-Anderson         0.0122 0.1765 0.3407
#> 
#> 
#> $call
#> distrib   level      cc 
#>   "bin"  "0.95"   "0.5" 
#> 
# With intermediate continuity adjustment
rdci(
  x1 = c(5, 7), n1 = c(56, 34),
  x2 = c(0, 1), n2 = c(29, 34),
  precis = 4, cc = 0.25
)
#> $estimates
#> , , 5/56 vs 0/29
#> 
#>                               lower    est  upper
#> SCAS_cc(0.25)               -0.0338 0.0893 0.1978
#> Gart-Nam_cc(0.25)           -0.0324 0.0893 0.1969
#> Miettinen-Nurminen_cc(0.25) -0.0474 0.0893 0.2041
#> Mee_cc(0.25)                -0.0461 0.0893 0.2033
#> MOVER-W_cc(0.25)            -0.0527 0.0893 0.1981
#> MOVER-J_cc(0.25)            -0.0264 0.0893 0.1904
#> Wald_cc(0.25)                0.0015 0.0893 0.1771
#> 
#> , , 7/34 vs 1/34
#> 
#>                              lower    est  upper
#> SCAS_cc(0.25)               0.0175 0.1765 0.3506
#> Gart-Nam_cc(0.25)           0.0188 0.1765 0.3492
#> Miettinen-Nurminen_cc(0.25) 0.0181 0.1765 0.3532
#> Mee_cc(0.25)                0.0195 0.1765 0.3519
#> MOVER-W_cc(0.25)            0.0074 0.1765 0.3486
#> MOVER-J_cc(0.25)            0.0159 0.1765 0.3395
#> Wald_cc(0.25)               0.0145 0.1765 0.3385
#> 
#> 
#> $call
#> distrib   level      cc 
#>   "bin"  "0.95"  "0.25" 
#> 
```
