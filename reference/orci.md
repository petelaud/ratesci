# Confidence intervals for odds ratio with independent binomial proportions.

Confidence intervals for comparisons of two binomial rates. This
convenience wrapper function produces a selection of the methods below
as appropriate for the odds ratio contrast (OR), with or without
continuity adjustment.

- SCAS (skewness-corrected asymptotic score)

- Miettinen-Nurminen, Gart and unadjusted Asymptotic Score methods

- MOVER-W

- MOVER-J (based on Jeffreys intervalse)

- Approximate normal (Woolf/Gart) methods (strongly advise this is not
  used for any purpose but included for reference)

## Usage

``` r
orci(x1, n1, x2, n2, level = 0.95, std_est = TRUE, cc = FALSE, precis = 6)
```

## Arguments

- x1, x2:

  Numeric vectors of numbers of events in group 1 & group 2
  respectively.

- n1, n2:

  Numeric vectors of sample sizes (for binomial rates) or exposure times
  (for Poisson rates) in each group.

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

  Number (default 6) specifying precision (i.e. number of decimal
  places) to be used in root-finding subroutine for the score confidence
  intervals. (Note other methods use closed-form calculations so are not
  affected.)

## Value

A list containing the following components:

- estimates:

  an array containing the confidence interval for RR using various
  methods. The methods shown depends on the cc argument (if cc = TRUE
  then the continuity-adjusted methods are given).

- call:

  details of the function call.

## References

Laud PJ. Equal-tailed confidence intervals for comparison of rates.
Pharmaceutical Statistics 2017; 16:334-348.

## Author

Pete Laud, <pete@sheffstat.co.uk>

## Examples

``` r
# Selected example datasets from Newcombe 1998 and Fagerland et al. 2011
# (note Fagerland et al. appear to have the Miettinen-Nurminen method
#  labelled as Koopman)
orci(
  x1 = c(5, 7), n1 = c(56, 34),
  x2 = c(0, 1), n2 = c(29, 34),
  precis = 4
)
#> $estimates
#> , , 5/56 vs 0/29
#> 
#>                               lower est    upper
#> SCAS                         0.7587 Inf      Inf
#> Gart                         0.7720 Inf      Inf
#> Miettinen-Nurminen           0.6959 Inf      Inf
#> Uncorrected Asymptotic Score 0.7046 Inf      Inf
#> MOVER-W                      0.5982 Inf      Inf
#> MOVER-J                      0.8639 Inf      Inf
#> Woolf logit                  0.0000 Inf      Inf
#> Gart adjusted logit          0.3364 Inf 118.0279
#> 
#> , , 7/34 vs 1/34
#> 
#>                               lower    est    upper
#> SCAS                         1.2714 8.5556 169.0681
#> Gart                         1.2883 8.5556 171.0103
#> Miettinen-Nurminen           1.2466 8.5556  56.4616
#> Uncorrected Asymptotic Score 1.2613 8.5556  55.8483
#> MOVER-W                      1.1799 8.5556  56.4428
#> MOVER-J                      1.3273 8.5556  87.4833
#> Woolf logit                  0.9905 8.5556  73.9003
#> Gart adjusted logit          0.9828 8.5556  37.7486
#> 
#> 
#> $call
#> distrib   level      cc 
#>   "bin"  "0.95" "FALSE" 
#> 
# With conventional continuity adjustment
orci(
  x1 = c(5, 7), n1 = c(56, 34),
  x2 = c(0, 1), n2 = c(29, 34),
  precis = 4, cc = TRUE
)
#> $estimates
#> , , 5/56 vs 0/29
#> 
#>                                  lower est upper
#> SCAS_cc                         0.4399 Inf   Inf
#> Gart_cc                         0.4452 Inf   Inf
#> Miettinen-Nurminen_cc           0.4440 Inf   Inf
#> Uncorrected Asymptotic Score_cc 0.4486 Inf   Inf
#> MOVER-W_cc                      0.4380 Inf   Inf
#> MOVER-J_cc                      0.5336 Inf   Inf
#> 
#> , , 7/34 vs 1/34
#> 
#>                                  lower    est      upper
#> SCAS_cc                         0.9280 8.5556 37072.9285
#> Gart_cc                         0.9388 8.5556 58169.6654
#> Miettinen-Nurminen_cc           0.9332 8.5556   199.1816
#> Uncorrected Asymptotic Score_cc 0.9433 8.5556   196.7604
#> MOVER-W_cc                      0.9633 8.5556   176.8584
#> MOVER-J_cc                      1.0428 8.5556   360.6648
#> 
#> 
#> $call
#> distrib   level      cc 
#>   "bin"  "0.95"   "0.5" 
#> 
# With intermediate continuity adjustment
orci(
  x1 = c(5, 7), n1 = c(56, 34),
  x2 = c(0, 1), n2 = c(29, 34),
  precis = 4, cc = 0.25
)
#> $estimates
#> , , 5/56 vs 0/29
#> 
#>                                        lower est upper
#> SCAS_cc(0.25)                         0.5664 Inf   Inf
#> Gart_cc(0.25)                         0.5743 Inf   Inf
#> Miettinen-Nurminen_cc(0.25)           0.5498 Inf   Inf
#> Uncorrected Asymptotic Score_cc(0.25) 0.5560 Inf   Inf
#> MOVER-W_cc(0.25)                      0.5080 Inf   Inf
#> MOVER-J_cc(0.25)                      0.6592 Inf   Inf
#> 
#> , , 7/34 vs 1/34
#> 
#>                                        lower    est    upper
#> SCAS_cc(0.25)                         1.0820 8.5556 560.3375
#> Gart_cc(0.25)                         1.0954 8.5556 579.6362
#> Miettinen-Nurminen_cc(0.25)           1.0754 8.5556  94.9040
#> Uncorrected Asymptotic Score_cc(0.25) 1.0875 8.5556  93.8207
#> MOVER-W_cc(0.25)                      1.0638 8.5556  88.2634
#> MOVER-J_cc(0.25)                      1.1715 8.5556 156.2726
#> 
#> 
#> $call
#> distrib   level      cc 
#>   "bin"  "0.95"  "0.25" 
#> 
```
