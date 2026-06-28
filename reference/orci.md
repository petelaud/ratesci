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
#> SCAS                         0.7696 Inf      Inf
#> Gart                         0.7822 Inf      Inf
#> Miettinen-Nurminen           0.7175 Inf      Inf
#> Uncorrected Asymptotic Score 0.7257 Inf      Inf
#> MOVER-W                      0.6292 Inf      Inf
#> MOVER-J                      0.8753 Inf      Inf
#> Woolf logit                  0.0000 Inf      Inf
#> Gart adjusted logit          0.3312 Inf 101.2048
#> 
#> , , 7/34 vs 1/34
#> 
#>                               lower    est    upper
#> SCAS                         1.2383 8.5556 129.3294
#> Gart                         1.2532 8.5556 130.8308
#> Miettinen-Nurminen           1.2086 8.5556  43.0330
#> Uncorrected Asymptotic Score 1.2209 8.5556  42.5757
#> MOVER-W                      1.1537 8.5556  41.9763
#> MOVER-J                      1.2771 8.5556  67.1689
#> Woolf logit                  0.9096 8.5556  53.8695
#> Gart adjusted logit          0.9233 8.5556  27.0778
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
#> SCAS_cc                         0.4628 Inf   Inf
#> Gart_cc                         0.4678 Inf   Inf
#> Miettinen-Nurminen_cc           0.4765 Inf   Inf
#> Uncorrected Asymptotic Score_cc 0.4810 Inf   Inf
#> MOVER-W_cc                      0.4779 Inf   Inf
#> MOVER-J_cc                      0.5654 Inf   Inf
#> 
#> , , 7/34 vs 1/34
#> 
#>                                  lower    est      upper
#> SCAS_cc                         0.9360 8.5556 27593.3534
#> Gart_cc                         0.9456 8.5556 43279.8893
#> Miettinen-Nurminen_cc           0.9424 8.5556   149.5435
#> Uncorrected Asymptotic Score_cc 0.9512 8.5556   147.7432
#> MOVER-W_cc                      0.9712 8.5556   136.6421
#> MOVER-J_cc                      1.0356 8.5556   282.5531
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
#> SCAS_cc(0.25)                         0.5853 Inf   Inf
#> Gart_cc(0.25)                         0.5928 Inf   Inf
#> Miettinen-Nurminen_cc(0.25)           0.5788 Inf   Inf
#> Uncorrected Asymptotic Score_cc(0.25) 0.5847 Inf   Inf
#> MOVER-W_cc(0.25)                      0.5444 Inf   Inf
#> MOVER-J_cc(0.25)                      0.6838 Inf   Inf
#> 
#> , , 7/34 vs 1/34
#> 
#>                                        lower    est    upper
#> SCAS_cc(0.25)                         1.0724 8.5556 422.9320
#> Gart_cc(0.25)                         1.0843 8.5556 437.4584
#> Miettinen-Nurminen_cc(0.25)           1.0643 8.5556  71.7594
#> Uncorrected Asymptotic Score_cc(0.25) 1.0746 8.5556  70.9530
#> MOVER-W_cc(0.25)                      1.0563 8.5556  66.8112
#> MOVER-J_cc(0.25)                      1.1453 8.5556 121.2444
#> 
#> 
#> $call
#> distrib   level      cc 
#>   "bin"  "0.95"  "0.25" 
#> 
```
