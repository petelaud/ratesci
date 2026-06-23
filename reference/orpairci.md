# Confidence intervals for conditional odds ratio (OR) with paired binomial rates.

Confidence intervals for comparisons of two binomial rates from paired
data. This convenience wrapper function produces a selection of the
methods below for the conditional odds ratio (OR) contrast, with or
without optional continuity adjustment (where available).

- Transformed SCAS (skewness-corrected asymptotic score)

- Transformed Wilson Score method

- Transformed mid-P

- Transformed Jeffreys

- Approximate log-normal (Wald) method

## Usage

``` r
orpairci(x, level = 0.95, std_est = TRUE, cc = FALSE, precis = 8)
```

## Arguments

- x:

  A numeric vector object specified as c(a, b, c, d) where:  
  a is the number of pairs with the event (e.g. success) under both
  conditions (e.g. treated/untreated, or case/control)  
  b is the count of the number with the event on condition 1 only (=
  x12)  
  c is the count of the number with the event on condition 2 only (=
  x21)  
  d is the number of pairs with no event under both conditions  
  (Note the order of a and d is only important for contrast="RR".)

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
  places) to be used in output.

## Value

A list containing the following components:

- data:

  the input data in 2x2 matrix form.

- estimates:

  an array containing the confidence interval for paired OR using
  various methods. The methods shown depends on the cc argument (if cc =
  TRUE then the continuity-adjusted methods are given).

- call:

  details of the function call.

## References

Fagerland MW, Lydersen S, Laake P. Recommended tests and confidence
intervals for paired binomial proportions. Statistics in Medicine 2014;
33(16):2850-2875

Laud PJ. Improved confidence intervals and tests for paired binomial
proportions. (2026, Under review)

## Author

Pete Laud, <p.j.laud@sheffield.ac.uk>

## Examples

``` r
# Example data from Fagerland et al 2014
orpairci(x = c(1, 1, 7, 12), precis = 3)
#> [[1]]
#>          Test_2
#> Test_1    Success Failure
#>   Success       1       1
#>   Failure       7      12
#> 
#> $estimates
#>                      lower   est upper
#> Transformed SCASp    0.008 0.143 0.912
#> Transformed midp     0.006 0.143 0.924
#> Transformed Wilson   0.023 0.143 0.890
#> Transformed Jeffreys 0.014 0.143 0.831
#> Wald                 0.018 0.143 1.161
#> 
#> $call
#> level    cc 
#>  0.95  0.00 
#> 
# with conventional continuity adjustment
orpairci(x = c(1, 1, 7, 12), precis = 3, cc = TRUE)
#> [[1]]
#>          Test_2
#> Test_1    Success Failure
#>   Success       1       1
#>   Failure       7      12
#> 
#> $estimates
#>                             lower   est upper
#> Transformed SCASp_cc        0.000 0.143 1.204
#> Transformed Clopper-Pearson 0.003 0.143 1.112
#> Transformed Wilson_cc       0.007 0.143 1.142
#> Transformed Jeffreys_cc     0.003 0.143 1.112
#> 
#> $call
#> level    cc 
#>  0.95  0.50 
#> 
# with intermediate continuity adjustment
orpairci(x = c(1, 1, 7, 12), precis = 3, cc = 0.25)
#> [[1]]
#>          Test_2
#> Test_1    Success Failure
#>   Success       1       1
#>   Failure       7      12
#> 
#> $estimates
#>                               lower   est upper
#> Transformed SCASp_cc(0.25)    0.002 0.143 1.052
#> Transformed midp_cc(0.25)     0.004 0.143 1.029
#> Transformed Wilson_cc(0.25)   0.014 0.143 1.009
#> Transformed Jeffreys_cc(0.25) 0.008 0.143 0.966
#> 
#> $call
#> level    cc 
#>  0.95  0.25 
#> 
```
