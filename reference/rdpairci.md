# Confidence intervals for rate difference (RD) with paired binomial rates.

Confidence intervals for comparisons of two binomial rates from paired
data. This convenience wrapper function produces a selection of the
methods below for the rate difference (RD) contrast, with or without
optional continuity adjustment (where available).

- SCAS (skewness-corrected asymptotic score)

- SCASu (omitting the 'N-1' adjustment)

- Tango Asymptotic Score method

- MOVER-NW (aka Newcombe Hybrid Score or square-and-add)

- MOVER-NJ (based on Jeffreys intervals)

- Agresti-Min

- Bonett-Price

- Approximate normal (Wald) method (strongly advise this is not used for
  any purpose but included for reference)

## Usage

``` r
rdpairci(x, level = 0.95, std_est = TRUE, cc = FALSE, precis = 8)
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

  Number (default 8) specifying precision (i.e. number of decimal
  places) to be used in root-finding subroutine for the score confidence
  intervals. (Note other methods use closed-form calculations so are not
  affected.)

## Value

A list containing the following components:

- data:

  the input data in 2x2 matrix form.

- estimates:

  an array containing the confidence interval for paired RD using
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

Pete Laud, <pete@sheffstat.co.uk>

## Examples

``` r
# Example data from Fagerland et al 2014
rdpairci(x = c(1, 1, 7, 12), precis = 3)
#> [[1]]
#>          Test_2
#> Test_1    Success Failure
#>   Success       1       1
#>   Failure       7      12
#> 
#> $estimates
#>               lower    est  upper
#> SCAS         -0.528 -0.286 -0.018
#> SCASu        -0.523 -0.286 -0.026
#> Tango score  -0.517 -0.286 -0.026
#> MOVER-NW     -0.507 -0.286 -0.026
#> MOVER-NJ     -0.511 -0.286 -0.032
#> Wald         -0.520 -0.286 -0.052
#> Agresti-Min  -0.493 -0.286 -0.029
#> Bonett-Price -0.508 -0.286 -0.013
#> 
#> $call
#> level    cc 
#>  0.95  0.00 
#> 
# with conventional continuity adjustment
rdpairci(x = c(1, 1, 7, 12), precis = 3, cc = TRUE)
#> [[1]]
#>          Test_2
#> Test_1    Success Failure
#>   Success       1       1
#>   Failure       7      12
#> 
#> $estimates
#>                 lower    est  upper
#> SCAS_cc        -0.573 -0.286  0.039
#> SCASu_cc       -0.568 -0.286  0.031
#> Tango score_cc -0.561 -0.286  0.031
#> MOVER-NW_cc    -0.531 -0.286  0.008
#> MOVER-NJ_cc    -0.535 -0.286  0.003
#> Wald_cc        -0.529 -0.286 -0.042
#> 
#> $call
#> level    cc 
#>  0.95  0.50 
#> 
# with intermediate continuity adjustment
rdpairci(x = c(1, 1, 7, 12), precis = 3, cc = 0.25)
#> [[1]]
#>          Test_2
#> Test_1    Success Failure
#>   Success       1       1
#>   Failure       7      12
#> 
#> $estimates
#>                       lower    est  upper
#> SCAS_cc(0.25)        -0.551 -0.286  0.010
#> SCASu_cc(0.25)       -0.546 -0.286  0.003
#> Tango score_cc(0.25) -0.539 -0.286  0.003
#> MOVER-NW_cc(0.25)    -0.519 -0.286 -0.009
#> MOVER-NJ_cc(0.25)    -0.523 -0.286 -0.014
#> Wald_cc(0.25)        -0.525 -0.286 -0.047
#> 
#> $call
#> level    cc 
#>  0.95  0.25 
#> 
```
