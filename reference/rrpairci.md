# Confidence intervals for rate ratio (RR) with paired binomial rates.

Confidence intervals for comparisons of two binomial rates from paired
data. This convenience wrapper function produces a selection of the
methods below for the rate ratio (RR) contrast, with or without optional
continuity adjustment (where available).

- SCAS (skewness-corrected asymptotic score)

- SCASu (omitting the 'N-1' adjustment)

- Tang Asymptotic Score method

- MOVER-W (based on Wilson method without Newcombe correlation
  adjustment)

- MOVER-NW (based on Wilson method with Newcombe correlation adjustment)

- MOVER-NJ (based on Jeffreys method with correlation adjustment)

- Bonett-Price hybrid method

- Bonett-Price-J variant using Jeffreys intervals

- Approximate log-normal (Wald) method (strongly advise this is not used
  for any purpose but included for reference)

## Usage

``` r
rrpairci(x, level = 0.95, std_est = TRUE, cc = FALSE, precis = 8)
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

  an array containing the confidence interval for paired RR using
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
rrpairci(x = c(1, 1, 7, 12), precis = 3)
#> [[1]]
#>          Test_2
#> Test_1    Success Failure
#>   Success       1       1
#>   Failure       7      12
#> 
#> $estimates
#>                lower  est upper
#> SCAS           0.043 0.25 0.928
#> SCASu          0.043 0.25 0.898
#> Tang score     0.065 0.25 0.907
#> MOVER-W        0.069 0.25 0.869
#> MOVER-NW       0.066 0.25 0.905
#> MOVER-NJ       0.051 0.25 0.873
#> Wald           0.063 0.25 1.000
#> Bonett-Price   0.068 0.25 0.923
#> Bonett-Price-J 0.054 0.25 0.885
#> 
#> $call
#> level    cc 
#>  0.95  0.00 
#> 
# with conventional continuity adjustment
rrpairci(x = c(1, 1, 7, 12), precis = 3, cc = TRUE)
#> [[1]]
#>          Test_2
#> Test_1    Success Failure
#>   Success       1       1
#>   Failure       7      12
#> 
#> $estimates
#>                   lower  est upper
#> SCAS_cc           0.022 0.25 1.166
#> SCASu_cc          0.022 0.25 1.131
#> Tang score_cc     0.040 0.25 1.120
#> MOVER-W_cc        0.044 0.25 0.986
#> MOVER-NW_cc       0.042 0.25 1.032
#> MOVER-NJ_cc       0.030 0.25 1.013
#> Bonett-Price_cc   0.042 0.25 1.127
#> Bonett-Price-J_cc 0.031 0.25 1.106
#> 
#> $call
#> level    cc 
#>  0.95  0.50 
#> 
# with intermediate continuity adjustment
rrpairci(x = c(1, 1, 7, 12), precis = 3, cc = 0.25)
#> [[1]]
#>          Test_2
#> Test_1    Success Failure
#>   Success       1       1
#>   Failure       7      12
#> 
#> $estimates
#>                         lower  est upper
#> SCAS_cc(0.25)           0.031 0.25 1.042
#> SCASu_cc(0.25)          0.032 0.25 1.010
#> Tang score_cc(0.25)     0.052 0.25 1.009
#> MOVER-W_cc(0.25)        0.056 0.25 0.927
#> MOVER-NW_cc(0.25)       0.054 0.25 0.967
#> MOVER-NJ_cc(0.25)       0.040 0.25 0.942
#> Bonett-Price_cc(0.25)   0.055 0.25 1.020
#> Bonett-Price-J_cc(0.25) 0.042 0.25 0.992
#> 
#> $call
#> level    cc 
#>  0.95  0.25 
#> 
```
