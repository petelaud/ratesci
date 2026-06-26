# MOVER confidence intervals for comparisons of paired binomial rates.

Confidence intervals for the rate (or risk) difference ("RD"), or rate
ratio ("RR"), for paired binomial data. This function applies the Method
of Variance Estimates Recovery (MOVER) for RD and RR, incorporating
Newcombe's correlation correction. All methods have options for
continuity adjustment, where the magnitude of adjustment can be
customised.

## Usage

``` r
moverpairci(
  x,
  level = 0.95,
  contrast = "RD",
  type = "jeff",
  corc = TRUE,
  cc = FALSE,
  precis = 6,
  warn = TRUE,
  ...
)
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

- contrast:

  Character string indicating the contrast of interest:  
  "RD" = rate difference (default);  
  "RR" = rate ratio.

- type:

  Character string indicating the method used for the intervals for the
  individual group rates.  
  "jeff" = Jeffreys equal-tailed intervals (default);  
  "SCASp" = skewness-corrected score,  
  "midp" = mid-p,  
  "wilson" = Wilson score (not recommended, known to be skewed).

- corc:

  Logical (default TRUE) indicating whether to apply adjustment to the
  correlation estimate from Newcombe.

- cc:

  Number or logical (default FALSE) specifying (amount of) continuity
  adjustment. When a score-based method is used, cc = 0.5 corresponds to
  the continuity-corrected McNemar test.

- precis:

  Number (default 6) specifying precision (i.e. number of decimal
  places) to be used in output.

- warn:

  Logical (default TRUE) giving the option to suppress warnings.

- ...:

  Other arguments.

## Value

A list containing the following components:

- data:

  the input data in 2x2 matrix form.

- estimates:

  the requested contrast, with its confidence interval and the specified
  confidence level, along with estimates of the marginal probabilities
  and the correlation coefficient (uncorrected and corrected).

- call:

  details of the function call.

## References

Newcombe RG. Improved confidence intervals for the difference between
binomial proportions based on paired data. Statistics in Medicine 1998;
17:2635-2650

Tang M-L, Li H-Q, Tang N-S. Confidence interval construction for
proportion ratio in paired studies based on hybrid method. Statistical
Methods in Medical Research 2010; 21(4):361-378

Tang N-S et al. Asymptotic confidence interval construction for
proportion difference in medical studies with bilateral data.
Statistical Methods in Medical Research. 2011; 20(3):233-259

Fagerland MW, Lydersen S, Laake P. Recommended tests and confidence
intervals for paired binomial proportions. Statistics in Medicine 2014;
33(16):2850-2875

Laud PJ. Equal-tailed confidence intervals for comparison of rates.
Pharmaceutical Statistics 2017; 16:334-348.

DelRocco N et al. New Confidence Intervals for Relative Risk of Two
Correlated Proportions. Statistics in Biosciences 2023; 15:1–30

Chang P et al. Continuity corrected score confidence interval for the
difference in proportions in paired data. Journal of Applied Statistics
2024; 51-1:139-152

Laud PJ. Comments on "New Confidence Intervals for Relative Risk of Two
Correlated Proportions" (2023). Statistics in Biosciences 2025;
https://doi.org/10.1007/s12561-025-09479-4

Laud PJ. Improved confidence intervals and tests for paired binomial
proportions. (2026, Under review)

## Author

Pete Laud, <pete@sheffstat.co.uk>

## Examples

``` r
# Example data from Fagerland et al 2014
# MOVER-NJ method
moverpairci(x = c(1, 1, 7, 12), contrast = "RD", corc = TRUE, type = "jeff")
#> $data
#>          Test_2
#> Test_1    Success Failure
#>   Success       1       1
#>   Failure       7      12
#> 
#> $estimates
#>          lower       est     upper level    p1hat    p2hat phi_hat
#> [1,] -0.510506 -0.285714 -0.032389  0.95 0.095238 0.380952       0
#> 
#> $call
#> contrast     type    level       cc 
#>     "RD"   "jeff"   "0.95"  "FALSE" 
#> 
# MOVER-NJ
moverpairci(x = c(1, 1, 7, 12), contrast = "RR", corc = TRUE, type = "jeff")
#> $data
#>          Test_2
#> Test_1    Success Failure
#>   Success       1       1
#>   Failure       7      12
#> 
#> $estimates
#>         lower  est    upper level    p1hat    p2hat phi_hat
#> [1,] 0.051297 0.25 0.873051  0.95 0.095238 0.380952       0
#> 
#> $call
#> contrast     type    level       cc 
#>     "RR"   "jeff"   "0.95"  "FALSE" 
#> 
```
