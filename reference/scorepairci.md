# Score confidence intervals and tests for comparisons of paired binomial rates.

Confidence intervals for the rate (or risk) difference ("RD"), rate
ratio ("RR") or conditional odds ratio ("OR"), for paired binomial data.
(For paired Poisson rates, suggest use the tdasci function with
`distrib = "poi"`, and `weighting = "MH"`, with pairs as strata.) This
function applies the score-based Tango and Tang methods for RD and RR
respectively, with iterative and closed-form versions, and an added
skewness correction for improved one-sided coverage. For OR, intervals
are produced based on transforming the SCASp interval for the single
proportion. Includes options for continuity adjustment, where the
magnitude of adjustment can be customised.

## Usage

``` r
scorepairci(
  x,
  level = 0.95,
  contrast = "RD",
  bcf = TRUE,
  skew = TRUE,
  closedform = FALSE,
  cc = FALSE,
  theta0 = NULL,
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
  "RR" = rate ratio;  
  "OR" = conditional odds ratio.

- bcf:

  Logical (default TRUE) indicating whether to apply 'N-1' variance bias
  correction in the score denominator. (Under evaluation, manuscript
  under review.)

- skew:

  Logical (default TRUE) indicating whether to apply skewness correction
  or not. (Under evaluation, manuscript under review.)

- closedform:

  Logical (default FALSE) indicating whether to use closed form
  calculation (only available if `skew = FALSE`)

- cc:

  Number or logical (default FALSE) specifying (amount of) continuity
  adjustment. cc = 0.5 corresponds to the continuity-corrected McNemar
  test.

- theta0:

  Number to be used in a one-sided significance test (e.g.
  non-inferiority margin). 1-sided p-value will be \< 0.025 iff 2-sided
  95\\ excludes theta0. NB: can also be used for a superiority test by
  setting theta0 = 0.

- precis:

  Number (default 6) specifying precision (i.e. number of decimal
  places) to be used in optimisation subroutine for the confidence
  interval.

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

- pval:

  the corresponding 2-sided significance test against the null
  hypothesis that p_1 = p_2, and one-sided significance tests against
  the null hypothesis that theta \>= or \<= theta0 as specified.

- call:

  details of the function call.

## References

Tango T. Equivalence test and confidence interval for the difference in
proportions for the paired-sample design. Statistics in Medicine 1998;
17:891-908

Newcombe RG. Improved confidence intervals for the difference between
binomial proportions based on paired data. Statistics in Medicine 1998;
17:2635-2650

Tango T. Improved confidence intervals for the difference between
binomial proportions based on paired data by Robert G. Newcombe,
Statistics in Medicine, 17, 2635-2650 (1998). Statistics in Medicine
1999; 18(24):3511-3513

Nam J-M, Blackwelder WC. Analysis of the ratio of marginal probabilities
in a matched-pair setting. Stat Med 2002; 21(5):689–699

Tang N-S, Tang M-L, Chan ISF. On tests of equivalence via non-unity
relative risk for matched-pair design. Statistics in Medicine 2003;
22:1217-1233

Yang Z, Sun X and Hardin JW. A non-iterative implementation of Tango's
score confidence interval for a paired difference of proportions.
Statistics in Medicine 2013; 32:1336-1342

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

Pete Laud, <p.j.laud@sheffield.ac.uk>

## Examples

``` r
# Example data from Fagerland et al 2014
# SCAS method for RD
scorepairci(x = c(1, 1, 7, 12), contrast = "RD")
#> $data
#>          Test_2
#> Test_1    Success Failure
#>   Success       1       1
#>   Failure       7      12
#> 
#> $estimates
#>          lower       est     upper level    p1hat    p2hat    p1mle    p2mle
#> [1,] -0.528113 -0.285861 -0.018422  0.95 0.095238 0.380952 0.095201 0.381062
#>       phi_hat phi_c  psi_hat
#> [1,] 0.079536     0 1.714286
#> 
#> $pval
#>         chisq pval2sided theta0 scorenull  pval_left pval_right
#> [1,] 4.285714 0.03843393      0 -2.070197 0.01921697   0.980783
#> 
#> $call
#>   contrast      level        bcf       skew         cc closedform 
#>       "RD"     "0.95"     "TRUE"     "TRUE"    "FALSE"    "FALSE" 
#> 
# Tango method
scorepairci(x = c(1, 1, 7, 12), contrast = "RD", skew = FALSE, bcf = FALSE)
#> $data
#>          Test_2
#> Test_1    Success Failure
#>   Success       1       1
#>   Failure       7      12
#> 
#> $estimates
#>          lower       est     upper level    p1hat    p2hat    p1mle    p2mle
#> [1,] -0.517232 -0.285714 -0.026003  0.95 0.095238 0.380952 0.095238 0.380952
#>       phi_hat phi_c  psi_hat
#> [1,] 0.079536     0 1.714286
#> 
#> $pval
#>      chisq pval2sided theta0 scorenull  pval_left pval_right
#> [1,]   4.5 0.03389485      0  -2.12132 0.01694743  0.9830526
#> 
#> $call
#>   contrast      level        bcf       skew         cc closedform 
#>       "RD"     "0.95"    "FALSE"    "FALSE"    "FALSE"    "FALSE" 
#> 
# SCAS for RR
scorepairci(x = c(1, 1, 7, 12), contrast = "RR")
#> $data
#>          Test_2
#> Test_1    Success Failure
#>   Success       1       1
#>   Failure       7      12
#> 
#> $estimates
#>         lower      est    upper level    p1hat    p2hat    p1mle    p2mle
#> [1,] 0.042901 0.262723 0.928199  0.95 0.095238 0.380952 0.099445 0.378517
#>       phi_hat phi_c  psi_hat
#> [1,] 0.079536     0 1.714286
#> 
#> $pval
#>         chisq pval2sided theta0 scorenull  pval_left pval_right
#> [1,] 4.285714 0.03843393      1 -2.070197 0.01921697   0.980783
#> 
#> $call
#>   contrast      level        bcf       skew         cc closedform 
#>       "RR"     "0.95"     "TRUE"     "TRUE"    "FALSE"    "FALSE" 
#> 
# Tang method
scorepairci(x = c(1, 1, 7, 12), contrast = "RR", skew = FALSE, bcf = FALSE)
#> $data
#>          Test_2
#> Test_1    Success Failure
#>   Success       1       1
#>   Failure       7      12
#> 
#> $estimates
#>         lower  est    upper level    p1hat    p2hat    p1mle    p2mle  phi_hat
#> [1,] 0.065279 0.25 0.906881  0.95 0.095238 0.380952 0.095238 0.380952 0.079536
#>      phi_c  psi_hat
#> [1,]     0 1.714286
#> 
#> $pval
#>      chisq pval2sided theta0 scorenull  pval_left pval_right
#> [1,]   4.5 0.03389485      1  -2.12132 0.01694743  0.9830526
#> 
#> $call
#>   contrast      level        bcf       skew         cc closedform 
#>       "RR"     "0.95"    "FALSE"    "FALSE"    "FALSE"    "FALSE" 
#> 
# Transformed SCASp method for OR
scorepairci(x = c(1, 1, 7, 12), contrast = "OR")
#> $data
#>          Test_2
#> Test_1    Success Failure
#>   Success       1       1
#>   Failure       7      12
#> 
#> $estimates
#>         lower      est    upper
#> [1,] 0.007702 0.161863 0.912316
#> 
#> $pval
#>         chisq pval2sided theta0 scorenull  pval_left pval_right
#> [1,] 4.285714 0.03843393      1 -2.070197 0.01921697   0.980783
#> 
#> $call
#>   contrast      level        bcf       skew         cc closedform 
#>       "OR"     "0.95"     "TRUE"     "TRUE"    "FALSE"    "FALSE" 
#> 
# Transformed Wilson score method
scorepairci(x = c(1, 1, 7, 12), contrast = "OR", skew = FALSE, bcf = FALSE)
#> $data
#>          Test_2
#> Test_1    Success Failure
#>   Success       1       1
#>   Failure       7      12
#> 
#> $estimates
#>         lower      est    upper
#> [1,] 0.022931 0.142857 0.889959
#> 
#> $pval
#>      chisq pval2sided theta0 scorenull  pval_left pval_right
#> [1,]   4.5 0.03389485      1  -2.12132 0.01694743  0.9830526
#> 
#> $call
#>   contrast      level        bcf       skew         cc closedform 
#>       "OR"     "0.95"    "FALSE"    "FALSE"    "FALSE"    "FALSE" 
#> 
```
