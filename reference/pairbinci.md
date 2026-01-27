# Confidence intervals for comparisons of paired binomial rates.

Confidence intervals for the rate (or risk) difference ("RD"), rate
ratio ("RR") or conditional odds ratio ("OR"), for paired binomial data.
(For paired Poisson rates, suggest use the tdasci function with
`distrib = "poi"`, and `weighting = "MH"`, with pairs as strata.) This
function applies the score-based Tango and Tang methods for RD and RR
respectively, with iterative and closed-form versions, and an added
skewness correction for improved one-sided coverage. Also includes MOVER
options using the Method of Variance Estimates Recovery for paired RD
and RR, incorporating Newcombe's correlation correction, and some
simpler methods by Bonett & Price for RD and RR. For OR, intervals are
produced based on transforming various intervals for the single
proportion, including SCASp, mid-p and Jeffreys. All methods have
options for continuity adjustment, and the magnitude of adjustment can
be customised.

## Usage

``` r
pairbinci(
  x,
  level = 0.95,
  contrast = "RD",
  method = ifelse(contrast == "OR", "SCASp", "Score"),
  moverbase = ifelse(method %in% c("MOVER", "MOVER_newc", "BP"), "jeff", NULL),
  bcf = TRUE,
  skew = TRUE,
  cc = FALSE,
  theta0 = NULL,
  precis = 6,
  warn = TRUE,
  method_RD = NULL,
  method_RR = NULL,
  method_OR = NULL,
  cctype = NULL,
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

- method:

  Character string indicating the confidence interval method to be used.
  The following are available for `contrast = "RD"` or `"RR"`:  
  "Score" = (default) asymptotic score class of methods including Tango
  (for RD) / Tang (for RR), by iterative calculations, with optional
  skewness correction;  
  "Score_closed" = closed form solution for Tango/Tang intervals
  (without skewness correction);  
  "MOVER" = hybrid MOVER method (as per "method 8" in Newcombe, but with
  a choice of input methods - see moverbase);  
  "MOVER_newc" = hybrid MOVER methods with correction to correlation
  estimate (Newcombe's "method 10");  
  "TDAS" = t-distribution asymptotic score (experimental method, now
  deprecated);  
  "BP" = Wald with Bonett-Price adjustment for RD, or Hybrid
  Bonett-Price method for RR.  
  For `contrast = "OR"`, one of the following methods may be selected,
  all of which are based on transformation of an interval for a single
  proportion `b/(b+c)`:  
  "SCASp" = transformed skewness-corrected score (default);  
  "jeff" = transformed Jeffreys;  
  "midp" = transformed mid-p;  
  "wilson" = transformed Wilson score - included for reference only, not
  recommended.

- moverbase:

  Character string indicating the base method used as input for the
  MOVER methods for RD or RR (when method = "MOVER" or "MOVER_newc"),
  and for the Hybrid BP method for RR: "jeff" = Jeffreys equal-tailed
  interval (default), "SCASp" = skewness-corrected score, "midp" =
  mid-p, "wilson" = Wilson score (not recommended, known to be skewed).

- bcf:

  Logical (default FALSE) indicating whether to apply variance bias
  correction in the score denominator. (Under evaluation, manuscript
  under review.)

- skew:

  Logical (default TRUE) indicating whether to apply skewness correction
  or not. (Under evaluation, manuscript under review.)

  - Only applies for the iterative `method = "Score"`.

- cc:

  Number or logical (default FALSE) specifying (amount of) continuity
  adjustment. When a score-based method is used, cc = 0.5 corresponds to
  the continuity-corrected McNemar test.

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

- method_RD:

  (deprecated: parameter renamed to method)

- method_RR:

  (deprecated: parameter renamed to method)

- method_OR:

  (deprecated: parameter renamed to method)

- cctype:

  (deprecated: new equivariant cc method implemented instead.)

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

Agresti A, Min Y. Simple improved confidence intervals for comparing
matched proportions. Statistics in Medicine 2005; 24:729-740

Bonett DG, Price RM. Confidence intervals for a ratio of binomial
proportions based on paired data. Statistics in Medicine 2006;
25:3039-3047

Tang M-L, Li H-Q, Tang N-S. Confidence interval construction for
proportion ratio in paired studies based on hybrid method. Statistical
Methods in Medical Research 2010; 21(4):361-378

Tang N-S et al. Asymptotic confidence interval construction for
proportion difference in medical studies with bilateral data.
Statistical Methods in Medical Research. 2011; 20(3):233-259

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
proportions. (2025, Under review)

## Author

Pete Laud, <p.j.laud@sheffield.ac.uk>

## Examples

``` r
# Example from Fagerland et al 2014
# SCAS method for RD
pairbinci(x = c(1, 1, 7, 12), contrast = "RD", method = "Score")
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
#> contrast   method    level      bcf     skew       cc 
#>     "RD"  "Score"   "0.95"   "TRUE"   "TRUE"  "FALSE" 
#> 
# Tango method
pairbinci(x = c(1, 1, 7, 12), contrast = "RD", method = "Score", skew = FALSE, bcf = FALSE)
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
#> contrast   method    level      bcf     skew       cc 
#>     "RD"  "Score"   "0.95"  "FALSE"  "FALSE"  "FALSE" 
#> 
# MOVER-NJ method
pairbinci(x = c(1, 1, 7, 12), contrast = "RD", method = "MOVER_newc", moverbase = "jeff")
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
#>     contrast       method    moverbase        level          bcf         skew 
#>         "RD" "MOVER_newc"       "jeff"       "0.95"       "TRUE"       "TRUE" 
#>           cc 
#>      "FALSE" 
#> 
# SCAS for RR
pairbinci(x = c(1, 1, 7, 12), contrast = "RR", method = "Score")
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
#> contrast   method    level      bcf     skew       cc 
#>     "RR"  "Score"   "0.95"   "TRUE"   "TRUE"  "FALSE" 
#> 
# Tang method
pairbinci(x = c(1, 1, 7, 12), contrast = "RR", method = "Score", skew = FALSE, bcf = FALSE)
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
#> contrast   method    level      bcf     skew       cc 
#>     "RR"  "Score"   "0.95"  "FALSE"  "FALSE"  "FALSE" 
#> 
# MOVER-NJ
pairbinci(x = c(1, 1, 7, 12), contrast = "RR", method = "MOVER_newc", moverbase = "jeff")
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
#>     contrast       method    moverbase        level          bcf         skew 
#>         "RR" "MOVER_newc"       "jeff"       "0.95"       "TRUE"       "TRUE" 
#>           cc 
#>      "FALSE" 
#> 
# Transformed SCASp method for OR
pairbinci(x = c(1, 1, 7, 12), contrast = "OR", method = "SCASp")
#> $data
#>          Test_2
#> Test_1    Success Failure
#>   Success       1       1
#>   Failure       7      12
#> 
#> $estimates
#>         lower      est    upper
#> [1,] 0.007702 0.161863 0.912315
#> 
#> $pval
#>         chisq pval2sided theta0 scorenull  pval_left pval_right
#> [1,] 4.285714 0.03843393      1 -2.070197 0.01921697   0.980783
#> 
#> $call
#> contrast   method    level      bcf     skew       cc 
#>     "OR"  "SCASp"   "0.95"   "TRUE"   "TRUE"  "FALSE" 
#> 
# Transformed Wilson method
pairbinci(x = c(1, 1, 7, 12), contrast = "OR", method = "wilson")
#> $data
#>          Test_2
#> Test_1    Success Failure
#>   Success       1       1
#>   Failure       7      12
#> 
#> $estimates
#>         lower      est   upper
#> [1,] 0.022932 0.142857 0.88996
#> 
#> $call
#> contrast   method    level      bcf     skew       cc 
#>     "OR" "wilson"   "0.95"   "TRUE"   "TRUE"  "FALSE" 
#> 
```
