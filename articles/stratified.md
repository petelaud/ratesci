# 3: Stratified contrasts

``` r

library(ratesci)
```

## Confidence intervals for stratified comparisons of independent binomial or Poisson rates

Stratified analysis might be required for combining results from
multiple studies in a meta-analysis, or adjusting clinical trial results
for factors used in a stratified randomisation. Associated tests against
a null hypothesis of common effects over strata (homogeneity) can also
be useful for examining subgroup effects in a clinical trial. The
illustration below uses data from a meta-analysis of 9 trials studying
the effectiveness of graduated compression stockings for prevention of
postoperative deep vein thrombosis (DVT) (from [Roderick et al.
2005](#ref-roderick2005)):

``` r

data(compress, package = "ratesci")
strat_rd <- scoreci(x1 = compress$event.gcs, 
                   n1 = compress$n.gcs, 
                   x2 = compress$event.control, 
                   n2 = compress$n.control, 
                   contrast = "RD", 
                   stratified = TRUE,
                   precis = 4)
strat_rd$estimates
#>        lower     est   upper level  p1hat p2hat  p1mle  p2mle
#> [1,] -0.1632 -0.1242 -0.0866  0.95 0.0858  0.21 0.0936 0.2179
```

``` r

strat_rd$pval
#>         chisq  pval2sided theta0 scorenull   pval_left pval_right
#> [1,] 42.19806 8.24816e-11      0 -6.496004 4.12408e-11          1
```

The Qtest object provides a heterogeneity test and related quantities,
including an assessment of qualitative heterogeneity (indicating whether
treatment effects are in opposite directions in different strata). All
of which are derived using the score statistic and its estimated
variance, for a consistent approach - see ([Laud 2017](#ref-laud2017)
Supplementary Appendix S4).

``` r

strat_rd$Qtest
#>            Q         Q_df     pval_het           I2         tau2           Qc 
#> 14.415376544  8.000000000  0.071560846 44.503704254  0.001756709  0.000000000 
#> pval_qualhet 
#>  0.996093750
```

It is important to note that heterogeneity tests often have low power,
so it is recommended to inspect the stratum-specific estimates as
well([Greenland 1982](#ref-greenland1982)). The stratdata object
provides estimates and intervals for each stratum, along with details of
stratum weights and contributions to the Q-statistic (these may be
further examined using the `hetplot` argument).

``` r

strat_rd$stratdata
#>       x1j n1j x2j n2j p1hatj p2hatj wt_fixed wtpct_fixed wtpct_rand theta_j
#>  [1,]  15  97  37 103 0.1546 0.3592  49.9550     15.6159    15.6159 -0.2044
#>  [2,]   0   8   5  10 0.0000 0.5000   4.4444      1.3893     1.3893 -0.5000
#>  [3,]  11  50  23  48 0.2200 0.4792  24.4898      7.6555     7.6555 -0.2585
#>  [4,]   4 110  16 110 0.0364 0.1455  55.0000     17.1930    17.1930 -0.1096
#>  [5,]   7  65   7  32 0.1077 0.2188  21.4433      6.7032     6.7032 -0.1129
#>  [6,]   8  25   8  25 0.3200 0.3200  12.5000      3.9075     3.9075  0.0000
#>  [7,]   5 126  17 126 0.0397 0.1349  63.0000     19.6938    19.6938 -0.0957
#>  [8,]   0 104   4  92 0.0000 0.0435  48.8163     15.2600    15.2600 -0.0451
#>  [9,]   7  80  16  81 0.0875 0.1975  40.2484     12.5817    12.5817 -0.1103
#>       lower_j upper_j    V_j Stheta_j    Q_j
#>  [1,] -0.3211 -0.0849 0.0037  -0.0804 1.7506
#>  [2,] -0.7868 -0.1032 0.0411  -0.3758 3.4351
#>  [3,] -0.4354 -0.0710 0.0091  -0.1350 1.9917
#>  [4,] -0.1901 -0.0353 0.0015   0.0151 0.1517
#>  [5,] -0.2901  0.0407 0.0070   0.0132 0.0246
#>  [6,] -0.2598  0.2598 0.0176   0.1242 0.8753
#>  [7,] -0.1694 -0.0270 0.0013   0.0290 0.6349
#>  [8,] -0.1012 -0.0047 0.0012   0.0807 5.4844
#>  [9,] -0.2211 -0.0017 0.0030   0.0142 0.0671
```

### Random effects

In this instance, there is weak evidence ($`p = 0.07`$) that the
treatment effect (on the RD scale) varies across strata. (Estimation on
the RR contrast scale gives more consistency across strata.) For the
sake of illustration, a random effects analysis of RD would be obtained
as follows, giving a wider interval that incorporates the stratum
variability. The `random = TRUE` option invokes the t-distribution
asymptotic score (TDAS) method ([Laud 2017](#ref-laud2017)), which
modifies the asymptotic score methodology using formulae from the
Hartung-Knapp-Sidik-Jonkman (HKSJ) random effects meta-analysis method
(with superior performance compared to the DerSimonian-Laird method).

``` r

strat_rd_rand <- scoreci(x1 = compress$event.gcs, 
                   n1 = compress$n.gcs, 
                   x2 = compress$event.control, 
                   n2 = compress$n.control, 
                   contrast = "RD", 
                   stratified = TRUE,
                   random = TRUE,
                   prediction = TRUE,
                   precis = 4)
strat_rd_rand$estimates
#>        lower     est upper level  p1hat p2hat  p1mle  p2mle
#> [1,] -0.1884 -0.1242 -0.06  0.95 0.0858  0.21 0.0936 0.2178
```

Note that the TDAS method does not include a skewness correction, so
when `random = TRUE`, the `skew` argument only affects the `stratdata`
output element.

A prediction interval for the treatment effect in a new study ([Higgins
et al. 2008](#ref-higgins2008)) can be obtained using
`prediction = TRUE`:

``` r

strat_rd_rand$prediction
#>        lower  upper
#> [1,] -0.2436 0.0158
```

## References

Greenland, Sander. 1982. “Interpretation and Estimation of Summary
Ratios Under Heterogeneity.” *Statistics in Medicine* 1 (3): 217–27.
<https://doi.org/10.1002/sim.4780010304>.

Higgins, Julian P. T., Simon G. Thompson, and David J. Spiegelhalter.
2008. “A Re-Evaluation of Random-Effects Meta-Analysis.” *Journal of the
Royal Statistical Society Series A: Statistics in Society* 172 (1):
137–59. <https://doi.org/10.1111/j.1467-985x.2008.00552.x>.

Laud, Peter J. 2017. “Equal-Tailed Confidence Intervals for Comparison
of Rates.” *Pharmaceutical Statistics* 16 (5): 334–48.
<https://doi.org/10.1002/pst.1813>.

Roderick, P, G Ferris, K Wilson, et al. 2005. “Towards Evidence-Based
Guidelines for the Prevention of Venous Thromboembolism: Systematic
Reviews of Mechanical Methods, Oral Anticoagulation, Dextran and
Regional Anaesthesia as Thromboprophylaxis.” *Health Technology
Assessment* 9 (49). <https://doi.org/10.3310/hta9490>.
