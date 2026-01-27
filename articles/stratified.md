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
                   stratified = TRUE)
strat_rd$estimates
#>       lower    est   upper level  p1hat p2hat  p1mle p2mle
#> [1,] -0.163 -0.124 -0.0866  0.95 0.0858  0.21 0.0936 0.218
```

``` r
strat_rd$pval
#>      chisq pval2sided theta0 scorenull pval_left pval_right
#> [1,]  42.2   8.25e-11      0      -6.5  4.12e-11          1
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
#>     14.41529      8.00000      0.07156     44.50337      0.00176      0.00000 
#> pval_qualhet 
#>      0.99609
```

It is important to note that heterogeneity tests often have low power,
so it is recommended to inspect the stratum-specific estimates as
well([Greenland 1982](#ref-greenland1982)). The stratdata object
provides estimates and intervals for each stratum, along with details of
stratum weights and contributions to the Q-statistic (these may be
further examined using the `hetplot` argument).

``` r
strat_rd$stratdata
#>       x1j n1j x2j n2j p1hatj p2hatj wt_fixed wtpct_fixed wtpct_rand   theta_j
#>  [1,]  15  97  37 103 0.1546 0.3592    49.95       15.62      15.62 -2.04e-01
#>  [2,]   0   8   5  10 0.0000 0.5000     4.44        1.39       1.39 -5.00e-01
#>  [3,]  11  50  23  48 0.2200 0.4792    24.49        7.66       7.66 -2.59e-01
#>  [4,]   4 110  16 110 0.0364 0.1455    55.00       17.19      17.19 -1.10e-01
#>  [5,]   7  65   7  32 0.1077 0.2188    21.44        6.70       6.70 -1.13e-01
#>  [6,]   8  25   8  25 0.3200 0.3200    12.50        3.91       3.91 -2.98e-08
#>  [7,]   5 126  17 126 0.0397 0.1349    63.00       19.69      19.69 -9.57e-02
#>  [8,]   0 104   4  92 0.0000 0.0435    48.82       15.26      15.26 -4.51e-02
#>  [9,]   7  80  16  81 0.0875 0.1975    40.25       12.58      12.58 -1.10e-01
#>       lower_j  upper_j     V_j Stheta_j    Q_j
#>  [1,]  -0.321 -0.08495 0.00369  -0.0804 1.7506
#>  [2,]  -0.787 -0.10319 0.04111  -0.3758 3.4351
#>  [3,]  -0.435 -0.07103 0.00914  -0.1350 1.9917
#>  [4,]  -0.190 -0.03530 0.00151   0.0151 0.1517
#>  [5,]  -0.290  0.04065 0.00704   0.0132 0.0246
#>  [6,]  -0.260  0.25979 0.01763   0.1242 0.8753
#>  [7,]  -0.169 -0.02700 0.00132   0.0290 0.6349
#>  [8,]  -0.101 -0.00470 0.00119   0.0807 5.4843
#>  [9,]  -0.221 -0.00172 0.00300   0.0142 0.0671
```

### Random effects

In this instance, there is weak evidence ($p = 0.07$) that the treatment
effect (on the RD scale) varies across strata. (Estimation on the RR
contrast scale gives more consistency across strata.) For the sake of
illustration, a random effects analysis of RD would be obtained as
follows, giving a wider interval that incorporates the stratum
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
                   prediction = TRUE)
strat_rd_rand$estimates
#>       lower    est upper level  p1hat p2hat  p1mle p2mle
#> [1,] -0.188 -0.124 -0.06  0.95 0.0858  0.21 0.0936 0.218
```

Note that the TDAS method does not include a skewness correction, so
when `random = TRUE`, the `skew` argument only affects the `stratdata`
output element.

A prediction interval for the treatment effect in a new study ([Higgins,
Thompson, and Spiegelhalter 2008](#ref-higgins2008)) can be obtained
using `prediction = TRUE`:

``` r
strat_rd_rand$prediction
#>       lower  upper
#> [1,] -0.244 0.0158
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

Roderick, P, G Ferris, K Wilson, H Halls, D Jackson, R Collins, and C
Baigent. 2005. “Towards Evidence-Based Guidelines for the Prevention of
Venous Thromboembolism: Systematic Reviews of Mechanical Methods, Oral
Anticoagulation, Dextran and Regional Anaesthesia as
Thromboprophylaxis.” *Health Technology Assessment* 9 (49).
<https://doi.org/10.3310/hta9490>.
