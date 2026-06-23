# 2: Basic contrasts (single stratum)

``` r

library(ratesci)
```

## Confidence intervals for comparisons of independent binomial or Poisson rates

### SCAS and other asymptotic score methods

For comparison of two groups when the outcome is a binomial proportion,
you may wish to quantify the effect size using the difference
($`\hat \theta_{RD} = \hat p_1 - \hat p_2`$, where
$`\hat p_i = x_i/n_i`$) or the ratio of the two proportions
($`\hat \theta_{RR} = \hat p_1 / \hat p_2`$). Another option is the odds
ratio
($`\hat \theta_{OR} = \hat p_1 (1 - \hat p_2) / \hat p_2(1 - \hat p_1)`$),
which is less easy to interpret, but is widely used due to its useful
properties for logistic regression modelling or case-control studies.

To calculate a confidence interval (CI) for any of these contrasts, the
skewness-corrected asymptotic score (SCAS) method is recommended, as one
that succeeds, on average, at containing the true parameter $`\theta`$
with the appropriate nominal probability (e.g. 95%), and has evenly
distributed tail probabilities ([Laud 2017](#ref-laud2017)). It is a
modified version of the Miettinen-Nurminen (MN) asymptotic score method
([Miettinen and Nurminen 1985](#ref-miettinen1985)). MN and SCAS each
represent a family of methods encompassing analysis of RD and RR for
binomial proportions or Poisson rates (e.g. exposure-adjusted incidence
rate), and binomial OR.

The plots below illustrate the one-sided and two-sided interval coverage
probabilities[^1] achieved by SCAS compared to some other popular
methods[^2] with $`(n_1, n_2) = (45, 15)`$, with the second row of each
plot using moving average smoothing. A selection of coverage probability
plots for other sample size combinations and other contrasts can be
found in the “plots” folder of the [ratesci GitHub
repository](https://github.com/petelaud/ratesci/tree/master/plots).

![](images/summaryRD95_45%2C15os_small.jpg)

Differences in two-sided coverage can be more subtle:

![](images/summaryRD95_45%2C15_small.jpg)

The skewness correction, introduced by ([Gart and Nam
1988](#ref-gart1988)) should be considered essential for analysis of
relative risk, but can also have a substantial effect for RD as seen
above. The one-sided coverage of MN for RR can be poor even in large
samples with equal-sized groups, which means that the confidence
interval tends to be located too close to 1:

![](images/binRR_95_100%2C100os_small.jpg)

\[Note: ([Fagerland et al. 2011](#ref-fagerland2011)) dismissed the MN
method, saying that it can be somewhat liberal. In fact, based on
reconstruction of their results, it appears that their evaluation was of
the closely-related but inferior Mee method, which omits the N/(N-1)
variance bias correction. This error was corrected in their 2017 book.
Essentially, it seems they reject the MN method on the basis of some
very small regions of the parameter space where coverage is slightly
below nominal. These dips in coverage do occur, but fluctuate above and
below the nominal level in a pattern of ridges and furrows aligned along
the diagonal of a fixed value of $`\theta_{RD}`$. The most severe dips
tend to occur when sample sizes are multiples of 10
(e.g. $`(n_1, n_2) = (40, 10)`$), which might be relatively rare in
practice. The moving average smoothing in the above surface plots is
designed to even out these fluctuations in coverage. The ‘N-1’
adjustment in the MN method elevates coverage probabilities slightly,
and the skewness correction gives a further improvement overall, as well
as substantially improving one-sided coverage when group sizes are
unequal.\]

[`scoreci()`](https://petelaud.github.io/ratesci/reference/scoreci.md)
is used as follows, for example to obtain a SCAS 95% CI for the
difference `5/56 - 0/29`:

``` r

out <- scoreci(x1 = 5, n1 = 56, x2 = 0, n2 = 29, precis = 4)
out$estimates
#>        lower    est  upper level x1 n1 x2 n2  p1hat p2hat  p1mle p2mle
#> [1,] -0.0186 0.0917 0.1867  0.95  5 56  0 29 0.0893     0 0.0917     0
```

The underlying z-statistic is used to obtain a two-sided hypothesis test
for association (`pval2sided`). Note that under certain conditions
(i.e. $`n_1 = n_2`$) this is equivalent to the Egon Pearson ‘N-1’
chi-squared test ([Campbell 2007](#ref-campbell2007)). The facility is
also provided for a custom one-sided test against any specified null
hypothesis value $`\theta_0`$, e.g. for non-inferiority testing
(`pval_left` and `pval_right`). See the [tests
vignette](https://petelaud.github.io/ratesci/articles/tests.html) for
more details.

``` r

out$pval
#>         chisq pval2sided theta0 scorenull pval_left pval_right
#> [1,] 3.024862 0.08199729      0  1.739213 0.9590014 0.04099865
```

The
[`scoreci()`](https://petelaud.github.io/ratesci/reference/scoreci.md)
function provides options to omit the skewness correction
(`skew = FALSE` for the MN method) and the variance bias correction
(`bcf = FALSE` for the Mee method, which would be consistent with the
Karl Pearson unadjusted chi-squared test). To keep things simple when
choosing the SCAS method, you may instead use the
[`scasci()`](https://petelaud.github.io/ratesci/reference/scasci.md)
wrapper function instead of
[`scoreci()`](https://petelaud.github.io/ratesci/reference/scoreci.md):

``` r

scasci(x1 = 5, n1 = 56, x2 = 0, n2 = 29, precis = 4)$estimates
#>        lower    est  upper level x1 n1 x2 n2  p1hat p2hat  p1mle p2mle
#> [1,] -0.0186 0.0917 0.1867  0.95  5 56  0 29 0.0893     0 0.0917     0
```

For analysis of relative risk, use:

``` r

scasci(x1 = 5, n1 = 56, x2 = 0, n2 = 29, contrast = "RR", precis = 4)$estimates
#>       lower     est upper level x1 n1 x2 n2  p1hat p2hat  p1mle  p2mle
#> [1,] 0.7696 16.3634   Inf  0.95  5 56  0 29 0.0893     0 0.0868 0.0053
```

And for odds ratio:

``` r

scasci(x1 = 5, n1 = 56, x2 = 0, n2 = 29, contrast = "OR", precis = 4)$estimates
#>       lower    est upper level x1 n1 x2 n2  p1hat p2hat  p1mle  p2mle
#> [1,] 0.7587 17.665   Inf  0.95  5 56  0 29 0.0893     0 0.0865 0.0053
```

For analysis of Poisson incidence rates, for example exposure-adjusted
adverse event rates, the same methodology is adapted for the Poisson
distribution, with the `distrib` argument. So, for example if the
denominators 56 and 29 represent the number of patient-years at risk in
each group, instead of the number of patients:

``` r

scasci(x1 = 5, n1 = 56, x2 = 0, n2 = 29, contrast = "RR", distrib = "poi", precis = 4)$estimates
#>       lower     est upper level x1 n1 x2 n2  p1hat p2hat  p1mle  p2mle
#> [1,] 0.7264 16.0536   Inf  0.95  5 56  0 29 0.0893     0 0.0865 0.0054
```

### MOVER methods

An alternative family of methods to produce confidence intervals for the
same set of binomial and Poisson contrasts is the Method of Variance
Estimates Recovery (MOVER), also known as Square-and-Add (For technical
details see chapter 7 of [Newcombe 2012](#ref-newcombe2012)). Originally
labelled as a “Score” method ([Newcombe 1998](#ref-newcombe1998)) due to
the involvement of the Wilson Score interval, this involves a
combination of intervals calculated separately for $`p_1`$ and $`p_2`$.
Newcombe based his method on Wilson intervals (this version is labelled
as “MOVER-W” in the ratesci package), but coverage can be slightly
improved by using equal-tailed Jeffreys intervals instead (“MOVER-J”)
([Laud 2017](#ref-laud2017))[^3]. Coverage properties remain generally
inferior to SCAS - with larger sample sizes, MOVER-J has two-sided
coverage close to (but slightly below) the nominal level, but central
location is not achieved (except for the RR contrast).

There is no corresponding hypothesis test for the MOVER methods.

A MOVER-J interval for binomial RD from the same observed data as above
would be obtained using:

``` r

moverci(x1 = 5, n1 = 56, x2 = 0, n2 = 29)$estimates
#>             lower        est     upper level x1 n1 x2 n2      p1hat       p2hat
#> [1,] -0.009734968 0.08403347 0.1772258  0.95  5 56  0 29 0.09177968 0.007746206
```

For comparison, Newcombe’s version using Wilson intervals is:

``` r

moverci(x1 = 5, n1 = 56, x2 = 0, n2 = 29, type = "wilson")$estimates
#>            lower        est   upper level x1 n1 x2 n2      p1hat p2hat
#> [1,] -0.03813715 0.08928571 0.19256  0.95  5 56  0 29 0.08928571     0
```

The MOVER approach may be further modified by adapting the Jeffreys
interval (which uses a non-informative $`Beta(0.5, 0.5)`$ prior for the
group rates) to incorporate prior information about $`p_1`$ and $`p_2`$
for an approximate Bayesian interval. For example, given data from a
pilot study with 1/10 events in the first group and 0/10 in the second
group, the non-informative $`Beta(0.5, 0.5)`$ prior distributions for
$`p_1`$ and $`p_2`$ might be updated as follows:

``` r

moverbci(x1 = 5, n1 = 56, x2 = 0, n2 = 29, a1 = 1.5, b1 = 9.5, a2 = 0.5, b2 = 10.5)
#> $estimates
#>            lower        est     upper level x1 n1 x2 n2     p1hat       p2hat
#> [1,] 0.009168536 0.08723163 0.1722089  0.95  5 56  0 29 0.0930101 0.005778473
#> 
#> $call
#>  distrib contrast    level     type      adj       cc       a1       b1 
#>    "bin"     "RD"   "0.95"   "jeff"  "FALSE"      "0"    "1.5"    "9.5" 
#>       a2       b2 
#>    "0.5"   "10.5"
```

WARNING: Note that confidence intervals for ratio measures should be
equivariant, in the sense that:

1.  interchanging the groups produces reciprocal intervals for RR or OR,
    i.e. the lower and upper confidence limits for $`x_1/n_1`$ vs
    $`x_2/n_2`$ are the reciprocal of the upper and lower limits for
    $`x_2/n_2`$ vs $`x_1/n_1`$
2.  in addition for OR, interchanging events and non-events produces
    reciprocal intervals, i.e. the lower and upper confidence limits for
    the odds ratio of $`(n_1 - x_1)/n_1`$ vs $`(n_2 - x_2)/n_2`$ are the
    reciprocal of the upper and lower limits for $`x_1/n_1`$ vs
    $`x_2/n_2`$.

Unfortunately, the MOVER method for contrast = “OR”, as described in
([Fagerland and Newcombe 2012](#ref-fagerland2012)), does not satisfy
the second of these requirements, unlike the asymptotic score and other
methods:

``` r

orci(x1 = 5, n1 = 56, x2 = 1, n2 = 29)$estimates
#> , , 5/56 vs 1/29
#> 
#>                                 lower      est    upper
#> SCAS                         0.385475 2.745098 56.08070
#> Gart                         0.389345 2.745098 56.59553
#> Miettinen-Nurminen           0.392931 2.745098 18.56785
#> Uncorrected Asymptotic Score 0.396590 2.745098 18.40530
#> MOVER-W                      0.363784 2.745098 18.78207
#> MOVER-J                      0.399469 2.745098 28.76882
#> Woolf logit                  0.305395 2.745098 24.67482
#> Gart adjusted logit          0.315100 2.745098 13.06682
```

``` r

1 / orci(x1 = 51, n1 = 56, x2 = 28, n2 = 29)$estimates
#> , , 51/56 vs 28/29
#> 
#>                                 lower      est     upper
#> SCAS                         56.08210 2.745096 0.3854748
#> Gart                         56.59630 2.745096 0.3893453
#> Miettinen-Nurminen           18.56769 2.745096 0.3929306
#> Uncorrected Asymptotic Score 18.40536 2.745096 0.3965898
#> MOVER-W                      20.71423 2.745096 0.4012105
#> MOVER-J                      32.91856 2.745096 0.4219566
#> Woolf logit                  24.67491 2.745096 0.3053949
#> Gart adjusted logit          13.06677 2.745096 0.3150998
```

## Technical details

([Miettinen and Nurminen 1985](#ref-miettinen1985)) took the concept of
the Wilson score interval for a single proportion and applied it to all
contrasts of two proportions (not just RD) and Poisson rates. The
procedure involves finding, at any possible value of the contrast
parameter $`\theta`$, the maximum likelihood estimates of the two
proportions, restricted to the given value of $`\theta`$. The resulting
estimates are then included in the score statistic, and the confidence
interval found by an iterative root-finding algorithm. Miettinen and
Nurminen’s statistic was defined to follow a $`\chi^2_1`$ distribution,
but the SCAS formula uses an equivalent statistic with a N(0, 1)
distribution, to facilitate one-sided tests.

Gart and Nam, in a series of papers between 1985 and 1990, proposed
similar confidence interval methods that included a skewness correction
in the score statistic (and also a bias correction for OR). Their
formulae omitted the ‘N-1’ bias correction used by Miettinen and
Nurminen.

Farrington and Manning published an article on non-inferiority tests for
RD and RR, using the same restricted maximum likelihood methodology as
Miettinen and Nurminen (but omitting the ‘N-1’ correction).

The SCAS method combines the features of all the above methods, to give
a comprehensive family of asymptotic score methods with optimised
performance for universal use for confidence intervals for a wide range
of sample sizes, with coherent tests for association or for
non-inferiority. For further details, see ([Laud 2017](#ref-laud2017)).

## References

Campbell, Ian. 2007. “Chi-Squared and Fisherirwin Tests of Two-by-Two
Tables with Small Sample Recommendations.” *Statistics in Medicine* 26
(19): 3661–75. <https://doi.org/10.1002/sim.2832>.

Fagerland, Morten W, Stian Lydersen, and Petter Laake. 2011.
“Recommended Confidence Intervals for Two Independent Binomial
Proportions.” *Statistical Methods in Medical Research* 24 (2): 224–54.
<https://doi.org/10.1177/0962280211415469>.

Fagerland, Morten W., and Robert G. Newcombe. 2012. “Confidence
Intervals for Odds Ratio and Relative Risk Based on the Inverse
Hyperbolic Sine Transformation.” *Statistics in Medicine* 32 (16):
2823–36. <https://doi.org/10.1002/sim.5714>.

Gart, John J., and Jun-mo Nam. 1988. “Approximate Interval Estimation of
the Ratio of Binomial Parameters: A Review and Corrections for
Skewness.” *Biometrics* 44 (2): 323. <https://doi.org/10.2307/2531848>.

Laud, Peter J. 2017. “Equal-Tailed Confidence Intervals for Comparison
of Rates.” *Pharmaceutical Statistics* 16 (5): 334–48.
<https://doi.org/10.1002/pst.1813>.

Miettinen, Olli, and Markku Nurminen. 1985. “Comparative Analysis of Two
Rates.” *Statistics in Medicine* 4 (2): 213–26.
<https://doi.org/10.1002/sim.4780040211>.

Newcombe, Robert G. 1998. “Interval Estimation for the Difference
Between Independent Proportions: Comparison of Eleven Methods.”
*Statistics in Medicine* 17 (8): 873–90.
<https://doi.org/10.1002/(sici)1097-0258(19980430)17:8%3C873::aid-sim779%3E3.0.co;2-i>.

Newcombe, Robert G. 2012. *Confidence Intervals for Proportions and
Related Measures of Effect Size*. CRC Press.
<https://doi.org/10.1201/b12670>.

[^1]: RNCP = the probability that the true value of $`\theta`$ falls
    outside to the right of the confidence interval. CP = the
    probability that the confidence interval contains the true value of
    $`\theta`$. Both calculated by precise enumeration of bivariate
    binomial probabilities.

[^2]: AN = “Wald” Approximate Normal, MOVER-J = Method of Variance
    Estimates Recovery, based on Newcombe’s method but using Jeffreys
    equal-tailed intervals instead of Wilson

[^3]: Newcombe considered a MOVER method based on equal-tailed Jeffreys
    intervals (“MOVER-R Jeffreys”) in chapter 11 of ([Newcombe
    2012](#ref-newcombe2012)). However, the example intervals shown in
    his Table 11.6 do not match the output of
    `moverci(2, 14, 1, 11, contrast = "RR", type = "jeff")`. The results
    for “MOVER-R Wilson” are matched by
    `moverci(2, 14, 1, 11, contrast = "RR", type = "wilson")`,
    suggesting that the discrepancy is in the calculation of the
    Jeffreys intervals themselves rather than the MOVER formula.
    However, Jeffreys intervals for a single proportion (“method 12”) in
    Newcombe’s Table 3.2 (chapter 3), are matched by
    `jeffreysci(1, 10)`, so the source of the discrepancy remains a
    mystery. Similar discrepancies are observed for contrast = “OR”.
