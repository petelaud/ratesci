
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ratesci <a href="https://petelaud.github.io/ratesci/"><img src="man/figures/logo.png" alt="ratesci website" align="right" height="139"/></a>

<!-- badges: start -->

[![R-CMD-check](https://github.com/petelaud/ratesci/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/petelaud/ratesci/actions/workflows/R-CMD-check.yaml)
[![License: GPL
v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![CRAN_Version](https://www.r-pkg.org/badges/version/ratesci)](https://cran.r-project.org/package=ratesci)
[![Total
Downloads](https://cranlogs.r-pkg.org/badges/grand-total/ratesci)](https://cranlogs.r-pkg.org/badges/grand-total/ratesci)

<!-- badges: end -->

ratesci is an [R](https://www.r-project.org) package to compute
confidence intervals for rate (or risk) difference (‘RD’), rate ratio
(‘RR’, also known as relative risk), or odds ratio (‘OR’). All three
contrasts apply for binomial proportions, and the first two may also be
used for the comparison of Poisson ‘exposure-adjusted’ incidence rates.
`scoreci()` incorporates ‘skewness-corrected’ asymptotic score (‘SCAS’)
methods, which ensure equal-tailed coverage (or central location), in
other words for a nominal 95% confidence interval, the one-sided
non-coverage probability is (on average) close to 2.5% on each side.
Stratified calculations are also catered for (e.g. meta-analysis,
including random effects), as well as confidence intervals for the
single binomial or Poisson rate, and for clustered binomial proportion
(with `clusterpci()`) and binomial matched pairs (with `pairbinci()`).
Corresponding hypothesis tests against any specified null parameter
value are provided in each case. Omission of the skewness correction is
also allowed, resulting in the often-recommended ‘Miettinen-Nurminen’
asymptotic score methods, which can have inferior one-sided coverage,
especially for RR. The hypothesis test for binomial RD or RR when the
skewness correction is omitted corresponds to the Farrington-Manning
test.

The stratified (fixed effects) version without skewness correction
produces a hypothesis test which is equivalent to the
Cochran-Mantel-Haenszel (CMH) test, when MH weighting is used for RD or
RR, or IVS weighting for OR, and the corresponding confidence intervals
are guaranteed to be coherent with the test. In the single-stratum case,
the hypothesis tests are equivalent to a Chi-squared test, and for
paired proportions, a McNemar test (in both cases with extension to null
hypotheses for equivalence/non-inferiority tests).

For large (single-stratum) sample sizes, the ‘MOVER’ methods
(`moverci()`) improve on traditional approximate methods with respect to
one-sided and two-sided coverage, particularly in the case of RR. Being
based on Bayesian methods, these also allow the option to incorporate
prior beliefs about the rates in each group - by default, the
‘non-informative’ Jeffreys priors are used. These methods are adapted
from the Newcombe ‘square-and-add’ method, which is also included for
reference.

For those wishing to achieve strictly conservative coverage, so-called
‘continuity corrections’ are provided as approximations to ‘exact’
methods, with the option to adjust the strength of the correction. The
performance of these adjustments has not been extensively evaluated, but
they appear to be more successful for SCAS than for MOVER, in terms of
achieving conservative coverage.

An online calculator based on this package is available
[here](https://ssu.shef.ac.uk/ratesci/calc.php)

## Installation

The current official
(i.e. [CRAN](https://CRAN.R-project.org/package=ratesci)) release can be
installed within R with:

``` r
install.packages("ratesci")
```

The latest development version of the package can be installed with:

``` r
# install.packages("pak")
pak::pak("petelaud/ratesci")
```

This builds the package from source based on the current version on
[GitHub](https://github.com/petelaud/ratesci)

## Example

Below is a basic example which shows you how to request a confidence
interval for the difference between proportions 5/56 - 0/29. The `$call`
output element shows that the default settings give an interval for the
risk difference (`contrast` = “RD”), for binomial proportions (`distrib`
= “bin”), at a 95% confidence level. Variance bias correction (`bcf`)
and skewness correction (`skew`) are applied, continuity correction
(`cc`) is not. This is the skewness-corrected asymptotic score (“SCAS”)
confidence interval. (For Miettinen-Nurminen, use `skew` = FALSE, for
Gart-Nam, use `bcf` = FALSE.)

``` r
library(ratesci)
scoreci(x1 = 5, n1 = 56, x2 = 0, n2 = 29, precis = 4)
#> $estimates
#>            Lower        MLE    Upper level x1 n1 x2 n2      p1hat p2hat
#> [1,] -0.01861954 0.09168625 0.186718  0.95  5 56  0 29 0.08928571     0
#>           p1mle p2mle
#> [1,] 0.09168625     0
#> 
#> $pval
#>         chisq pval2sided theta0 scorenull pval_left pval_right
#> [1,] 3.024862 0.08199729      0  1.739213 0.9590014 0.04099865
#> 
#> $call
#>  distrib contrast    level      bcf     skew       cc 
#>    "bin"     "RD"   "0.95"   "TRUE"   "TRUE"  "FALSE"
```

An example of a paired analysis follows, using the data from Table II of
Fagerland et al. Here the bias and skewness corrections are again
applied by default. Omitting both would produce the Tango asymptotic
score interval for `contrast = RD`, or the Tang method for
`contrast = RR`.

``` r
pairbinci(x = c(1, 1, 7, 12), precis = 4)
#> $data
#>    x2i
#> x1i  0  1
#>   0 12  7
#>   1  1  1
#> 
#> $estimates
#>        Lower     MLE   Upper level  p1hat p2hat  p1mle  p2mle phi_hat phi_c
#> [1,] -0.5281 -0.2859 -0.0184  0.95 0.0952 0.381 0.0952 0.3811  0.0795     0
#>      psi_hat
#> [1,]  1.7143
#> 
#> $pval
#>         chisq pval2sided theta0 scorenull  pval_left pval_right
#> [1,] 4.285714 0.03843393      0 -2.070197 0.01921697   0.980783
#> 
#> $call
#> contrast   method    level      bcf     skew       cc 
#>     "RD"  "Score"   "0.95"   "TRUE"   "TRUE"  "FALSE"
```

#### Overview

ratesci contains the following functions:

For comparisons of rates (contrasts RD, RR and OR):

- `scoreci()`: for score-based confidence intervals including SCAS,
  Miettinen-Nurminen and Gart-Nam, with or without stratification.
- `scasci()`: wrapper function to compute SCAS intervals.
- `tdasci()`: wrapper function to compute TDAS stratified intervals
  incorporating random effects.
- `moverci()`: for the MOVER methods, including Newcombe and MOVER-J.
- `moverbci()`: wrapper function to compute MOVER-B intervals.
- `pairbinci()`: for paired binomial data, including SCAS, asymptotic
  score and MOVER methods for RD and RR, and transformed binomial
  intervals for conditional OR.

For single binomial or Poisson rates:

- `scaspci()`: non-iterative SCAS method for a single rate. For
  stratified calculations use `scoreci()` with contrast = “p”.
- `jeffreysci()`: wrapper function to compute Jeffreys interval for a
  single rate (with option to incorporate prior information).
- `rateci()`: wrapper function for selected methods for a single rate,
  including SCAS, Jeffreys, midp and Clopper-Pearson/Garwood.
- `clusterpci()`: Saha’s Wilson-based interval for a single proportion
  based on clustered data, with a skewness-corrected version.
