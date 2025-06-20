# ratesci 1.0.0 (2025-06-20)

## New features
* New function `clusterpci()` for CI and test for a single binomial proportion 
  from clustered data.
* Example datasets are now included.
* Improved documentation with pkgdown website & vignettes.

### In `pairbinci()`:
* `skew` for skewness correction. 
* `bcf` for variance bias correction.
* Default paired RD and RR method changed to new SCAS method 
  (i.e. including skewness correction - manuscript under review).
* `method_RD`, `method_RR` and `method_OR` are replaced with `method`.
* Bonett-Price methods for RD and RR (including proposed Jeffreys variant 
  option for RR).
* TDAS method is deprecated.
* Default for MOVER method changed to Jeffreys.
* MOVER calculations now use x/N as point estimate instead of median from the
  Beta distribution.
* `cc` uses a new form of correction for RR giving equivariant intervals. Also 
  allows consistency with the continuity-corrected McNemar test (or an 
  intermediate correction of the user's choosing). `cctype` is deprecated.
* Default conditional odds ratio method changed to SCASp (with closed-form 
  calculation).
* Output object now includes estimates of p1, p2, phi (correlation) 
  and psi (odds ratio used by Fagerland et al).
* Output object now includes function call.

### In `scaspci()`:
* `bcf` option now implemented for contrast = "p" (default = FALSE).
* `bign` allows a different sample size to be used in the bias correction 
  (used within transformed SCASp method for paired OR in `pairbinci`, for 
  consistency with 'N-1' test).

### In `scoreci()`:
* `bcf` option now implemented for contrast = "p" (default = FALSE).
* (Note adjusted sample size for bias correction can be achieved by including a 
   non-zero value for n2.)
* ORbias, RRtang, and MNtol arguments renamed as or_bias, rr_tang and mn_tol.
* Implementation of `precis` argument is improved for RR and OR contrasts.
* For contrast = "RD", weighting = "Tang" provides optimal test if RR is 
  constant across strata.

## Bug fixes
### In `exactci()`:
* Corrected duplicate point estimate reported for vector inputs.
* Derive point estimates to match LCL and UCL with level = 0.
* Corrected LCL for Poisson mid-p method.

## Other
* Output object column names are updated (lower, est, upper) for consistency 
  and style conformity.
* Tests added to confirm consistency of score methods vs McNemar test.
* Dependence on polynom package removed.
* Edition 3 of testthat implemented.

# ratesci 0.5.0 (2025-01-10)

## New features
### In `pairbinci()`:
* `cc` continuity correction is now available for all methods for all contrasts. 
* `cctype` controls the type of correction to apply for `contrast` = "RR".
* New default `method_RD` = "Score_closed" for non-iterative calculation of 
  the Tango score interval for `contrast` = "RD". Thanks to Tony Yang for 
  permission to use the code in his 2013 paper.
* New default `method_RR` = "Score_closed" for non-iterative calculation of 
  the Tang score interval for `contrast` = "RR". Thanks to Guogen Shan for 
  contributing code via email.
* Added paired MOVER methods with `method_RD` = "MOVER" and `method_RR` = "MOVER".
  Also "MOVER_newc" incorporates Newcombe's correlation correction.
* Added `moverbase`, for specifying different versions of the MOVER methods 
  (Wilson, Jeffreys, midp or SCAS).
* Added "jeff" and "wilson" `method_OR` options for transformed binomial 
  methods for OR.
* Confirmed and documented that the 2-sided significance test is equivalent 
  to the McNemar test (with or without continuity correction).

### In `scoreci()`:
* Confirmed that continuity corrections for all stratified (fixed-effects) 
  binomial contrasts are consistent with the Mantel-Haenszel correction.
* Updated heterogeneity test to consistently omit non-informative 
  (but non-empty) strata, and output the degrees of freedom.

### In `moverci()`:
* Added continuity correction for `type` = "wilson".
* Added options for `type` = "SCAS" and "midp" intervals.
* Standardised output to include lower CL, midpoint, upper CL, in that order.

## Bug fixes
### In `scoreci()`:
* Improved handling of special cases for MN weighting (#25, thanks to 
  Vincent Jaquet for reporting the issue and proposed solution. 
  Also #27 for RR, thanks to @lovestat.) 
  As a result, double-zero strata need not be excluded when weighting = "MN".

### In `moverci()`:
* Corrected calculation of score intervals for single Poisson rate, using 
  Rao score interval.
* Same correction affects MOVER method for comparison of Poisson rates
 [i.e. `moverci()` with `distrib` = "poi" and `type` = "wilson"]

## Other
* Improved documentation of hypothesis tests and continuity corrections, clarifying links to Chi-squared tests and CMH test with selected weights.
* Correction to documentation of default weights for OR.
* Added tests confirming equivalence of iterative and closed-form methods in pairbinci.

# ratesci 0.4-0 (2021-12-04)

## New features
### In `scoreci()`:
* MN weighting now iterates to convergence (@jonjvallejo, #20).
* Added optional prediction interval for random effects method (also in `tdasci()`).
* Added xlim and ylim arguments to control plot output.
* Added sda & fda arguments for optional sparse/full data adjustment 
    when x1 + x2 = 0 or x1 + x2 = n1 + n2 in a stratum.
* Added INV option for weights that omit the variance bias correction.
* Added RRtang argument to apply Tang's alternative score for RR (recommended 
for stratified analysis with INV/IVS weights. Experimental for Poisson RR).
    `Stheta = (p1hat - p2hat * theta) / p2d`  (see Tang 2020)
* Added simplified skewness correction option (causes p-values to be omitted, see Tang 2021 & Laud 2021).
* Introduced warning and plot features for very rare occasions when quadratic 
  skewness correction cannot be calculated due to a negative discriminant.
* p-value suppressed where affected by negative discriminants.
* Changed ORbias default to TRUE (see Laud 2018).
* Changed weighting default to MH for RD & RR, INV for OR (for consistency with CMH test).
* Added hetplot argument to separate heterogeneity plots from score function plot.
* Uninformative strata are now retained in the analysis except if: 
  * contrast = OR with MH weighting;
  * contrast = RR with IVS/INV weighting if RRtang = FALSE;
  * random = TRUE (needs further evaluation);
  * excluded using new option dropzeros = TRUE.

### In `tdasci()`:
* Default uses skew = TRUE for stratum CIs.

## Bug fixes
* MN weighting in `scoreci()` corrected for distrib="poi".
* Fixed bug in `scoreci()` for calculation of stratum CIs with random=TRUE.
* Fixed bug in `scoreci()` for distrib = "poi" and contrast = "p" (#7).
* Fixed finite precision bug in `scaspci()`.
* Fixed bug in `rateci()` for closed-form calculation of continuity-corrected SCAS.
* Fixed bug in `scoreci()` for stratified zero scores calculated as NA, resulting in UL = 0. (Thanks to Lidia Mukina for reporting the bug.)
* Fixed variable plot ranges for vectorised inputs.

## Other
* Renamed tdas argument to 'random'.
* Removed redundant t2 variable.

# ratesci 0.3-0 (2018-02-15)

## New features
* Added bias correction in `scoreci()` for OR SCAS method (derived from Gart 1985).
* Added score methods (Tango & Tang) as default for paired binomial RD and RR in `pairbinci()`.
* Added transformed mid-p method for paired OR for comparison with transformed SCAS.
* Added `scaspci()` for non-iterative SCAS methods for single binomial or Poisson rate.
* Added `rateci()` for selected methods for single binomial or Poisson rate.

## Bug fixes
* Fixed bug in `pairbinci()` for contrast="OR".
* Fixed bug in `moverci()` for contrast="p" and type="wilson".
* Corrected error in cc for stratified SCAS method for OR.
* Clarified documentation regarding continuity corrections.
* Set Stheta to 0 if |Stheta|<cc in `scoreci()`
* Fixed stratified calulations for contrast = "p" in `scoreci()`.

# ratesci 0.2-0 (2017-04-21)

## New features
* Added `pairbinci()` for all comparisons of paired binomial rates.
* Added option to suppress warnings in scoreci.
* Added Galbraith plot (for assessing stratum heterogeneity) to `scoreci()`.
* Added qualitative interaction test to `scoreci()`.
* Added stratum estimates & CIs to `scoreci()` output when stratified = TRUE.

## Bug fixes
* Fixed bug for contrast = "p" in `moverci()`.
* Fixed bug in `tdasci()` wrapper function.
* Fixed bug for stratified OR.
* Altered adjustment options for boundary cases in `moverci()`.
* Changed point estimate used in `moverci()` to posterior median for type = "jeff",
  to ensure consistent calculations with informative priors.
