# Approximate Bayesian ("MOVER-B") confidence intervals for comparisons of independent binomial or Poisson rates.

Wrapper function for the MOVER-B methods. Approximate Bayesian
confidence intervals for the rate (or risk) difference ("RD") or ratio
("RR") for independent binomial or Poisson rates, or for odds ratio
("OR", binomial only). (developed from Newcombe, Donner & Zou, Li et al,
and Fagerland & Newcombe, and generalised as "MOVER-B" in Laud 2017)
including special case "MOVER-J" using non-informative priors with
optional continuity adjustment. This function is vectorised in x1, x2,
n1, and n2.

## Usage

``` r
moverbci(
  x1,
  n1,
  x2,
  n2,
  a1 = 0.5,
  b1 = 0.5,
  a2 = 0.5,
  b2 = 0.5,
  distrib = "bin",
  contrast = "RD",
  level = 0.95,
  cc = 0,
  ...
)
```

## Arguments

- x1, x2:

  Numeric vectors of numbers of events in group 1 & group 2
  respectively.

- n1, n2:

  Numeric vectors of sample sizes (for binomial rates) or exposure times
  (for Poisson rates) in each group.

- a1, b1, a2, b2:

  Numbers defining the Beta(ai,bi) prior distributions for each group
  (default ai = bi = 0.5 for Jeffreys uninformative priors). Gamma
  priors for Poisson rates require only a1, a2.

- distrib:

  Character string indicating distribution assumed for the input data:  
  "bin" = binomial (default);  
  "poi" = Poisson.

- contrast:

  Character string indicating the contrast of interest:  
  "RD" = rate difference (default);  
  "RR" = rate ratio;  
  "OR" = odds ratio;  
  "p" gives an interval for the single proportion `x1/n1`.

- level:

  Number specifying confidence level (between 0 and 1, default 0.95).

- cc:

  Number or logical specifying (amount of) continuity adjustment
  (default FALSE). Numeric value is taken as the gamma parameter in Laud
  2017, Appendix S2 (default 0.5 if `cc = TRUE`). Forced equal to 0.5 if
  `type = "exact"`.

- ...:

  Additional arguments.

## Value

A list containing the following components:

- estimates:

  a matrix containing estimates of the rates in each group and of the
  requested contrast, with its confidence interval

- call:

  details of the function call

## Author

Pete Laud, <p.j.laud@sheffield.ac.uk>
