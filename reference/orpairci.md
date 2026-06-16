# Confidence intervals for conditional odds ratio (OR) with paired binomial rates.

Confidence intervals for comparisons of two binomial rates from paired
data. This convenience wrapper function produces a selection of the
methods below for the conditional odds ratio (OR) contrast, with or
without optional continuity adjustment (where available).

- Transformed SCAS (skewness-corrected asymptotic score)

- Transformed Wilson Score method

- Transformed mid-P

- Transformed Jeffreys

- Approximate log-normal (Wald) method

## Usage

``` r
orpairci(x, level = 0.95, std_est = TRUE, cc = FALSE, precis = 8)
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

  Number (default 6) specifying precision (i.e. number of decimal
  places) to be used in output.

## References

Laud PJ. Improved confidence intervals and tests for paired binomial
proportions. (2026, Under review)

## Author

Pete Laud, <p.j.laud@sheffield.ac.uk>
