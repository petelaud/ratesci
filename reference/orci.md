# Confidence intervals for odds ratio.

Confidence intervals for comparisons of two binomial rates. This
convenience wrapper function produces a selection of the methods below
as appropriate for the odds ratio contrast (OR), with or without
continuity adjustment.

- SCAS (skewness-corrected asymptotic score)

- Miettinen-Nurminen, Gart and unadjusted Asymptotic Score methods

- MOVER-R Wilson

- MOVER-R Jeffreys

- Approximate normal (Woolf/Gart) methods (strongly advise this is not
  used for any purpose but included for reference)

## Usage

``` r
orci(x1, n1, x2, n2, level = 0.95, std_est = TRUE, cc = FALSE, precis = 8)
```

## Arguments

- x1, x2:

  Numeric vectors of numbers of events in group 1 & group 2
  respectively.

- n1, n2:

  Numeric vectors of sample sizes (for binomial rates) or exposure times
  (for Poisson rates) in each group.

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

  Number (default 8) specifying precision (i.e. number of decimal
  places) to be used in root-finding subroutine for the score confidence
  intervals. (Note other methods use closed-form calculations so are not
  affected.)

## References

Laud PJ. Equal-tailed confidence intervals for comparison of rates.
Pharmaceutical Statistics 2017; 16:334-348. (Appendix A.4)

## Author

Pete Laud, <p.j.laud@sheffield.ac.uk>
