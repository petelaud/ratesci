% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/moverci.R
\name{moverbci}
\alias{moverbci}
\title{Approximate Bayesian ("MOVER-B") confidence intervals for
comparisons of independent binomial or Poisson rates.}
\usage{
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
}
\arguments{
\item{x1, x2}{Numeric vectors of numbers of events in group 1 & group 2
respectively.}

\item{n1, n2}{Numeric vectors of sample sizes (for binomial rates) or exposure
times (for Poisson rates) in each group.}

\item{a1, b1, a2, b2}{Numbers defining the Beta(ai,bi) prior distributions for
each group (default ai = bi = 0.5 for Jeffreys uninformative priors). Gamma
priors for Poisson rates require only a1, a2.}

\item{distrib}{Character string indicating distribution assumed for the input
data:\cr
"bin" = binomial (default);\cr
"poi" = Poisson.}

\item{contrast}{Character string indicating the contrast of interest: \cr
"RD" = rate difference (default); \cr
"RR" = rate ratio; \cr
"OR" = odds ratio; \cr
"p" gives an interval for the single proportion \code{x1/n1}.}

\item{level}{Number specifying confidence level (between 0 and 1, default
0.95).}

\item{cc}{Number or logical specifying (amount of) continuity adjustment
(default FALSE). Numeric value is taken as the gamma parameter in Laud
2017, Appendix S2 (default 0.5 if \code{cc = TRUE}). Forced equal to 0.5 if
\code{type = "exact"}.}

\item{...}{Additional arguments.}
}
\value{
A list containing the following components: \describe{
\item{estimates}{a matrix containing estimates of the rates in each group
and of the requested contrast, with its confidence interval}
\item{call}{details of the function call} }
}
\description{
Wrapper function for the MOVER-B methods.  Approximate Bayesian confidence
intervals for the rate (or risk) difference ("RD") or ratio ("RR") for
independent binomial or Poisson rates, or for odds ratio ("OR", binomial
only). (developed from Newcombe, Donner & Zou, Li et al, and Fagerland &
Newcombe, and generalised as "MOVER-B" in Laud 2017) including
special case "MOVER-J" using non-informative priors with optional continuity
adjustment.  This function is vectorised in x1, x2, n1, and n2.
}
\author{
Pete Laud, \email{p.j.laud@sheffield.ac.uk}
}
