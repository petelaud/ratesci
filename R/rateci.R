#' Skewness-corrected asymptotic score ("SCAS") confidence intervals for
#' single binomial or Poisson rate using closed-form calculations.
#'
#' Closed-form function for computing confidence intervals for a single rate.
#' Note: For associated hypothesis tests, use `scoreci()` with `contrast = "p"`.
#' This function is vectorised in x, n.
#'
#' @param x Numeric vector of number of events.
#' @param n Numeric vector of sample sizes (for binomial rates) or exposure
#'   times (for Poisson rates).
#' @param distrib Character string indicating distribution assumed for the input
#'   data:\cr
#'   "bin" = binomial (default);\cr
#'   "poi" = Poisson.
#' @param level Number specifying confidence level (between 0 and 1, default
#'   0.95).
#' @param bcf Logical (default TRUE) indicating whether to apply bias correction
#'   in the score denominator. Applicable to `distrib = "bin"` only.
#' @param bign Sample size N to be used in the calculation of bcf, if different
#'   from n. (Used by transformed SCASp method for paired conditional OR in
#'   `pairbinci()`.)
#' @param xihat Number specifying estimated variance inflation factor for a
#'   skewness corrected version of the Saha Wilson Score interval for clustered
#'   binomial proportions. Need to calculate using BMS and WMS as per Saha 2016.
#'   Used by `clusterpci()` function for data entered per cluster.
#' @param cc Number or logical (default FALSE) specifying (amount of) continuity
#'   adjustment. Numeric value is taken as the gamma parameter in Laud 2017,
#'   Appendix S2 (default 0.5 for 'conventional' adjustment if cc = TRUE).
#' @param ... Other arguments.
#' @return A list containing the following components: \describe{
#'   \item{estimates}{a matrix containing estimated rate(s), the
#'   SCAS confidence interval, and the input values x and n.}
#'   \item{call}{details of the function call.} }
#' @author Pete Laud, \email{p.j.laud@@sheffield.ac.uk}
#' @references
#'   Laud PJ. Equal-tailed confidence intervals for comparison of
#'   rates. Pharmaceutical Statistics 2017; 16:334-348. (Appendix A.4)
#'
#' @export
scaspci <- function(x,
                    n,
                    distrib = "bin",
                    level = 0.95,
                    bcf = FALSE,
                    bign = n,
                    xihat = 1,
                    cc = FALSE,
                    ...) {
  if (as.character(cc) == "TRUE") cc <- 0.5
  z <- qnorm(1 - (1 - level) / 2)
    # bcf for the non-iterative formula
    if (distrib == "bin") {
      lambda <- switch(as.character(bcf),
                       "TRUE" = bign / (bign - 1),
                       "FALSE" = 1
      )
    } else lambda <- 1
    # Apply optional variance bias correction factor
    # and optional variance inflation factor for clustered data
    za <- qnorm(1 - (1 - level) / 2) * sqrt(lambda * xihat)

  if (distrib == "poi") {
    Du <- (x + cc) / n - (z^2 - 1) / (6 * n)
    #    Dl <- Rmpfr::pmax(0, (x - cc) / n - (z^2 - 1) / (6 * n))
    Dl <- pmax(0, (x - cc) / n - (z^2 - 1) / (6 * n))
    A <- 1
    Bu <- -2 * Du - z^2 / n
    Bl <- -2 * Dl - z^2 / n
    Cu <- Du^2
    Cl <- Dl^2
    D0 <- -1 / (6 * n) - x / n
    B0 <- 2 * D0
    A0 <- 1
    C0 <- D0^2
  } else if (distrib == "bin") {
    E <- (z^2 - 1) / (3 * n * lambda * xihat) - 1
    # Alteration to published formula,
    # to deal with non-nested intervals when level > 0.99
    #    Du <- Rmpfr::pmax(0, (n - x - cc) / n - (z^2 - 1) / (6 * n))
    #    Dl <- Rmpfr::pmax(0, (x - cc) / n - (z^2 - 1) / (6 * n))
    Du <- pmax(0, (n - x - cc) / n - (z^2 - 1) / (6 * n * lambda * xihat))
    Dl <- pmax(0, (x - cc) / n - (z^2 - 1) / (6 * n * lambda * xihat))
    A <- za^2 / n + E^2
    Bu <- 2 * E * Du - za^2 / n
    Bl <- 2 * E * Dl - za^2 / n
    Cu <- Du^2
    Cl <- Dl^2
    E0 <- 1 + 1 / (3 * n * lambda * xihat)
    D0 <- -1 / (6 * n * lambda * xihat) - x / n
    A0 <- E0^2
    B0 <- 2 * E0 * D0
    C0 <- D0^2
  }

  CI <- (cbind(
    #    Lower = Rmpfr::asNumeric((-Bl - sqrt(Rmpfr::pmax(0, Bl^2 - 4 * A * Cl))) / (2 * A)),
    #    MLE = Rmpfr::asNumeric((-B0 - sqrt(Rmpfr::pmax(0, (B0^2 - 4 * A0 * C0)))) / (2 * A0)),
    lower = ifelse(x == 0, 0, ((-Bl - sqrt(pmax(0, Bl^2 - 4 * A * Cl))) / (2 * A))),
    est = ((-B0 - sqrt(pmax(0, (B0^2 - 4 * A0 * C0)))) / (2 * A0)),
    upper = if (distrib == "bin") {
      #      Rmpfr::asNumeric(1 - (-Bu - sqrt(Rmpfr::pmax(0, Bu^2 - 4 * A * Cu))) / (2 * A))
      ifelse(x == n, 1, (1 - (-Bu - sqrt(pmax(0, Bu^2 - 4 * A * Cu))) / (2 * A)))
    } else {
      #      Rmpfr::asNumeric((-Bu + sqrt(Rmpfr::pmax(0, Bu^2 - 4 * A * Cu))) / (2 * A))
      ((-Bu + sqrt(pmax(0, Bu^2 - 4 * A * Cu))) / (2 * A))
    },
    x = x,
    n = n
  ))
  outlist <- list(estimates = CI)
  call <- c(
    distrib = distrib, level = level, bcf = bcf,
    cc = cc
  )
  outlist <- append(outlist, list(call = call))
  return(outlist)
}

#' Selected confidence intervals for the single binomial or Poisson rate.
#'
#' Confidence intervals for the single binomial or Poisson rate. Including
#' SCAS or Jeffreys intervals, with or without continuity adjustment, and
#' 'exact' Clopper-Pearson/Garwood or mid-p intervals.
#' This function is vectorised in x, n.
#'
#' @param x Numeric vector of number of events.
#' @param n Numeric vector of sample size (for binomial rate) or exposure
#'   times (for Poisson rate).
#' @param distrib Character string indicating distribution assumed for the input
#'   data: "bin" = binomial (default), "poi" = Poisson.
#' @param level Number specifying confidence level (between 0 and 1, default
#'   0.95).
#' @param cc Number or logical (default FALSE) specifying continuity
#'   adjustment.
#' @return A list containing, for each method, a matrix containing lower and upper
#'   confidence limits and point estimate of p for each value of x and n.
#'   Methods shown depend on the cc
#'   parameter, which specifies whether the continuity adjustment is applied to
#'   the SCAS and Jeffreys methods. The corresponding 'exact' method is
#'   Clopper-Pearson/Garwood if cc = TRUE and mid-p if cc = FALSE.
#'   The last list item contains details of the function call.
#' @author Pete Laud, \email{p.j.laud@@sheffield.ac.uk}
#' @references
#'   Laud PJ. Equal-tailed confidence intervals for comparison of
#'   rates. Pharmaceutical Statistics 2017; 16:334-348. (Appendix A.4)
#'
#'   Brown LD, Cai TT and DasGupta A. Interval estimation for a binomial
#'   proportion. Statistical Science 2001; 16(2):101-133.
#'
#' @export
rateci <- function(x,
                   n,
                   distrib = "bin",
                   level = 0.95,
                   cc = FALSE
                   ) {
  # in case x is input as a vector but n is not
  if (length(n) < length(x) && length(x) > 1) {
    n <- rep(n, length.out = length(x))
  }
  if (as.character(cc) == "TRUE") cc <- 0.5
  if (!is.numeric(c(x, n, level))) {
    print("Non-numeric inputs!")
    stop()
  }
  if (distrib == "bin" && any(x > n + 0.001)) {
    print("x > n not possible for distrib = 'bin'")
    stop()
  }
  if (any(c(x, n) < 0)) {
    print("Negative inputs!")
    stop()
  }
  if (any(n == 0)) {
    print("Sample size is zero!")
    stop()
  }

  ci_scas <- scaspci(
    x = x,
    n = n,
    distrib = distrib,
    level = level,
    cc = cc
  )$estimates[, c(1:5), drop = FALSE]
  ci_jeff <- jeffreysci(
    x = x,
    n = n,
    ai = 0.5,
    bi = 0.5,
    cc = cc,
    level = level,
    distrib = distrib,
    adj = TRUE
  )$estimates[, c(1:5), drop = FALSE]
  ci_exact <- exactci(
    x = x,
    n = n,
    level = level,
    midp = 0.5 - cc,
    distrib = distrib
  )
  if (cc == 0) {
    outlist <- list(scas = ci_scas, jeff = ci_jeff, midp = ci_exact)
  } else if (cc == 0.5) {
    if (distrib == "bin") {
      outlist <- list(scas_cc = ci_scas, jeff_cc = ci_jeff, cp = ci_exact)
    } else {
      outlist <- list(scas_cc = ci_scas, jeff_cc = ci_jeff, garwood = ci_exact)
    }
  } else {
    outlist <- list(scas_cc = ci_scas, jeff_cc = ci_jeff)
    # exact method not applicable if using a compromise value of cc
  }
  call <- c(
    distrib = distrib, level = level, cc = cc
  )
  outlist <- append(outlist, list(call = call))
  return(outlist)
}

#' Clopper-Pearson/Garwood and mid-p intervals for single binomial or
#' Poisson rate
#'
#' to calculate exact 'exact' confidence interval for a single binomial or Poisson rate x/n
#'
#' @author Pete Laud, \email{p.j.laud@@sheffield.ac.uk}
#'
#' @noRd
exactci <- function(x,
                    n,
                    level = 0.95,
                    midp = TRUE,
                    distrib = "bin",
                    precis = 8
                    ) {
  alpha <- 1 - level
  if (as.character(midp) == "TRUE") midp <- 0.5
  if (distrib == "bin") {
    lowroot <- function(p) {
      pbinom(x - 1, n, p) + midp * dbinom(x, n, p) -
        (1 - alpha / 2)
    }
    midroot <- function(p) {
      pbinom(x - 1, n, p) + 0.5 * dbinom(x, n, p) - 0.5
    }
    uproot <- function(p) {
      pbinom(x, n, p) - midp * dbinom(x, n, p) - alpha / 2
    }
  } else if (distrib == "poi") {
    lowroot <- function(p) ppois(x - 1, p) + midp * dpois(x, p) - (1 - alpha / 2)
    midroot <- function(p) ppois(x - 1, p) + 0.5 * dpois(x, p) - 0.5
    uproot <- function(p) ppois(x, p) - midp * dpois(x, p) - alpha / 2
  }
  lower <- bisect(
    ftn = lowroot, precis = precis, uplow = "low",
    contrast = "p", distrib = distrib
  )
  est <- bisect(
    ftn = midroot, precis = precis, uplow = "low",
    contrast = "p", distrib = distrib
  )
  upper <- bisect(
    ftn = uproot, precis = precis, uplow = "up", contrast = "p",
    distrib = distrib
  )
  outdata <- cbind(lower = lower, est = est, upper = upper) /
    ifelse(distrib == "poi", n, 1)
  return(cbind(outdata, x = x, n = n))
}


