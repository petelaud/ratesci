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
      ifelse(n == 0, Inf, (-Bu + sqrt(pmax(0, Bu^2 - 4 * A * Cu))) / (2 * A))
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
#' @description
#' Confidence intervals for the single binomial or Poisson rate. This
#' convenience wrapper function produces a selection of alternative methods. The
#' first three are recommended for achieving 1-sided and 2-sided coverage
#' probability close to the nominal levels (see Laud 2017 and Laud 2018):
#' - SCAS (skewness-corrected asymptotic score)
#' - Jeffreys
#' - mid-p (two versions, using exact calculation or approximation via
#' Beta/Gamma distribution, see p.115 of Brown et al.))
#' The following more approximate methods are included for users wishing to use
#' a more established or commonly used method:
#' - Wilson score
#' - Agresti-Coull
#' - Wald (strongly advise this is not used for any purpose but included for reference)
#'
#' All methods can be made more conservative with a 'continuity adjustment',
#' which may either be specified as TRUE, or an intermediate 'compromise'
#' value between 0 and 0.5 may be selected. When `cc` is `TRUE` or `0.5`, the
#' mid-p method becomes the Clopper-Pearson interval (or Garwood for Poisson
#' rates).
#' Note that Brown et al's Beta formulation perfectly matches the exact interval
#' when `cc` is TRUE (i.e. for Clopper-Pearson) but not when `cc` is `FALSE`
#' (for mid-p)
#' All methods except Agresti-Coull have equivalent formulae for the Poisson
#' distribution: Garwood for Clopper-Pearson, Rao score for Wilson score.
#' Jeffreys has a Poisson equivalent using the Gamma distribution. e.g. See Brown
#' et al. 2003, Swift 2009 and Laud 2017. The formulation for the approximate
#' mid-p interval using Gamma distribution for a Poisson rate has been deduced
#' by the package author from the corresponding formulae from Brown et al., and
#' has not (to the best of my knowledge) been published.
#'
#' This function is vectorised in x, n.
#'
#' @param x Numeric vector of number of events.
#' @param n Numeric vector of sample size (for binomial rate) or exposure
#'   times (for Poisson rate).
#' @param distrib Character string indicating distribution assumed for the input
#'   data: "bin" = binomial (default), "poi" = Poisson.
#' @param level Number specifying confidence level (between 0 and 1, default
#'   0.95).
#' @param std_est logical, specifying if the crude point estimate for the
#'   proportion value x/n should be returned (TRUE, default) or the
#'   method-specific alternative point estimate consistent with a 0% confidence
#'   interval (FALSE).
#' @param cc Number or logical (default FALSE) specifying continuity
#'   adjustment.
#' @param precis Number (default 8) specifying precision (i.e. number of decimal
#'   places) to be used in root-finding subroutine for the exact confidence
#'   interval. (Note all other methods use closed-form calculations so are not
#'   affected.)
#' @return A list containing, for each method, a matrix containing lower and
#'   upper confidence limits and point estimate of p for each value of x and n.
#'   Methods shown depend on the cc parameter, which specifies whether the
#'   continuity adjustment is applied to the SCAS and Jeffreys methods. The
#'   corresponding 'exact' method is Clopper-Pearson/Garwood if cc = TRUE and
#'   mid-p if cc = FALSE.
#'   The last list item contains details of the function call.
#' @author Pete Laud, \email{p.j.laud@@sheffield.ac.uk}
#' @references
#'   Laud PJ. Equal-tailed confidence intervals for comparison of
#'   rates. Pharmaceutical Statistics 2017; 16:334-348. (Appendix A.4)
#'
#'   Brown LD, Cai TT and DasGupta A. Interval estimation for a binomial
#'   proportion. Statistical Science 2001; 16(2):101-133.
#'
#'   Garwood F. Fiducial limits for the Poisson distribution. Biometrika
#'   1936; 28(3-4):437, doi:10.1093/biomet/28.3-4.437.
#'
#' @export
rateci <- function(x,
                   n,
                   distrib = "bin",
                   level = 0.95,
                   std_est = TRUE,
                   cc = FALSE,
                   precis = 8
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
  if (distrib == "poi" && any(n == 0 & x > n + 0.001)) {
    print("x > 0 not possible for distrib = 'poi' if n = 0")
    stop()
  }
  if (any(c(x, n) < 0)) {
    print("Negative inputs!")
    stop()
  }

  ci_scas <- scaspci(
    x = x,
    n = n,
    distrib = distrib,
    level = level,
    cc = cc
  )$estimates[, c(1:5), drop = FALSE]
  if (std_est) ci_scas[, 2] <- x / n

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
  if (std_est) ci_jeff[, 2] <- x / n

  ci_exact <- exactci(
    x = x,
    n = n,
    level = level,
    midp = 0.5 - cc,
    beta = FALSE,
    distrib = distrib,
    precis = precis
  )
  if (std_est) ci_exact[, 2] <- x / n

  ci_beta <- exactci(
    x = x,
    n = n,
    level = level,
    midp = 0.5 - cc,
    beta = TRUE,
    distrib = distrib
  )
  if (std_est) ci_beta[, 2] <- x / n

  ci_wilson <- cbind(wilsonci(
    x = x,
    n = n,
    level = level,
    cc = cc,
    distrib = distrib
  ), x, n)

  ci_wald <- cbind(waldci(
    x = x,
    n = n,
    level = level,
    cc = cc,
    distrib = distrib
  ), x, n)


  z <- qnorm(1 - (1 - level)/2)
  ci_ac <-  cbind(waldci(
    x = x + z^2 / 2,
    n = n + z^2,
    level = level,
    cc = cc,
    distrib = distrib
  ), x, n)
  if (std_est) ci_ac[, 2] <- x / n

#  outarr <- 1
  mydimnames <- dimnames(ci_scas)
  if (distrib == "bin") {
    outarr <- array(c(ci_scas, ci_jeff, ci_exact, ci_beta,
                      ci_wilson, ci_wald, ci_ac),
                  dim <- c(dim(ci_scas), 7))
  } else if (distrib == "poi") {
    outarr <- array(c(ci_scas, ci_jeff, ci_exact, ci_beta,
                      ci_wilson, ci_wald),
                    dim <- c(dim(ci_scas), 6))
  }
  if (cc == 0) {
    mydimnames[[3]] <- c("SCAS", "Jeffreys", "mid-p", "mid-p(beta)",
                         "Wilson", "Wald", "Agresti-Coull")
    if (distrib == "bin") {
      outlist <- list(scas = ci_scas, jeff = ci_jeff, midp = ci_exact,
                      midp_beta = ci_beta, wilson = ci_wilson, wald = ci_wald)
    } else if (distrib == "poi") {
      outlist <- list(scas = ci_scas, jeff = ci_jeff, midp = ci_exact,
                      midp_gamma = ci_beta)
      mydimnames[[3]][4] <- "mid-p(gamma)"
    }
  } else if (cc == 0.5) {
    mydimnames[[3]] <- c("SCAS_cc", "Jeffreys_cc", "Clopper-Pearson", "CP(beta)",
                         "Wilson_cc", "Wald_cc", "Agresti-Coull_cc")
    if (distrib == "bin") {
      outlist <- list(scas_cc = ci_scas, jeff_cc = ci_jeff, cp = ci_exact,
                      cp_beta = ci_beta)
    } else if (distrib == "poi") {
      outlist <- list(scas_cc = ci_scas, jeff_cc = ci_jeff, garwood = ci_exact,
                      gw_gamma = ci_beta)
      mydimnames[[3]][c(3,4)] <- c("Garwood", "Garwood(gamma)")
    }
  } else {
    mydimnames[[3]] <- c("SCAS_cc", "Jeffreys_cc", "mid-p_cc", "mid-p(beta)_cc",
                         "Wilson_cc", "Wald_cc", "Agresti-Coull_cc")
    if (distrib == "bin") {
      outlist <- list(scas_cc = ci_scas, jeff_cc = ci_jeff, beta_cc = ci_beta)
    } else if (distrib == "poi") {
      outlist <- list(scas_cc = ci_scas, jeff_cc = ci_jeff, gamma_cc = ci_beta)
      mydimnames[[3]][4] <- "mid-p(gamma)_cc"
    }
    # exact method not applicable if using a compromise value of cc
    # - but experimental beta/gamma distribution version included
  }
  if (distrib == "poi") mydimnames[[3]] <- mydimnames[[3]][-7]
  dimnames(outarr) <- mydimnames
  outarr <- aperm(outarr, c(3,2,1))

  call <- c(
    distrib = distrib, level = level, cc = cc, std_est = std_est
  )
  outlist <- append(outlist, list(ciarray = outarr, call = call))
  return(outlist)
}

#' Clopper-Pearson/Garwood and mid-p intervals for single binomial or
#' Poisson rate
#'
#' to calculate exact 'exact' confidence interval for a single binomial
#' or Poisson rate x/n
#'
#' @author Pete Laud, \email{p.j.laud@@sheffield.ac.uk}
#'
#' @noRd
exactci <- function(x,
                    n,
                    level = 0.95,
                    midp = TRUE,
                    beta = FALSE,
                    distrib = "bin",
                    precis = 8
                    ) {
  alpha <- 1 - level
  if (length(n) < length(x)) n <- rep(n, length(x))
  if (as.character(midp) == "TRUE") midp <- 0.5
  if (beta == FALSE) {
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
    est[x == 0 & n == 0] <- NaN
    upper <- bisect(
      ftn = uproot, precis = precis, uplow = "up", contrast = "p",
      distrib = distrib
    )
  } else if (beta == TRUE) {
    if (distrib == "bin") {
      est <- 0.5 * (qbeta(0.5, x, n - x + 1) + qbeta(0.5, x + 1, n - x))
      lower <- (1 - midp) * qbeta(alpha / 2, x, n - x + 1) + midp *
        qbeta(alpha / 2, x + 1, n - x)
      lower[x == n] <- qbeta(alpha / (2 * (1 - midp)), x, n - x + 1)[x == n]
      lower[x == 0] <- 0
      est[x == 0] <- 0
      upper <- (1 - midp) * qbeta(1 - alpha / 2, x + 1, n - x) + midp *
        qbeta(1 - alpha / 2, x, n - x + 1)
      upper[x == 0] <- qbeta(1 - alpha / (2 * (1 - midp)), x + 1, n - x)[x == 0]
      upper[x == n] <- 1
      est[x == n] <- 1
    } else if (distrib == "poi") {
      est <- n * (0.5 * qgamma(0.5, x, scale = 1 / n) + 0.5 *
                    qgamma(0.5, x + 1, scale = 1 / n))
      lower <- n * ((1 - midp) * qgamma(alpha / 2, x, scale = 1 / n) + (midp) *
                      qgamma(alpha / 2, x + 1, scale = 1 / n))
      lower[x == 0] <- 0
      est[x == 0] <- 0
      upper <- n * ((midp) * qgamma(1 - alpha / 2, x, scale = 1 / n) + (1 - midp) *
                      qgamma(1 - alpha / 2, x + 1, scale = 1 / n))
      upper[x == 0] <- (n * qgamma(1 - alpha / (2 * (1 - midp)),
                                   x + 1, scale = 1 / n))[x == 0]
    }
  }
  outdata <- cbind(lower = lower, est = est, upper = upper) /
    ifelse(distrib == "poi", n, 1)
  if (distrib == 'poi') {
    outdata[x == 0, 'lower'] <- 0
  }
  return(cbind(outdata, x = x, n = n))
}

