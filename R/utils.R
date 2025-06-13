#' Bisection root-finding
#'
#' vectorized limit-finding routine - turns out not to be any quicker but is
#' neater. The bisection method is just as efficient as the secant method
#' suggested by G&N, and affords greater control over whether the final estimate
#' has score<z the secant method is better for RR and for Poisson rates, where
#' there is no upper bound for d, however it is not guaranteed to converge New
#' version not reliant on point estimate This could be modified to solve upper
#' and lower limits simultaneously
#'
#' @inheritParams scoreci
#'
#' @author Pete Laud, \email{p.j.laud@@sheffield.ac.uk}
#'
#' @noRd
bisect <- function(ftn,
                   contrast,
                   distrib,
                   precis,
                   max.iter = 100,
                   uplow = "low") {
  tiny <- (10^-(precis)) / 2
  nstrat <- length(eval(ftn(1)))
  hi <- rep(1, nstrat)
  lo <- rep(-1, nstrat)
  dp <- 2
  dpt <- 1E12
  lastmidt <- 1E12
  niter <- 1
  while (niter <= max.iter && any(dpt > tiny | is.na(hi))) {
#  while (niter <= max.iter && any(dp > tiny | is.na(hi))) {
    dp <- 0.5 * dp
    mid <- pmax(-1, pmin(1, round((hi + lo) / 2, 15)))
    # rounding avoids machine precision problem with, e.g. 7/10-6/10
    if (contrast == "RD" && distrib == "bin") {
      midt <- mid
    } else if (contrast == "RD" && distrib == "poi") {
      midt <- round(tan(pi * mid / 2), 15)
    } else if (contrast %in% c("RR", "OR") ||
               (contrast == "p" && distrib == "poi")) {
      midt <- round(tan(pi * (mid + 1) / 4), 15)
    } else if (contrast == "p" && distrib == "bin") {
      midt <- (mid + 1) / 2
    }
    scor <- ftn(midt)
    dpt <- abs(midt - lastmidt)
    check <- (scor <= 0) | is.na(scor)
    # ??scor=NA only happens when |p1-p2|=1 and |theta|=1 for RD
    # (in which case hi==lo anyway), or if p1=p2=0
    # also for RR when p1=0 and theta=0, or paired RR when theta=0 or Inf
    hi[check] <- mid[check]
    lo[!check] <- mid[!check]
    niter <- niter + 1
    lastmidt <- midt
  }
  if (uplow == "low") {
    best <- lo
  } else {
    best <- hi
  }
  if (contrast == "RD" && distrib == "bin") {
    return(best)
  } else if ((contrast %in% c("RD") && distrib == "poi")) {
    return(tan(best * pi / 2))
  } else if (contrast %in% c("RR", "OR") ||
             (contrast == "p" && distrib == "poi")) {
    return(tan((best + 1) * pi / 4))
  } else if (contrast == "p" && distrib == "bin") {
    return((best + 1) / 2)
  }
}

#' Internal function - solve a quadratic
#'
#' @author Pete Laud, \email{p.j.laud@@sheffield.ac.uk}
#'
#' @noRd
quadroot <- function(a, b, c_) {
  # GET ROOTS OF A QUADRATIC EQUATION
  r1x <- (-b + sqrt(b^2 - 4 * a * c_)) / (2 * a)
  r2x <- (-b - sqrt(b^2 - 4 * a * c_)) / (2 * a)
  r1 <- pmin(r1x, r2x)
  r2 <- pmax(r1x, r2x)
  cbind(r1, r2)
}

#' Wilson score interval, and equivalent Rao score for Poisson data
#'
#' with optional continuity adjustment
#'
#' @param xihat Number specifying estimated variance inflation factor for the
#'   Saha Wilson Score interval for clustered binomial proportions
#' @author Pete Laud, \email{p.j.laud@@sheffield.ac.uk}
#' @references
#' Altman DG, Machin D, Bryant TN et al (2000) Statistics with confidence,
#' 2nd edn. BMJ Books, Bristol
#'
#' Li et al 2014. Comput Stat (2014) 29:869–889
#'
#' Labelled as "Second Normal" in REVSTAT – Statistical Journal
#' Volume 10, Number 2, June 2012, 211–227
#' (which provides the continuity adjustment formula)
#'
#' Schwertman, N.C. and Martinez, R.A. (1994). Approximate Poisson confidence
#' limits, Communication in Statistics — Theory and Methods, 23(5), 1507-1529.
#'
#' @noRd
wilsonci <- function(x,
                     n,
                     level = 0.95,
                     cc = FALSE,
                     xihat = 1,
                     distrib = "bin") {
  if (as.character(cc) == "TRUE") cc <- 0.5
  z <- qnorm(1 - (1 - level) / 2)
  est <- x / n
  if (distrib == "bin") {
    # Apply optional variance inflation factor for clustered data
    za <- z * sqrt(xihat)
    # See Newcombe 1998, set LCL to 0 if x = 0
    lower <- ifelse(x == 0,
                    0,
                    (2 * (x - cc) + za^2 -
                       za * sqrt(za^2 - 2 * (2 * cc + cc / n) +
                                   4 * ((x / n) * (n * (1 - x / n) + 2 * cc)))
                     ) / (2 * (n + za^2))
                    )
    upper <- ifelse(x == n,
                    1,
                    (2 * (x + cc) + za^2 +
                       za * sqrt(za^2 + 2 * (2 * cc - cc / n) +
                                   4 * ((x / n) * (n * (1 - x / n) - 2 * cc)))
                     ) / (2 * (n + za^2))
                    )
  } else if (distrib == "poi") {
    lower <- ((x - cc) + z^2 / 2 - z * sqrt(x - cc + z^2 / 4)) / n
    lower[x == 0] <- 0
    upper <- ((x + cc) + z^2 / 2 + z * sqrt(x + cc + z^2 / 4)) / n
  }
  cbind(lower = lower, est = est, upper = upper)
}

#' Rounding with trailing zeros
#'
#' @author Pete Laud, \email{p.j.laud@@sheffield.ac.uk}
#'
#' @noRd
myround <- function(x,
                    dp) {
  format(round(x, dp), nsmall = dp)
}

#' Egon Pearson Chi-squared test
#' https://rpubs.com/seriousstats/epcs_test
#'
#' For reference against tests produced by scoreci(..., skew = FALSE),
#' i.e. MN method
#'
#' @noRd
epcs.test <- function(data.cells,
                      z.adjust.method = c("none", "hommel")) {
  ucs.test <- suppressWarnings(stats::chisq.test(data.cells,
    simulate.p.value = FALSE,
    correct = FALSE
  ))
  N <- sum(data.cells)
  nrows <- dim(data.cells)[1]
  ncols <- dim(data.cells)[2]
  corrected.stat <- ucs.test$stat[[1]] * (N - 1) / N
  pval <- pchisq(corrected.stat, ucs.test$par, lower.tail = FALSE)
  output.list <- list(
    "Uncorrected Pearson Chi-Square" = ucs.test$stat,
    "Egon Pearson Chi-Square" = corrected.stat,
    "df" = prod(dim(data.cells) - 1), "p-value" = pval,
    "Smallest expected value (should be greater than 1)" = min(ucs.test$expected),
    "Raw data" = data.cells
  )
  return(output.list)
}


