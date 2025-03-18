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
#' with optional continuity correction
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
#' (which provides the continuity correction formula)
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
  cbind(Lower = lower, MLE = est, Upper = upper)
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

