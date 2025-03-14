#' Confidence intervals for a single binomial rate from clustered data.
#'
#' Asymptotic Score confidence intervals for a proportion estimated from
#' a clustered sample, as decribed by Saha et al. 2016.
#' With optional skewness correction to improve interval location
#' (to be evaluated).
#'
#' @param x Numeric vector of number of events per cluster.
#' @param n Numeric vector of sample sizes per cluster.
#' @param level Number specifying confidence level (between 0 and 1, default
#'   0.95).
#' @param skew Logical (default TRUE) indicating whether to apply skewness
#'   correction or not. (To be evaluated)
#' @return A list containing the following components: \describe{
#'   \item{estimates}{the estimate and confidence interval for p and
#'   the specified confidence level, along with estimates of the ICC and
#'   the variance inflation factor, xihat}}
#' @examples
#'   # Data example from Liang 1992, used in Saha 2016:
#'   x <- c(rep(c(0, 1), c(36, 12)),
#'          rep(c(0, 1, 2), c(15, 7, 1)),
#'          rep(c(0, 1, 2, 3), c(5, 7, 3, 2)),
#'          rep(c(0, 1, 2), c(3, 3, 1)),
#'          c(0, 2, 3, 4, 6))
#'   n <- c(rep(1, 48),
#'          rep(2, 23),
#'          rep(3, 17),
#'          rep(4, 7),
#'          rep(6, 5))
#'   clusterpci(x, n, skew = FALSE)
#'   clusterpci(x, n, skew = TRUE)
#' @author Pete Laud, \email{p.j.laud@@sheffield.ac.uk}
#' @references
#'   Saha K, Miller D and Wang S. A comparison of some approximate confidence
#'   intervals for a single proportion for clustered binary outcome data.
#'   Int J Biostat 2016; 12: 1–18
#'
#'   Short MI et al. (2020) A novel confidence interval for a single proportion
#'   in the presence of clustered binary outcome data.
#'   Stat Meth Med Res 29(1) 111–121
#' @export
clusterpci <- function(x,
                       n,
                       level = 0.95,
                       skew = TRUE) {
  totx <- sum(x)
  # Notation as per Short et al.
  M1 <- sum(n)
  M2 <- sum(n^2)
  nn <- length(n)
  BMS <- (sum(x^2/n) - sum(x)^2 / M1) / (nn - 1)
  WMS <- (sum(x) - sum(x^2/n)) / sum(n-1)
  nstar <- (M1^2 - sum(n^2)) / ((nn - 1) * M1)
  phihat <- (BMS - WMS) / (BMS + (nstar - 1) * WMS) # ICC
  xihat <- 1 + phihat * (M2 - M1) / M1 # Variance inflation factor

  if (skew == TRUE) {
    out <- scaspci(totx, M1, xihat = xihat, level = level, bcf = FALSE)
  } else {
    out <- wilsonci(totx, M1, xihat = xihat, level = level)
  }
  newout <- cbind(out, totx = totx, totn = M1,  xihat = xihat, ICC = phihat)
  return(newout)
}

