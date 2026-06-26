#' Confidence intervals for rate ratio (RR) with independent binomial or Poisson rates.
#'
#' @description
#' Confidence intervals for comparisons of two binomial or Poisson rates.
#' This convenience wrapper function produces a selection of the methods below
#' as appropriate for the selected distribution (binomial or Poisson) for the
#' rate ratio (/relative risk) contrast (RR), with or without continuity adjustment.
#'
#' - SCAS (skewness-corrected asymptotic score)
#' - Miettinen-Nurminen, Koopman, Gart-Nam Asymptotic Score methods
#' - MOVER-W
#' - MOVER-J (based on Jeffreys intervals)
#' - Approximate normal (Katz log) methods
#'     (strongly advise this is not used for any purpose but included for reference)
#'
#' @param x1,x2 Numeric vectors of numbers of events in group 1 & group 2
#'   respectively.
#' @param n1,n2 Numeric vectors of sample sizes (for binomial rates) or exposure
#'   times (for Poisson rates) in each group.
#' @param distrib Character string indicating distribution assumed for the input
#'   data: \cr
#'   "bin" = binomial (default), \cr
#'   "poi" = Poisson.
#' @param level Number specifying confidence level (between 0 and 1, default
#'   0.95).
#' @param cc Number or logical (default FALSE) specifying (amount of) continuity
#'   adjustment. Numeric value between 0 and 0.5 is taken as the gamma parameter
#'   in Laud 2017, Appendix S2 (`cc = TRUE` translates to 0.5 for 'conventional'
#'   Yates adjustment). \cr
#' @param std_est logical, specifying if the crude point estimate for the
#'   contrast value should be returned (TRUE, default) or the
#'   method-specific alternative point estimate consistent with a 0% confidence
#'   interval (FALSE).
#' @param precis Number (default 6) specifying precision (i.e. number of decimal
#'   places) to be used in root-finding subroutine for the score confidence
#'   intervals. (Note other methods use closed-form calculations so are not
#'   affected.)
#' @return A list containing the following components: \describe{
#'   \item{estimates}{an array containing the confidence interval for RR
#'   using various methods. The methods shown depends on the cc argument
#'   (if cc = TRUE then the continuity-adjusted methods are given).}
#'   \item{call}{details of the function call.}
#'   }
#' @examples
#' # Selected example datasets from Newcombe 1998 and Fagerland et al. 2011
#' # (note Fagerland et al. appear to have the Miettinen-Nurminen method
#' #  labelled as Koopman)
#' rrci(
#'   x1 = c(5, 7), n1 = c(56, 34),
#'   x2 = c(0, 1), n2 = c(29, 34),
#'   precis = 4
#' )
#' # With conventional continuity adjustment
#' rrci(
#'   x1 = c(5, 7), n1 = c(56, 34),
#'   x2 = c(0, 1), n2 = c(29, 34),
#'   precis = 4, cc = TRUE
#' )
#' # With intermediate continuity adjustment
#' rrci(
#'   x1 = c(5, 7), n1 = c(56, 34),
#'   x2 = c(0, 1), n2 = c(29, 34),
#'   precis = 4, cc = 0.25
#' )
#'
#' @author Pete Laud, \email{pete@@sheffstat.co.uk}
#' @references
#'   Laud PJ. Equal-tailed confidence intervals for comparison of
#'   rates. Pharmaceutical Statistics 2017; 16:334-348.
#'
#' @export
rrci <- function(x1,
                 n1,
                 x2,
                 n2,
                 distrib = "bin",
                 level = 0.95,
                 std_est = TRUE,
                 cc = FALSE,
                 precis = 6) {
  if (length(x1) != length(x2)) {
    print("x1 and x2 must be the same length")
    stop()
  }
  nstrat <- length(x1)
  # in case x1,x2 are input as vectors but n1,n2 are not
  if (length(n1) < nstrat && nstrat > 1) n1 <- rep(n1, length.out = nstrat)
  if (length(n2) < nstrat && nstrat > 1) {
    n2 <- rep(n2, length.out = nstrat)
  }

  contrast <- "RR"
  est <- (x1 / n1) / (x2 / n2)
  if (as.character(cc) == "TRUE") cc <- 0.5

  ci_wald <- waldci(
    x1 = x1,
    n1 = n1,
    x2 = x2,
    n2 = n2,
    distrib = distrib,
    contrast = contrast,
    level = level,
    cc = cc
  )

  ci_adjwald <- rep(NA, 3)
  if (cc == FALSE) {
    ci_adjwald <- waldci(
      x1 = x1 + 0.5,
      n1 = n1 + 0.5,
      x2 = x2 + 0.5,
      n2 = n2 + 0.5,
      distrib = distrib,
      contrast = contrast,
      level = level,
      cc = cc
    )
  }

  ci_scas <- scasci(
    x1 = x1,
    n1 = n1,
    x2 = x2,
    n2 = n2,
    distrib = distrib,
    contrast = contrast,
    level = level,
    cc = cc,
    precis = precis
  )$estimates[, c(1:3), drop = FALSE]

  ci_gn <- scoreci(
    x1 = x1,
    n1 = n1,
    x2 = x2,
    n2 = n2,
    distrib = distrib,
    contrast = contrast,
    level = level,
    skew = TRUE,
    bcf = FALSE,
    or_bias = TRUE,
    cc = cc,
    precis = precis
  )$estimates[, c(1:3), drop = FALSE]

  ci_mn <- scoreci(
    x1 = x1,
    n1 = n1,
    x2 = x2,
    n2 = n2,
    distrib = distrib,
    contrast = contrast,
    level = level,
    skew = FALSE,
    bcf = TRUE,
    or_bias = FALSE,
    cc = cc,
    precis = precis
  )$estimates[, c(1:3), drop = FALSE]

  ci_mee <- scoreci(
    x1 = x1,
    n1 = n1,
    x2 = x2,
    n2 = n2,
    distrib = distrib,
    contrast = contrast,
    level = level,
    skew = FALSE,
    bcf = FALSE,
    or_bias = FALSE,
    cc = cc,
    precis = precis
  )$estimates[, c(1:3), drop = FALSE]

  ci_moverw <- moverci(
    x1 = x1,
    n1 = n1,
    x2 = x2,
    n2 = n2,
    distrib = distrib,
    contrast = contrast,
    level = level,
    type = "wilson",
    cc = cc
  )$estimates[, c(1:3), drop = FALSE]

  ci_moverj <- moverci(
    x1 = x1,
    n1 = n1,
    x2 = x2,
    n2 = n2,
    distrib = distrib,
    contrast = contrast,
    level = level,
    type = "jeff",
    adj = TRUE,
    cc = cc
  )$estimates[, c(1:3), drop = FALSE]

  mydimnames <- dimnames(ci_scas)
  mydimnames[[1]] <- paste0(x1, "/", n1, " vs ", x2, "/", n2)

  methodnames <- c(
    "SCAS", "Gart-Nam", "Miettinen-Nurminen", "Koopman",
    "MOVER-W", "MOVER-J",
    "Katz log", "Adjusted log"
  )
  if (distrib == "poi") methodnames[7] <- "Approximate Lognormal"

  mydimnames[[3]] <- methodnames

  outarr <- array(
    c(ci_scas, ci_gn, ci_mn, ci_mee, ci_moverw, ci_moverj, ci_wald, ci_adjwald),
    dim <- c(dim(ci_scas), 8)
  )[drop = FALSE]
  dimnames(outarr) <- mydimnames

  if (std_est) outarr[, 2, ] <- est
  if (cc != FALSE) {
    methodnames <- paste0(methodnames, "_cc")
    if (cc != 0.5) methodnames <- paste0(methodnames, "(", cc, ")")
    mydimnames[[3]] <- methodnames
    dimnames(outarr) <- mydimnames
    outarr <- outarr[, , 1:6, drop = FALSE]
  }
  if (distrib == "poi" && cc == FALSE) {
    mydimnames[[3]] <- methodnames
    dimnames(outarr) <- mydimnames
    outarr <- outarr[, , 1:7, drop = FALSE]
  }
  if (distrib == "poi") outarr <- outarr[, , -c(2, 4), drop = FALSE]
  # dimnames(outarr) <- mydimnames
  outarr <- aperm(round(outarr, precis), c(3, 2, 1))

  call <- c(
    distrib = distrib,
    level = level,
    cc = cc
  )

  outlist <- list(estimates = outarr, call = call)
  return(outlist)
}
