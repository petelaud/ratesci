#' Confidence intervals for conditional odds ratio (OR) with paired binomial rates.
#'
#' @description
#' Confidence intervals for comparisons of two binomial rates from paired data.
#' This convenience wrapper function produces a selection of the methods below
#' for the conditional odds ratio (OR) contrast, with or without optional continuity
#' adjustment (where available).
#'
#' - Transformed SCAS (skewness-corrected asymptotic score)
#' - Transformed Wilson Score method
#' - Transformed mid-P
#' - Transformed Jeffreys
#' - Approximate log-normal (Wald) method
#'
#' @param x A numeric vector object specified as c(a, b, c, d)
#'   where: \cr
#'   a is the number of pairs with the event (e.g. success) under both
#'     conditions (e.g. treated/untreated, or case/control) \cr
#'   b is the count of the number with the event on condition 1 only (= x12) \cr
#'   c is the count of the number with the event on condition 2 only (= x21) \cr
#'   d is the number of pairs with no event under both conditions \cr
#'   (Note the order of a and d is only important for contrast="RR".)
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
#'   places) to be used in output.
#' @return A list containing the following components: \describe{
#'   \item{data}{the input data in 2x2 matrix form.}
#'   \item{estimates}{an array containing the confidence interval for paired OR
#'   using various methods. The methods shown depends on the cc argument
#'   (if cc = TRUE then the continuity-adjusted methods are given).}
#'   \item{call}{details of the function call.}
#'   }
#' @examples
#' # Example data from Fagerland et al 2014
#' orpairci(x = c(1, 1, 7, 12), precis = 3)
#' # with conventional continuity adjustment
#' orpairci(x = c(1, 1, 7, 12), precis = 3, cc = TRUE)
#' # with intermediate continuity adjustment
#' orpairci(x = c(1, 1, 7, 12), precis = 3, cc = 0.25)
#'
#' @author Pete Laud, \email{pete@@sheffstat.co.uk}
#' @references
#'   Fagerland MW, Lydersen S, Laake P. Recommended tests and
#'   confidence intervals for paired binomial proportions.
#'   Statistics in Medicine 2014; 33(16):2850-2875
#'
#'   Laud PJ. Improved confidence intervals and tests for paired binomial
#'   proportions. (2026, Under review)
#'
#' @export
orpairci <- function(x,
                     level = 0.95,
                     std_est = TRUE,
                     cc = FALSE,
                     precis = 8) {
  # Convert input data into 2x2 table to ease interpretation of output
  x1i <- rep(c("Success", "Success", "Failure", "Failure"), x)
  x2i <- rep(c("Success", "Failure", "Success", "Failure"), x)
  xi <- table(
    Test_1 = factor(x1i, levels = c("Success", "Failure")),
    Test_2 = factor(x2i, levels = c("Success", "Failure"))
  )

  contrast <- "OR"
  N <- sum(x)
  est <- x[2] / x[3]
  if (as.character(cc) == "TRUE") cc <- 0.5
  alpha <- 1 - level
  z <- qnorm(1 - alpha / 2)


  # special case for OR, use conditional method based on transforming the
  # SCAS interval for a single proportion
  # Transformed SCAS method with bias correction factor based on
  # N/(N-1) i.e. full sample size, not the number of discordant pairs
  # This gives a p-value matching that for other bias-corrected methods
  x12 <- x[2]
  x21 <- x[3]
  trans_ci <- scaspci(
    x = x12, n = x12 + x21, distrib = "bin",
    level = level, cc = cc, bcf = TRUE, bign = N
  )$estimates[, c(1:3), drop = FALSE]
  ci_scasp <- (trans_ci / (1 - trans_ci))
  trans_ci <- exactci(x = x12, n = x12 + x21, midp = 0.5 - cc, level = level)[, c(1:3)]
  ci_midp <- (trans_ci / (1 - trans_ci))
  trans_ci <- wilsonci(x = x12, n = x12 + x21, cc = cc, level = level)
  ci_wilson <- (trans_ci / (1 - trans_ci))
  trans_ci <- jeffreysci(x = x12, n = x12 + x21, cc = cc, level = level)$estimates[, c(1:3), drop = FALSE]
  ci_jeff <- (trans_ci / (1 - trans_ci))
  ci_wald <- exp(log(est) + c(-1, 0, 1) * z * sqrt(1 / x12 + 1 / x21))

  mydimnames <- dimnames(ci_scasp)

  methodnames <- c(
    "Transformed SCASp", "Transformed midp", "Transformed Wilson", "Transformed Jeffreys",
    "Wald"
  )

  mydimnames[[3]] <- methodnames

  outarr <- array(
    c(
      ci_scasp,
      ci_midp,
      ci_wilson,
      ci_jeff,
      ci_wald
    ),
    dim <- c(dim(ci_scasp), 5)
  )[drop = FALSE]
  dimnames(outarr) <- mydimnames

  if (std_est) outarr[, 2, ] <- est
  if (cc != FALSE) {
    methodnames <- paste0(methodnames, "_cc")
    if (cc != 0.5) methodnames <- paste0(methodnames, "(", cc, ")")
    mydimnames[[3]] <- methodnames
    dimnames(outarr) <- mydimnames
    if (cc == 0.5) methodnames[2] <- "Transformed Clopper-Pearson"
    mydimnames[[3]] <- methodnames
    dimnames(outarr) <- mydimnames
    outarr <- outarr[, , c(1:4), drop = FALSE]
  }
  outarr <- aperm(round(outarr, precis), c(3, 2, 1))[, , 1]

  call <- c(
    level = level,
    cc = cc
  )

  outlist <- list(xi, estimates = outarr, call = call)
  return(outlist)
}
