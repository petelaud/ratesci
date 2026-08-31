#' Confidence intervals for rate difference (RD) with paired binomial rates.
#'
#' @description
#' Confidence intervals for comparisons of two binomial rates from paired data.
#' This convenience wrapper function produces a selection of the methods below
#' for the rate difference (RD) contrast, with or without optional continuity
#' adjustment (where available).
#'
#' - SCAS (skewness-corrected asymptotic score)
#' - SCASu (omitting the 'N-1' adjustment)
#' - Tango Asymptotic Score method
#' - MOVER-NW (aka Newcombe Hybrid Score or square-and-add)
#' - MOVER-NJ (based on Jeffreys intervals)
#' - Agresti-Min
#' - Bonett-Price
#' - Approximate normal (Wald) method
#'     (strongly advise this is not used for any purpose but included for reference)
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
#' @param precis Number (default 8) specifying precision (i.e. number of decimal
#'   places) to be used in root-finding subroutine for the score confidence
#'   intervals. (Note other methods use closed-form calculations so are not
#'   affected.)
#' @return A list containing the following components: \describe{
#'   \item{data}{the input data in 2x2 matrix form.}
#'   \item{estimates}{an array containing the confidence interval for paired RD
#'   using various methods. The methods shown depends on the cc argument
#'   (if cc = TRUE then the continuity-adjusted methods are given).}
#'   \item{call}{details of the function call.}
#'   }
#' @examples
#' # Example data from Fagerland et al 2014
#' rdpairci(x = c(1, 1, 7, 12), precis = 3)
#' # with conventional continuity adjustment
#' rdpairci(x = c(1, 1, 7, 12), precis = 3, cc = TRUE)
#' # with intermediate continuity adjustment
#' rdpairci(x = c(1, 1, 7, 12), precis = 3, cc = 0.25)
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
rdpairci <- function(x,
                     level = 0.95,
                     std_est = TRUE,
                     cc = FALSE,
                     precis = 8) {

  contrast <- "RD"

  # Input checks
  if (!is.numeric(c(x))) {
    print("Non-numeric inputs!")
    stop()
  }
  if (!(is.vector(x) && length(x) == 4)) {
    print("Input x must be a vector of length 4!")
    stop()
  }
  if (any(x < 0)) {
    print("Negative inputs!")
    stop()
  }
  if (sum(x) == 0) {
    print("Sample size is zero!")
    stop()
  }

  # Convert input data into 2x2 table to ease interpretation of output
  x1i <- rep(c("Success", "Success", "Failure", "Failure"), x)
  x2i <- rep(c("Success", "Failure", "Success", "Failure"), x)
  xi <- table(
    Test_1 = factor(x1i, levels = c("Success", "Failure")),
    Test_2 = factor(x2i, levels = c("Success", "Failure"))
  )

  N <- sum(x)
  est <- (x[2] - x[3]) / N
  if (as.character(cc) == "TRUE") cc <- 0.5

  ci_agrmin <- ci_bp <- rep(NA, 3)
  ci_wald <- waldpairci(
    x = x,
    contrast = contrast,
    level = level,
    cc = cc
  )
  if (cc == FALSE) {
    ci_agrmin <- waldpairci(
      x = x + c(0.5, 0.5, 0.5, 0.5),
      contrast = contrast,
      level = level,
      cc = cc
    )
    ci_bp <- waldpairci(
      x = x + c(0, 1, 1, 0),
      contrast = contrast,
      level = level,
      cc = cc
    )
  }

  ci_scas <- scorepairci(
    x = x,
    contrast = contrast,
    level = level,
    cc = cc,
    precis = precis
  )$estimates[, c(1:3), drop = FALSE]

  ci_scasu <- scorepairci(
    x = x,
    contrast = contrast,
    bcf = FALSE,
    level = level,
    cc = cc,
    precis = precis
  )$estimates[, c(1:3), drop = FALSE]

  ci_tango <- scorepairci(
    x = x,
    contrast = contrast,
    skew = FALSE,
    bcf = FALSE,
    level = level,
    cc = cc,
    precis = precis
  )$estimates[, c(1:3), drop = FALSE]

  ci_moverw <- moverpairci(
    x = x,
    contrast = contrast,
    level = level,
    corc = TRUE,
    type = "wilson",
    cc = cc
  )$estimates[, c(1:3), drop = FALSE]

  ci_moverj <- moverpairci(
    x = x,
    contrast = contrast,
    level = level,
    corc = TRUE,
    type = "jeff",
    cc = cc
  )$estimates[, c(1:3), drop = FALSE]

  mydimnames <- dimnames(ci_scas)

  methodnames <- c(
    "SCAS", "SCASu", "Tango score", "MOVER-NW", "MOVER-NJ",
    "Wald", "Agresti-Min", "Bonett-Price"
  )

  mydimnames[[3]] <- methodnames

  outarr <- array(
    c(ci_scas, ci_scasu, ci_tango, ci_moverw, ci_moverj, ci_wald, ci_agrmin, ci_bp),
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
  # dimnames(outarr) <- mydimnames
  outarr <- aperm(round(outarr, precis), c(3, 2, 1))[, , 1]

  call <- c(
    level = level,
    cc = cc
  )

  outlist <- list(xi, estimates = outarr, call = call)
  return(outlist)
}


#' Wald approximate normal confidence intervals for a paired difference (RD)
#' or ratio (RR) or conditional odds ratio (OR)
#'
#' R code to calculate Wald CIs for contrast = "RD" or "RR", reluctantly
#' included for reference
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
#' @param contrast Character string indicating the contrast of interest: \cr
#'   "RD" = rate difference (default); \cr
#'   "RR" = rate ratio; \cr
#'   "OR" = conditional odds ratio.
#' @param cc Number or logical (default FALSE) specifying (amount of) continuity
#'   adjustment. When a score-based method is used, cc = 0.5 corresponds to the
#'   continuity-corrected McNemar test.
#' @param laplace Logical indicating whether to apply the Laplace adjustment
#'   for the OR contrast. Has no effect for other contrasts.
#' @author Pete Laud, \email{p.j.laud@@sheffield.ac.uk}
#' @references
#'   Fagerland MW, Lydersen S, Laake P. Recommended tests and
#'   confidence intervals for paired binomial proportions.
#'   Statistics in Medicine 2014; 33(16):2850-2875
#'
#' @noRd
waldpairci <- function(x,
                       level = 0.95,
                       contrast = "RD",
                       laplace = FALSE,
                       cc = FALSE) {
  x11 <- x[1]
  x10 <- x[2]
  x01 <- x[3]
  z0 <- qnorm(1 - (1 - level) / 2)
  x1 <- x11 + x10
  x0 <- x11 + x01
  if (as.character(cc) == "TRUE") cc <- 0.5

  if (contrast == "RD") {
    n <- sum(x)
    p12 <- (x10) / (n)
    p21 <- (x01) / (n)
    estimate <- p12 - p21
    v <- max(0, (p12 + p21 - (abs(p12 - p21) - 2 * cc/n)^2) / (n))
    estimates <- cbind(
      lower = pmax(-1, estimate - z0 * sqrt(v)),
      est = estimate,
      upper = pmin(1, estimate + z0 * sqrt(v))
    )
  } else if (contrast == "RR") {
    logestimate <- log(x1 / x0)
    v <- (x10 + x01) / (x1 * x0)
    estimates <- cbind(
      lower = ifelse(x0 == 0, 0, exp(logestimate - z0 * sqrt(v))),
      est = exp(logestimate),
      upper = ifelse(x1 == 0, Inf, exp(logestimate + z0 * sqrt(v)))
    )
    row.names(estimates) <- NULL
  } else if (contrast == "OR") {
    lapadd <- 1 * laplace
    xx10 <- x10 + lapadd
    xx01 <- x01 + lapadd
    logestimate <- log(xx10 / xx01)
    v <- 1 / xx10 + 1 / xx01
    estimates <- cbind(
      lower = ifelse((xx10 == 0 | xx01 == 0), 0, exp(logestimate - z0 * sqrt(v))),
      est = exp(logestimate),
      upper = ifelse((xx01 == 0 | xx10 == 0), Inf, exp(logestimate + z0 * sqrt(v)))
    )
  }
  list(estimates = estimates)
}

