#' Confidence intervals for rate ratio (RR) with paired binomial rates.
#'
#' @description
#' Confidence intervals for comparisons of two binomial rates from paired data.
#' This convenience wrapper function produces a selection of the methods below
#' for the rate ratio (RR) contrast, with or without optional continuity
#' adjustment (where available).
#'
#' - SCAS (skewness-corrected asymptotic score)
#' - SCASu (omitting the 'N-1' adjustment)
#' - Tang Asymptotic Score method
#' - MOVER-W (based on Wilson method without Newcombe correlation adjustment)
#' - MOVER-NW (based on Wilson method with Newcombe correlation adjustment)
#' - MOVER-NJ (based on Jeffreys method with correlation adjustment)
#' - Bonett-Price hybrid method
#' - Bonett-Price-J variant using Jeffreys intervals
#' - Approximate log-normal (Wald) method
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
#'
#' @author Pete Laud, \email{p.j.laud@@sheffield.ac.uk}
#' @references
#'   Laud PJ. Improved confidence intervals and tests for paired binomial
#'   proportions. (2026, Under review)
#'
#' @export
rrpairci <- function(x,
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

  contrast <- "RR"
  N <- sum(x)
  est <- (x[1] + x[2]) / (x[1] + x[3])
  if (as.character(cc) == "TRUE") cc <- 0.5

  ci_wald <- rep(NA, 3)
  if (cc == FALSE) {
    ci_wald <- waldpairci(
      x = x,
      contrast = contrast,
      level = level,
      cc = cc
    )
  }
  ci_bp <- bpci(
    x = x,
    contrast = contrast,
    level = level,
    cc = cc
  )
  ci_bpj <- bpci(
    x = x,
    contrast = contrast,
    method = "jeff",
    level = level,
    cc = cc
  )

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

  ci_tang <- scorepairci(
    x = x,
    contrast = contrast,
    skew = FALSE,
    bcf = FALSE,
    level = level,
    cc = cc,
    precis = precis
  )$estimates[, c(1:3), drop = FALSE]

  ci_moveruw <- moverpairci(
    x = x,
    contrast = contrast,
    level = level,
    corc = FALSE,
    type = "wilson",
    cc = cc
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
    "SCAS", "SCASu", "Tang score", "MOVER-W", "MOVER-NW", "MOVER-NJ",
    "Wald", "Bonett-Price",  "Bonett-Price-J"
  )

  mydimnames[[3]] <- methodnames

  outarr <- array(
    c(
      ci_scas,
      ci_scasu,
      ci_tang,
      ci_moveruw,
      ci_moverw,
      ci_moverj,
      ci_wald,
      ci_bp,
      ci_bpj
    ),
    dim <- c(dim(ci_scas), 9)
  )[drop = FALSE]
  dimnames(outarr) <- mydimnames

  if (std_est) outarr[, 2, ] <- est
  if (cc != FALSE) {
    methodnames <- paste0(methodnames, "_cc")
    if (cc != 0.5) methodnames <- paste0(methodnames, "(", cc, ")")
    mydimnames[[3]] <- methodnames
    dimnames(outarr) <- mydimnames
    outarr <- outarr[, , c(1:6, 8), drop = FALSE]
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
