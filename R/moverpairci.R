#' MOVER confidence intervals for comparisons of paired binomial rates.
#'
#' Confidence intervals for the rate (or risk) difference ("RD"), or rate ratio
#' ("RR"), for paired binomial data.
#' This function applies the Method of Variance Estimates Recovery (MOVER)
#' for RD and RR, incorporating Newcombe's correlation correction.
#' All methods have options for continuity adjustment, where the magnitude of
#' adjustment can be customised.
#'
#' @param x A numeric vector object specified as c(a, b, c, d)
#'   where: \cr
#'   a is the number of pairs with the event (e.g. success) under both
#'     conditions (e.g. treated/untreated, or case/control) \cr
#'   b is the count of the number with the event on condition 1 only (= x12) \cr
#'   c is the count of the number with the event on condition 2 only (= x21) \cr
#'   d is the number of pairs with no event under both conditions \cr
#'   (Note the order of a and d is only important for contrast="RR".)
#' @param contrast Character string indicating the contrast of interest: \cr
#'   "RD" = rate difference (default); \cr
#'   "RR" = rate ratio.
#' @param level Number specifying confidence level (between 0 and 1, default
#'   0.95).
#' @param type Character string indicating the method used for the intervals for
#'   the individual group rates. \cr
#'   "jeff" = Jeffreys equal-tailed intervals (default); \cr
#'   "SCASp" = skewness-corrected score, \cr
#'   "midp" = mid-p, \cr
#'   "wilson" = Wilson score (not recommended, known to be skewed).
#' @param corc Logical (default TRUE) indicating whether to apply adjustment
#'   to the correlation estimate from Newcombe.
#' @param cc Number or logical (default FALSE) specifying (amount of) continuity
#'   adjustment. When a score-based method is used, cc = 0.5 corresponds to the
#'   continuity-corrected McNemar test.
#' @param precis Number (default 6) specifying precision (i.e. number of decimal
#'   places) to be used in output.
#' @param warn Logical (default TRUE) giving the option to suppress warnings.
#' @param ... Other arguments.
#' @importFrom stats uniroot pbinom ppois dpois
#' @return A list containing the following components: \describe{
#'   \item{data}{the input data in 2x2 matrix form.}
#'   \item{estimates}{the requested contrast, with its confidence interval and
#'   the specified confidence level, along with estimates of the marginal
#'   probabilities and the correlation coefficient (uncorrected and
#'   corrected).}
#'   \item{call}{details of the function call.}}
#' @examples
#' # Example data from Fagerland et al 2014
#' # MOVER-NJ method
#' moverpairci(x = c(1, 1, 7, 12), contrast = "RD", corc = TRUE, type = "jeff")
#' # MOVER-NJ
#' moverpairci(x = c(1, 1, 7, 12), contrast = "RR", corc = TRUE, type = "jeff")
#'
#' @author Pete Laud, \email{pete@@sheffstat.co.uk}
#' @references
#'   Newcombe RG. Improved confidence intervals for the difference between
#'   binomial proportions based on paired data.
#'   Statistics in Medicine 1998; 17:2635-2650
#'
#'   Tang M-L, Li H-Q, Tang N-S. Confidence interval construction for proportion
#'   ratio in paired studies based on hybrid method.
#'   Statistical Methods in Medical Research 2010; 21(4):361-378
#'
#'   Tang N-S et al. Asymptotic confidence interval construction for proportion
#'   difference in medical studies with bilateral data.
#'   Statistical Methods in Medical Research. 2011; 20(3):233-259
#'
#'   Fagerland MW, Lydersen S, Laake P. Recommended tests and
#'   confidence intervals for paired binomial proportions.
#'   Statistics in Medicine 2014; 33(16):2850-2875
#'
#'   Laud PJ. Equal-tailed confidence intervals for comparison of rates.
#'   Pharmaceutical Statistics 2017; 16:334-348.
#'
#'   DelRocco N et al. New Confidence Intervals for Relative Risk of Two
#'   Correlated Proportions.
#'   Statistics in Biosciences 2023; 15:1–30
#'
#'   Chang P et al. Continuity corrected score confidence interval for the
#'   difference in proportions in paired data.
#'   Journal of Applied Statistics 2024; 51-1:139-152
#'
#'   Laud PJ. Comments on "New Confidence Intervals for Relative Risk of Two
#'   Correlated Proportions" (2023).
#'   Statistics in Biosciences 2025; https://doi.org/10.1007/s12561-025-09479-4
#'
#'   Laud PJ. Improved confidence intervals and tests for paired binomial
#'   proportions. (2026, Under review)
#'
#' @export
moverpairci <- function(x,
                        level = 0.95,
                        contrast = "RD",
                        type = "jeff",
                        corc = TRUE,
                        cc = FALSE,
                        precis = 6,
                        warn = TRUE,
                        ...) {
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
  if (!(tolower(substr(contrast, 1, 2)) %in% c("rd", "rr", "or"))) {
    print("Contrast must be one of 'RD', 'RR' or 'OR'")
    stop()
  }
  if (!(tolower(substr(type, 1, 4)) %in%
    c("scas", "wils", "midp", "jeff"))) {
    print("type must be one of 'SCASp', 'wilson', 'midp' or 'jeff'")
    stop()
  }

  if (as.character(cc) == "TRUE") cc <- 0.5
  # Default correction aligned with cc'd McNemar test
  # Note that this is 0.5 for consistency with other functions in the ratesci
  # package, but for contrast = "RD", the formulation of the test statistic
  # means this is doubled.

  doublezero <- (x[1] + x[2] + x[3] == 0)
  nodiscord <- (x[2] == 0 && x[3] == 0)

  # Convert input data into 2x2 table to ease interpretation
  x1i <- rep(c("Success", "Success", "Failure", "Failure"), x)
  x2i <- rep(c("Success", "Failure", "Success", "Failure"), x)
  xi <- table(
    Test_1 = factor(x1i, levels = c("Success", "Failure")),
    Test_2 = factor(x2i, levels = c("Success", "Failure"))
  )
  x1 <- x[1] + x[2]
  x2 <- x[1] + x[3]
  N <- sum(x)
  p1hat <- x1 / N
  p2hat <- x2 / N

  # correlation estimate for reporting
  phi_hat <- (x[1] * x[4] - x[2] * x[3]) / sqrt(x1 * (N - x1) * x2 * (N - x2))
  # if (is.na(phi_hat) | is.infinite(phi_hat)) {
  if (is.na(phi_hat)) {
    phi_hat <- 0
  }
  # Newcombe's adjusted correlation estimate
  phi_c <- phi_hat
  if (x[1] * x[4] - x[2] * x[3] > 0) {
    phi_c <- (max(x[1] * x[4] - x[2] * x[3] - N / 2, 0) /
      sqrt(x1 * (N - x1) * x2 * (N - x2)))
  }

  # parameter used by Fagerland et al, for reporting
  psi_hat <- x[1] * x[4] / (x[2] * x[3])
  if (is.na(psi_hat)) psi_hat <- 0

  if (FALSE) {
    # MOVER methods for RD and RR
    if (method == "MOVER") {
      estimates <- moverpairci(
        x = x, contrast = contrast, level = level,
        method = type, cc = cc, corc = FALSE
      )
      #      outlist <- list(data = xi, estimates = estimates)
    }
    if (method == "MOVER_newc") {
      estimates <- moverpairci(
        x = x, contrast = contrast, level = level,
        method = type, cc = cc, corc = TRUE
      )
      #      outlist <- list(data = xi, estimates = estimates)
    }
  }

  x1p <- x[1] + x[2]
  xp1 <- x[1] + x[3]
  z <- qnorm(1 - (1 - level) / 2)
  ## First calculate l1, u1, l2, u2 based on an interval for p1p and pp1
  if (type == "SCASp") {
    j1 <- rateci(x = x1p, n = N, distrib = "bin", level = level, cc = cc)$scas
    j2 <- rateci(x = xp1, n = N, distrib = "bin", level = level, cc = cc)$scas
  } else if (type == "wilson") {
    j1 <- wilsonci(x = x1p, n = N, distrib = "bin", level = level, cc = cc)
    j2 <- wilsonci(x = xp1, n = N, distrib = "bin", level = level, cc = cc)
  } else if (type == "jeff") {
    j1 <- rateci(x = x1p, n = N, distrib = "bin", level = level, cc = cc)$jeff
    j2 <- rateci(x = xp1, n = N, distrib = "bin", level = level, cc = cc)$jeff
  } else if (type == "midp") {
    j1 <- rateci(x = x1p, n = N, distrib = "bin", level = level, cc = cc)$midp
    j2 <- rateci(x = xp1, n = N, distrib = "bin", level = level, cc = cc)$midp
  }
  l1 <- j1[, 1]
  u1 <- j1[, 3]
  l2 <- j2[, 1]
  u2 <- j2[, 3]
  #  Early version using medians for p1 & p2 for type="jeff" instead of x/N
  #  Slightly improves coverage, but has inferior location
  #  p1phat <- j1[, 2]
  #  pp1hat <- j2[, 2]
  p1phat <- x1p / N
  pp1hat <- xp1 / N

  n2p <- x[3] + x[4]
  np2 <- x[2] + x[4]
  cor_hat <- (x[1] * x[4] - x[2] * x[3]) / sqrt(x1p * n2p * xp1 * np2)
  if (corc == TRUE) {
    if (x[1] * x[4] - x[2] * x[3] > 0) {
      cor_hat <- (max(x[1] * x[4] - x[2] * x[3] - N / 2, 0) /
        sqrt(x1p * n2p * xp1 * np2))
    }
  }
  if (is.na(cor_hat) | is.infinite(cor_hat)) {
    cor_hat <- 0
  }

  if (contrast == "RR") {
    estimate <- p1phat / pp1hat
    A <- (p1phat - l1) * (u2 - pp1hat) * cor_hat
    B <- (u1 - p1phat) * (pp1hat - l2) * cor_hat
    lower <- (A - p1phat * pp1hat +
      sqrt(pmax(0, (A - p1phat * pp1hat)^2 -
        l1 * (2 * p1phat - l1) * u2 * (2 * pp1hat - u2)))) /
      (u2 * (u2 - 2 * pp1hat))
    upper <- (B - p1phat * pp1hat -
      sqrt(pmax(0, (B - p1phat * pp1hat)^2 -
        u1 * (2 * p1phat - u1) * l2 * (2 * pp1hat - l2)))) /
      (l2 * (l2 - 2 * pp1hat))
    if (is.na(upper) | upper < 0) {
      upper <- Inf
    }
  } else if (contrast == "RD") {
    estimate <- p1phat - pp1hat
    lower <- p1phat - pp1hat -
      sqrt((p1phat - l1)^2 + (u2 - pp1hat)^2 -
        2 * (p1phat - l1) * (u2 - pp1hat) * cor_hat)
    upper <- p1phat - pp1hat +
      sqrt((u1 - p1phat)^2 + (pp1hat - l2)^2 -
        2 * (u1 - p1phat) * (pp1hat - l2) * cor_hat)
  }

  estimates <- cbind(
    lower = lower, est = estimate, upper = upper, level = level,
    p1hat = p1phat, p2hat = pp1hat, phi_hat = cor_hat
  )
  row.names(estimates) <- NULL

  #  outlist <- list(data = xi, estimates = round(estimates, precis))

  call <- c(
    contrast = contrast, type = type,
    level = level, cc = cc
  )
  #  outlist <- append(outlist, list(call = call))
  outlist <- list(data = xi, estimates = round(estimates, precis), call = call)

  return(outlist)
}


#' MOVER interval for paired RR or RD (binomial only)
#'
#' Method of Variance Estimates Recovery, applied to paired RD and RR.
#' With various options for the marginal rates, and with optional continuity
#' adjustment, and Newcombe's correction to the Pearson correlation estimate,
#' applied to both contrasts.
#'
#' @author Pete Laud, \email{pete@@sheffstat.co.uk}
#' @references
#'   Newcombe RG. Improved confidence intervals for the difference between
#'   binomial proportions based on paired data.
#'   Statistics in Medicine 1998; 17:2635-2650
#'
#'   Tang M-L, Li H-Q, Tang N-S. Confidence interval construction for proportion
#'   ratio in paired studies based on hybrid method.
#'   Statistical Methods in Medical Research 2010; 21(4):361-378
#'
#'   Tang N-S et al. Asymptotic confidence interval construction for proportion
#'   difference in medical studies with bilateral data.
#'   Statistical Methods in Medical Research. 2011; 20(3):233-259
#'
#' @inheritParams pairbinci
#' @param corc Logical (default TRUE) indicating whether to apply adjustment
#'   to the correlation estimate from Newcombe.
#'
#' @noRd
xmoverpairci <- function(x,
                         level = 0.95,
                         contrast = "RD",
                         method = "jeff",
                         corc = TRUE,
                         cc = FALSE) {
  if (!is.numeric(c(x))) {
    print("Non-numeric inputs!")
    stop()
  }

  x1p <- x[1] + x[2]
  xp1 <- x[1] + x[3]
  N <- sum(x)
  z <- qnorm(1 - (1 - level) / 2)
  ## First calculate l1, u1, l2, u2 based on an interval for p1p and pp1
  if (method == "SCASp") {
    j1 <- rateci(x = x1p, n = N, distrib = "bin", level = level, cc = cc)$scas
    j2 <- rateci(x = xp1, n = N, distrib = "bin", level = level, cc = cc)$scas
  } else if (method == "wilson") {
    j1 <- wilsonci(x = x1p, n = N, distrib = "bin", level = level, cc = cc)
    j2 <- wilsonci(x = xp1, n = N, distrib = "bin", level = level, cc = cc)
  } else if (method == "jeff") {
    j1 <- rateci(x = x1p, n = N, distrib = "bin", level = level, cc = cc)$jeff
    j2 <- rateci(x = xp1, n = N, distrib = "bin", level = level, cc = cc)$jeff
  } else if (method == "midp") {
    j1 <- rateci(x = x1p, n = N, distrib = "bin", level = level, cc = cc)$midp
    j2 <- rateci(x = xp1, n = N, distrib = "bin", level = level, cc = cc)$midp
  }
  l1 <- j1[, 1]
  u1 <- j1[, 3]
  l2 <- j2[, 1]
  u2 <- j2[, 3]
  #  Early version using medians for p1 & p2 for method="jeff" instead of x/N
  #  Slightly improves coverage, but has inferior location
  #  p1phat <- j1[, 2]
  #  pp1hat <- j2[, 2]
  p1phat <- x1p / N
  pp1hat <- xp1 / N

  n2p <- x[3] + x[4]
  np2 <- x[2] + x[4]
  cor_hat <- (x[1] * x[4] - x[2] * x[3]) / sqrt(x1p * n2p * xp1 * np2)
  if (corc == TRUE) {
    if (x[1] * x[4] - x[2] * x[3] > 0) {
      cor_hat <- (max(x[1] * x[4] - x[2] * x[3] - N / 2, 0) /
        sqrt(x1p * n2p * xp1 * np2))
    }
  }
  if (is.na(cor_hat) | is.infinite(cor_hat)) {
    cor_hat <- 0
  }

  if (contrast == "RR") {
    estimate <- p1phat / pp1hat
    A <- (p1phat - l1) * (u2 - pp1hat) * cor_hat
    B <- (u1 - p1phat) * (pp1hat - l2) * cor_hat
    lower <- (A - p1phat * pp1hat +
      sqrt(pmax(0, (A - p1phat * pp1hat)^2 -
        l1 * (2 * p1phat - l1) * u2 * (2 * pp1hat - u2)))) /
      (u2 * (u2 - 2 * pp1hat))
    upper <- (B - p1phat * pp1hat -
      sqrt(pmax(0, (B - p1phat * pp1hat)^2 -
        u1 * (2 * p1phat - u1) * l2 * (2 * pp1hat - l2)))) /
      (l2 * (l2 - 2 * pp1hat))
    if (is.na(upper) | upper < 0) {
      upper <- Inf
    }
  } else if (contrast == "RD") {
    estimate <- p1phat - pp1hat
    lower <- p1phat - pp1hat -
      sqrt((p1phat - l1)^2 + (u2 - pp1hat)^2 -
        2 * (p1phat - l1) * (u2 - pp1hat) * cor_hat)
    upper <- p1phat - pp1hat +
      sqrt((u1 - p1phat)^2 + (pp1hat - l2)^2 -
        2 * (u1 - p1phat) * (pp1hat - l2) * cor_hat)
  }

  estimates <- cbind(
    lower = lower, est = estimate, upper = upper, level = level,
    p1hat = p1phat, p2hat = pp1hat, phi_hat = cor_hat
  )
  row.names(estimates) <- NULL
  return(estimates)
}
