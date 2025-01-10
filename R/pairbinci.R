#' Confidence intervals for comparisons of paired binomial rates.
#'
#' Confidence intervals for the rate (or risk) difference ('RD'),
#' rate ratio ('RR') or odds ratio ('OR'), for paired binomial data.
#' (For paired Poisson rates, suggest use the tdasci function with distrib = "poi",
#' and weighting = "MH", with pairs as strata.)
#' This function applies the score-based Tango and Tang methods for RD and
#' RR respectively, with iterative and closed-form versions, as well as an
#' experimental method using the stratified TDAS method with pairs as strata.
#' Also includes MOVER options using the Method of Variance Estimates Recovery
#' for paired RD and RR.
#' For OR, intervals are produced based on transforming various intervals for
#' the single proportion, including SCAS, mid-p and Jeffreys.
#' All methods have options for continuity correction, and the magnitude of
#' correction can be customised.
#'
#' @param x A numeric vector object specified as c(a, b, c, d)
#'   where:
#'   a is the number of pairs with the event (e.g. success) under both
#'     conditions (e.g. treated/untreated, or case/control)
#'   b is the count of the number with the event on condition 1 only (= n12)
#'   c is the count of the number with the event on condition 2 only (= n21)
#'   d is the number of pairs with no event under both conditions
#'   (Note the order of a and d is only important for contrast="RR".)
#' @param contrast Character string indicating the contrast of interest:
#'   "RD" = rate difference (default), "RR" = rate ratio, "OR" = odds ratio.
#' @param level Number specifying confidence level (between 0 and 1, default
#'   0.95).
#' @param cc Number or logical (default FALSE) specifying (amount of) continuity
#'   correction. When a score-based method is used, cc = 0.5 corresponds to the
#'   continuity-corrected McNemar test.
#' @param cctype Character string indicating the type of continuity correction
#'   ("constant" or "delrocco") to be applied for contrast = "RR".
#'   (Note: both options produce non-equivariant intervals, an improved
#'   correction is currently under evaluation - watch this space.)
#' @param theta0 Number to be used in a one-sided significance test (e.g.
#'   non-inferiority margin). 1-sided p-value will be < 0.025 iff 2-sided 95\% CI
#'   excludes theta0. NB: can also be used for a superiority test by setting
#'   theta0 = 0.
#' @param method_RD Character string indicating the confidence interval method
#'   to be used for contrast = "RD". "Score" = iterative Tango asymptotic score,
#'   "Score_closed" = closed form solution for Tango interval (default),
#'   "MOVER" = hybrid MOVER method for paired RD (as per "method 8" in
#'             Newcombe, but with a choice of input methods - see moverbase),
#'   "MOVER_newc" = hybrid MOVER method with "correction" to correlation
#'                  estimate (Newcombe's "method 10"),
#'   "TDAS" = t-distribution asymptotic score (experimental method, seems to
#'   struggle with low numbers).
#' @param method_RR Character string indicating the confidence interval method
#'   to be used for contrast = "RR". "Score" = iterative Tang asymptotic score,
#'   "Score_closed" = closed form solution for Tang interval (default),
#'   "MOVER" = hybrid MOVER method for paired RR,
#'   "MOVER_newc" = hybrid MOVER method with "correction" to correlation
#'                  estimate from Newcombe's RD method,
#'   "TDAS" = t-distribution asymptotic score (experimental method, seems to
#'   struggle with low numbers).
#' @param method_OR Character string indicating the confidence interval method
#'   to be used for contrast = "OR", all of which are based on transformation of
#'   an interval for a single proportion b/(b+c):
#'   "SCAS" = transformed skewness-corrected score (default),
#'   "jeff" = transformed Jeffreys,
#'   "midp" = transformed mid-p,
#'   "wilson" = transformed Wilson score - included for reference only, not
#'   recommended.
#' @param moverbase Character string indicating the base method used as input
#'   for the MOVER methods (when method_RR or method_RD = "MOVER"):
#'   "wilson" = Wilson score (not recommended, known to be skewed),
#'   "jeff" = Jeffreys equal-tailed interval,
#'   "midp" = mid-p,
#'   "SCAS" = skewness-corrected score (default)
#' @param precis Number (default 6) specifying precision (i.e. number of decimal
#'   places) to be used in optimisation subroutine for the confidence interval.
#' @importFrom stats uniroot pbinom ppois dpois
#' @importFrom polynom polynomial
#' @return A list containing the following components: \describe{
#'   \item{data}{the input data in 2x2 matrix form}
#'   \item{estimates}{the requested contrast, with its confidence interval and
#'   the specified confidence level}
#'   \item{pval}{the corresponding 2-sided significance test
#'   against the null hypothesis that p_1 = p_2 (McNemar's test), and one-sided
#'   significance tests against the null hypothesis that theta >= or <= theta0
#'   as specified}}
#' @examples
#' # Data example from Agresti-Min 2005
#' pairbinci(x = c(53, 16, 8, 9), contrast = "RD", method_RD = "Score")
#' pairbinci(x = c(53, 16, 8, 9), contrast = "RD", method_RD = "TDAS")
#' pairbinci(x = c(53, 16, 8, 9), contrast = "RR", method_RR = "Score")
#' pairbinci(x = c(53, 16, 8, 9), contrast = "RR", method_RR = "TDAS")
#' pairbinci(x = c(53, 16, 8, 9), contrast = "OR", method_OR = "SCAS")
#' # Example from Fagerland et al 2014
#' pairbinci(x = c(1, 1, 7, 12), contrast = "RD", method_RD = "Score_closed")
#' pairbinci(x = c(1, 1, 7, 12), contrast = "RD", method_RD = "MOVER_newc", moverbase = "wilson")
#' pairbinci(x = c(1, 1, 7, 12), contrast = "RD", method_RD = "MOVER_newc", moverbase = "SCAS")
#' pairbinci(x = c(1, 1, 7, 12), contrast = "RR", method_RR = "Score_closed")
#' pairbinci(x = c(1, 1, 7, 12), contrast = "RR", method_RR = "MOVER_newc", moverbase = "wilson")
#' pairbinci(x = c(1, 1, 7, 12), contrast = "RR", method_RR = "MOVER_newc", moverbase = "SCAS")
#' pairbinci(x = c(1, 1, 7, 12), contrast = "OR", method_OR = "wilson")
#' pairbinci(x = c(1, 1, 7, 12), contrast = "OR", method_OR = "midp")
#' pairbinci(x = c(1, 1, 7, 12), contrast = "OR", method_OR = "SCAS")
#'
#' @author Pete Laud, \email{p.j.laud@@sheffield.ac.uk}
#' @references
#'   Tango T. Equivalence test and confidence interval for the difference
#'   in proportions for the paired-sample design.
#'   Statistics in Medicine 1998; 17:891-908
#'
#'   Newcombe RG. Improved confidence intervals for the difference between
#'   binomial proportions based on paired data.
#'   Statistics in Medicine 1998; 17:2635-2650
#'
#'   Tango T. Improved confidence intervals for the difference between binomial
#'   proportions based on paired data by Robert G. Newcombe, Statistics in
#'   Medicine, 17, 2635-2650 (1998).
#'   Statistics in Medicine 1999; 18(24):3511-3513
#'
#'   Nam J-M, Blackwelder WC. Analysis of the ratio of marginal
#'   probabilities in a matched-pair setting.
#'   Stat Med 2002; 21(5):689–699
#'
#'   Tang N-S, Tang M-L, Chan ISF. On tests of equivalence via non-unity
#'   relative risk for matched-pair design.
#'   Statistics in Medicine 2003; 22:1217-1233
#'
#'   Agresti A, Min Y. Simple improved confidence intervals for
#'   comparing matched proportions.
#'   Statistics in Medicine 2005; 24:729-740
#'
#'   Tang M-L, Li H-Q, Tang N-S. Confidence interval construction for proportion
#'   ratio in paired studies based on hybrid method.
#'   Statistical Methods in Medical Research 2010; 21(4):361-378
#'
#'   Tang N-S et al. Asymptotic confidence interval construction for proportion
#'   difference in medical studies with bilateral data.
#'   Statistical Methods in Medical Research. 2011; 20(3):233-259
#'
#'   Yang Z, Sun X and Hardin JW. A non-iterative implementation of Tango's
#'   score confidence interval for a paired difference of proportions.
#'   Statistics in Medicine 2013; 32:1336-1342
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
#' @export
pairbinci <- function(x,
                      contrast = "RD",
                      level = 0.95,
                      cc = FALSE,
                      cctype = "constant",
                      method_RD = "Score_closed",
                      method_RR = "Score_closed",
                      method_OR = "SCAS",
                      moverbase = "SCAS",
                      theta0 = NULL,
                      precis = 6) {
  if (!(tolower(substr(contrast, 1, 2)) %in% c("rd", "rr", "or"))) {
    print("Contrast must be one of 'RD', 'RR' or 'OR'")
    stop()
  }
  if (!(tolower(substr(method_RD, 1, 4)) %in%
    c("tdas", "scor", "move"))) {
    print("Method_RD must be one of 'Score_closed', 'Score', 'TDAS', 'MOVER' or 'MOVER_newc'")
    stop()
  }
  if (!(tolower(substr(method_RR, 1, 4)) %in%
    c("tdas", "scor", "move"))) {
    print("Method_RR must be one of 'Score_closed', 'Score', 'TDAS', 'MOVER' or 'MOVER_newc'")
    stop()
  }
  if (!(tolower(substr(method_OR, 1, 4)) %in%
    c("scas", "wils", "midp", "jeff"))) {
    print("Method_OR must be one of 'SCAS', 'wilson', 'midp' or 'jeff'")
    stop()
  }
  if (!(tolower(substr(moverbase, 1, 4)) %in%
    c("scas", "wils", "midp", "jeff"))) {
    print("moverbase must be one of 'SCAS', 'wilson', 'midp' or 'jeff'")
    stop()
  }
  if (!is.numeric(c(x))) {
    print("Non-numeric inputs!")
    stop()
  }
  if (as.character(cc) == "TRUE") cc <- 0.5
  # Default correction aligned with cc'd McNemar test
  # Note that this is 0.5 for consistency with other functions in the ratesci
  # package, but for contrast = "RD", the formulation of the test statistic
  # means this is doubled.

  doublezero <- (x[1] + x[2] + x[3] == 0)
  nodiscord <- (x[2] == 0 && x[3] == 0)

  # Convert the data into 2 columns of 0s and 1s for use in TDAS-based method
  # and output 2x2 table for validation
  x1i <- rep(c(1, 1, 0, 0), x)
  x2i <- rep(c(1, 0, 1, 0), x)
  xi <- table(x1i, x2i)

  if (contrast == "OR") {
    # special case for OR, use conditional method based on transforming the
    # SCAS interval for a single proportion
    b <- x[2]
    c <- x[3]
    if (method_OR == "SCAS") {
      trans_th0 <- NULL
      if (is.null(theta0)) theta0 <- 1
      trans_th0 <- theta0 / (1 + theta0)
      OR_ci <- scaspci(
        x = b, n = b + c, distrib = "bin",
        level = level, cc = cc
      )
      estimates <- OR_ci / (1 - OR_ci)
      scorezero <- scoretheta(
        theta = 0.5, x1 = b, n1 = b + c,
        contrast = "p", distrib = "bin",
        skew = TRUE, cc = cc, level = level
      )
      chisq_zero <- scorezero$score^2
      pval2sided <- pchisq(chisq_zero, 1, lower.tail = FALSE)
      scoreth0 <- scoretheta(
        theta = trans_th0, x1 = b, n1 = b + c,
        contrast = "p", distrib = "bin",
        skew = TRUE, cc = cc, level = level
      )
      pval_left <- scoreth0$pval
      pval_right <- 1 - pval_left
      scorenull <- scoreth0$score
      pval <- cbind(
        chisq = chisq_zero, pval2sided, theta0 = theta0,
        scorenull, pval_left, pval_right
      )
      outlist <- list(xi, estimates = estimates, pval = pval)
    } else if (method_OR == "midp") {
      trans_ci <- exactci(x = b, n = b + c, midp = 0.5 - cc, level = level)
      estimates <- (trans_ci / (1 - trans_ci))
      outlist <- list(xi, estimates = estimates)
    } else if (method_OR == "wilson") {
      trans_ci <- wilsonci(x = b, n = b + c, cc = cc, level = level)
      estimates <- (trans_ci / (1 - trans_ci))
      outlist <- list(xi, estimates = estimates)
    } else if (method_OR == "jeff") {
      trans_ci <- jeffreysci(x = b, n = b + c, cc = cc, level = level)
      estimates <- (trans_ci / (1 - trans_ci))
      outlist <- list(xi, estimates = estimates)
    }
  } else if (contrast != "OR") {
    if ((contrast == "RD" && method_RD == "TDAS") ||
      (contrast == "RR" && method_RR == "TDAS")) {
      # stratified TDAS method for paired data as suggested in Laud 2017
      n1i <- n2i <- rep(1, sum(x))
      if (doublezero && contrast == "RR") {
        outlist <- list(
          data = xi,
          estimates = cbind(Lower = 0, MLE = NA, Upper = Inf),
          pval = NA
        )
      } else if (nodiscord && contrast == "RR") {
        out <- scoreci(
          stratified = TRUE, random = FALSE,
          x1 = x1i, n1 = n1i, x2 = x2i, n2 = n2i, weighting = "MH",
          contrast = contrast, distrib = "bin", level = level, cc = cc,
          theta0 = theta0, warn = FALSE
        )
        outlist <- list(data = xi, estimates = out$estimates[, 1:3], pval = out$pval)
      } else {
        out <- tdasci(
          x1 = x1i, n1 = n1i, x2 = x2i, n2 = n2i, weighting = "MH",
          contrast = contrast, distrib = "bin", level = level, cc = cc,
          theta0 = theta0, warn = FALSE
        )
        outlist <- list(
          data = xi,
          estimates = out$estimates[, 1:3, drop = FALSE],
          pval = out$pval
        )
      }
    }
    # Iterative Score methods by Tango (for RD) & Tang (for RR):
    if ((contrast == "RD" && method_RD == "Score") ||
      (contrast == "RR" && method_RR == "Score")) {
      myfun <- function(theta) {
        scorepair(
          theta = theta, x = x, contrast = contrast, cc = cc,
          cctype = cctype
        )$score
      }
      myfun0 <- function(theta) {
        scorepair(
          theta = theta, x = x, contrast = contrast, cc = 0,
          cctype = cctype
        )$score
      }
      # Use bisection routine to locate lower and upper confidence limits
      qtnorm <- qnorm(1 - (1 - level) / 2)
      MLE <- bisect(
        ftn = function(theta) myfun0(theta) - 0, distrib = "bin",
        contrast = contrast, precis = precis + 1, uplow = "low"
      )
      lower <- bisect(
        ftn = function(theta) myfun(theta) - qtnorm,
        distrib = "bin", precis = precis + 1, contrast = contrast,
        uplow = "low"
      )
      upper <- bisect(
        ftn = function(theta) myfun(theta) + qtnorm,
        distrib = "bin", precis = precis + 1, contrast = contrast,
        uplow = "up"
      )
      estimates <- cbind(
        Lower = lower, MLE = MLE, Upper = upper,
        level = level
      )
      # Closed-form versions of Score methods by Chang (for RD Tango method)
      #  & DelRocco (for RR Tang method):
    } else if (contrast == "RD" && method_RD == "Score_closed") {
      estimates <- tangoci(x = x, level = level, cc = cc)
    } else if (contrast == "RR" && method_RR == "Score_closed") {
      estimates <- tangci(x = x, level = level, cc = cc, cctype = cctype)
    }
    # MOVER methods for RD and RR
    if ((contrast == "RD" && method_RD == "MOVER")) {
      estimates <- moverpair(
        x = x, contrast = "RD", level = level,
        method = moverbase, cc = cc
      )
      outlist <- list(data = xi, estimates = estimates)
    }
    if ((contrast == "RD" && method_RD == "MOVER_newc")) {
      estimates <- moverpair(
        x = x, contrast = "RD", level = level,
        method = moverbase, cc = cc, corc = TRUE
      )
      outlist <- list(data = xi, estimates = estimates)
    }
    if ((contrast == "RR" && method_RR == "MOVER")) {
      estimates <- moverpair(
        x = x, contrast = "RR", level = level,
        method = moverbase, cc = cc
      )
      outlist <- list(data = xi, estimates = estimates)
    }
    if ((contrast == "RR" && method_RR == "MOVER_newc")) {
      estimates <- moverpair(
        x = x, contrast = "RR", level = level,
        method = moverbase, cc = cc, corc = TRUE
      )
      outlist <- list(data = xi, estimates = estimates)
    }

    if ((contrast == "RD" && method_RD == "Score") ||
      (contrast == "RR" && method_RR == "Score") ||
      (contrast == "RD" && method_RD == "Score_closed") ||
      (contrast == "RR" && method_RR == "Score_closed")) {
      # optionally add p-value for a test of null hypothesis: theta<=theta0
      # default value of theta0 depends on contrast
      if (contrast == "RD") {
        theta00 <- 0
      } else {
        theta00 <- 1
      }
      if (is.null(theta0)) {
        # If null value for theta is not provided, use 0 for RD or 1 for RR
        theta0 <- theta00
      }
      scorezero <- scorepair(
        theta = theta00, x = x, contrast = contrast,
        cc = cc, cctype = cctype
      )
      scorenull <- scorepair(
        theta = theta0, x = x, contrast = contrast,
        cc = cc, cctype = cctype
      )
      pval_left <- scorenull$pval
      pval_right <- 1 - pval_left
      chisq_zero <- scorezero$score^2
      pval2sided <- pchisq(chisq_zero, 1, lower.tail = FALSE)
      pval <- cbind(
        chisq = chisq_zero, pval2sided, theta0 = theta0,
        scorenull = scorenull$score, pval_left, pval_right
      )

      outlist <- list(data = xi, estimates = estimates, pval = pval)
    }
  }
  return(outlist)
}


#' Internal function to evaluate the score at a given value of theta
#'
#' For iterative calculations:
#' function to evaluate the score at a given value of theta, given the observed
#' data for paired binomial RD and RR. uses the MLE solution (and notation)
#' given in Fagerland from Tango (1998/1999) & Tang (2003)
#' This function is not vectorised
#'
#' @author Pete Laud, \email{p.j.laud@@sheffield.ac.uk}
#' @references
#'   Tango T. Equivalence test and confidence interval for the difference
#'   in proportions for the paired-sample design.
#'   Statistics in Medicine 1998; 17:891-908
#'
#'   Nam J-M, Blackwelder WC. Analysis of the ratio of marginal
#'   probabilities in a matched-pair setting.
#'   Stat Med 2002; 21(5):689–699
#'
#'   Tang N-S, Tang M-L, Chan ISF. On tests of equivalence via non-unity
#'   relative risk for matched-pair design.
#'   Statistics in Medicine 2003; 22:1217-1233
#'
#' @inheritParams pairbinci
#'
#' @noRd
scorepair <- function(theta,
                      x,
                      contrast = "RD",
                      cc = FALSE,
                      cctype = "constant",
                      ...) {
  N <- sum(x)
  if (as.character(cc) == "TRUE") cc <- 0.5

  if (contrast == "RD") {
    # notation per Tango 1999 letter
    Stheta <- ((x[2] - x[3]) - N * theta)
    A <- 2 * N
    B <- -x[2] - x[3] + (2 * N - x[2] + x[3]) * theta
    C_ <- -x[3] * theta * (1 - theta)
    num <- (-B + Re(sqrt(as.complex(B^2 - 4 * A * C_))))
    p2d <- ifelse(num == 0, 0, num / (2 * A))
    V <- pmax(0, N * (2 * p2d + theta * (1 - theta)))
    corr <- 2 * cc * sign(Stheta)
  }
  if (contrast == "RR") {
    # per Tang 2003
    Stheta <- ((x[2] + x[1]) - (x[3] + x[1]) * theta)
    A <- N * (1 + theta)
    B <- (x[1] + x[3]) * theta^2 - (x[1] + x[2] + 2 * x[3])
    C_ <- x[3] * (1 - theta) * (x[1] + x[2] + x[3]) / N
    num <- (-B + Re(sqrt(as.complex(B^2 - 4 * A * C_))))
    q21 <- ifelse(num == 0, 0, num / (2 * A))
    V <- pmax(0, N * (1 + theta) * q21 + (x[1] + x[2] + x[3]) * (theta - 1))
    if (cctype == "constant") corr <- cc * 2 * sign(Stheta)
    if (cctype == "delrocco") corr <- cc * (x[1] + x[3]) / N * sign(Stheta)
  }
  score <- ifelse(Stheta == 0, 0, (Stheta - corr) / sqrt(V))
  score[abs(Stheta) < abs(corr)] <- 0
  pval <- pnorm(score)
  outlist <- list(score = score, pval = pval)
  return(outlist)
}

#' Closed form Tango asymptotic score confidence intervals for a paired
#' difference of proportions (RD).
#'
#' R code to calculate Tango's score-based CI using a non-iterative method.
#' For contrast = "RD" only.
#' Code originates from Appendix B of Yang 2013,
#' with updates to include continuity correction from Chang 2024.
#'
#' @author Pete Laud, \email{p.j.laud@@sheffield.ac.uk}
#' @references
#'   Yang Z, Sun X and Hardin JW. A non-iterative implementation of Tango's
#'   score confidence interval for a paired difference of proportions.
#'   Statistics in Medicine 2013; 32:1336-1342
#'
#'   Chang P et al. Continuity corrected score confidence interval for the
#'   difference in proportions in paired data.
#'   Journal of Applied Statistics 2024; 51-1:139-152
#'
#' @inheritParams pairbinci
#'
#' @noRd
tangoci <- function(x,
                    level = 0.95,
                    cc = FALSE) {
  options(digits = 12)
  if (as.character(cc) == "TRUE") {
    cc <- 0.5 # Default correction for paired RD aligned with cc'd McNemar test
  }

  # output 2x2 table for validation
  x1i <- rep(c(1, 1, 0, 0), x)
  x2i <- rep(c(1, 0, 1, 0), x)
  xi <- table(x1i, x2i)

  N <- sum(x)
  b <- x[2]
  c <- x[3]
  alpha <- 1 - level
  g <- qnorm(1 - alpha / 2)^2

  # Updated to include continuity correction from Chang 2024
  keep1 <- keep2 <- NULL
  for (uplow in c(-1, 1)) {
    corr <- 2 * cc * uplow
    G <- N^2 + g * N
    H <- 0.5 * g * (2 * N - b + c) - 2 * N * (b - c + corr) - g * N
    I <- (b - c + corr)^2 - 0.5 * g * (b + c)
    m0 <- G^2
    m1 <- 2 * G * H
    m2 <- H^2 + 2 * G * I - 0.25 * g^2 * ((2 * N - b + c)^2 - 8 * N * c)
    m3 <- 2 * H * I - 0.25 * g^2 * (8 * N * c - 2 * (b + c) * (2 * N - b + c))
    m4 <- I^2 - 0.25 * g^2 * (b + c)^2

    u1 <- m1 / m0
    u2 <- m2 / m0
    u3 <- m3 / m0
    u4 <- m4 / m0

    if (b != c | cc > 0) {
      #  if (TRUE) {
      nu1 <- -u2
      nu2 <- u1 * u3 - 4 * u4
      nu3 <- -(u3^2 + u1^2 * u4 - 4 * u2 * u4)

      d1 <- nu2 - nu1^2 / 3
      d2 <- 2 * nu1^3 / 27 - nu1 * nu2 / 3 + nu3
      CritQ <- d2^2 / 4 + d1^3 / 27

      y1A <- h1 <- h2 <- h3 <- h4 <- root1 <- root2 <- root3 <- root4 <- c()
      if (CritQ > 0) { ## keep one real root;
        BigA <- -d2 / 2 + Re(sqrt(as.complex(CritQ)))
        BigB <- -d2 / 2 - Re(sqrt(as.complex(CritQ)))
        x1 <- sign(BigA) * abs(BigA)^(1 / 3) + sign(BigB) * abs(BigB)^(1 / 3)
        y1A <- x1 - nu1 / 3
      }

      if (CritQ == 0) { ## keep two real roots;
        BigA <- -d2 / 2 + Re(sqrt(as.complex(CritQ)))
        BigB <- -d2 / 2 - Re(sqrt(as.complex(CritQ)))
        Omega <- complex(real = -1 / 2, imaginary = sqrt(3) / 2)
        Omega2 <- complex(real = -1 / 2, imaginary = -sqrt(3) / 2)
        x1 <- sign(BigA) * abs(BigA)^(1 / 3) + sign(BigB) * abs(BigB)^(1 / 3)
        x2 <- Omega * sign(BigA) * abs(BigA)^(1 / 3) +
          Omega * sign(BigB) * abs(BigB)^(1 / 3)
        y1A[1] <- x1 - nu1 / 3
        y1A[2] <- x2 - nu1 / 3
      }

      if (CritQ < 0) { ## keep three real roots;
        BigA <- -d2 / 2 + sqrt(as.complex(CritQ))
        BigB <- -d2 / 2 - sqrt(as.complex(CritQ))
        Omega <- complex(real = -1 / 2, imaginary = sqrt(3) / 2)
        Omega2 <- complex(real = -1 / 2, imaginary = -sqrt(3) / 2)
        x1 <- BigA^(1 / 3) + BigB^(1 / 3)
        x2 <- Omega * BigA^(1 / 3) + Omega2 * BigB^(1 / 3)
        x3 <- Omega2 * BigA^(1 / 3) + Omega * BigB^(1 / 3)
        y1A[1] <- x1 - nu1 / 3
        y1A[2] <- x2 - nu1 / 3
        y1A[3] <- x3 - nu1 / 3
      }

      y1 <- Re(y1A) # keep the real part;
      ny <- length(y1)

      for (i in 1:ny) {
        h1[i] <- Re(sqrt(as.complex(u1^2 / 4 - u2 + y1[i])))
        h2[i] <- (u1 * y1[i] / 2 - u3) / (2 * h1[i])
        h3[i] <- (u1 / 2 + h1[i])^2 - 4 * (y1[i] / 2 + h2[i])
        h4[i] <- (h1[i] - u1 / 2)^2 - 4 * (y1[i] / 2 - h2[i])
        if (h3[i] >= 0) {
          root1[i] <- (-(u1 / 2 + h1[i]) + sqrt(h3[i])) / 2
          root2[i] <- (-(u1 / 2 + h1[i]) - sqrt(h3[i])) / 2
        }
        if (h4[i] >= 0) {
          root3[i] <- ((h1[i] - u1 / 2) + sqrt(h4[i])) / 2
          root4[i] <- ((h1[i] - u1 / 2) - sqrt(h4[i])) / 2
        }
      }

      lower <- max(-1, min(root1, root2, root3, root4, na.rm = TRUE))
      upper <- min(max(root1, root2, root3, root4, na.rm = TRUE), 1)

      if (b == N & c == 0) {
        root <- c(lower, 1)
      } else if (b == 0 & c == N) {
        root <- c(-1, upper)
      } else {
        root <- c(lower, upper)
      }
    }
    if (b == c & cc == 0) { # Only works without cc
      root1 <- -sqrt(-u2)
      root2 <- sqrt(-u2)
      root <- c(root1, root2)
    }
    if (uplow == -1) {
      keep1 <- root[1]
    } else if (uplow == 1) keep2 <- root[2]
  }
  if (b == c & b == N / 2) {
    # Rare special case fails to find solution to the quartic
    # and we have to resort to iterative method
    fixint <- pairbinci(x = x, contrast = "RD", method_RD = "Score", level = level, cc = cc)$estimates[c(1, 3)]
    keep1 <- fixint[1]
    keep2 <- fixint[2]
  }

  estimates <- cbind(
    Lower = keep1, MLE = (b - c) / N, Upper = keep2,
    level = level
  )

  return(estimates)
}

#' Closed form Tang asymptotic score confidence intervals for a paired
#' ratio of proportions (RR).
#'
#' R code to calculate Tang's score-based CI using a non-iterative method.
#' For contrast = "RR" only.
#' # Adapted from code kindly provided by Guogen Shan for the closed-form
#' ASCC method proposed in DelRocco et al. 2023.
#' with modified form of continuity correction
#' for consistency with McNemar test, and unified code with/without cc.
#' cctype = "constant" uses correction of cc * (2) with e.g. cc=0.5
#' instead of xp1 / (m * N) with m=2 as used by DelRocco
#'   - Note the latter produces a relatively much smaller correction
#'
#' @author Pete Laud, \email{p.j.laud@@sheffield.ac.uk}, Guogen Shan
#' @references
#'   DelRocco N et al. New Confidence Intervals for Relative Risk of Two
#'   Correlated Proportions.
#'   Statistics in Biosciences 2023; 15:1–30
#'
#' @inheritParams pairbinci
#'
#' @noRd
tangci <- function(x,
                   cc = FALSE,
                   cctype = "constant",
                   level = 0.95,
                   ctrl = 1e-10) {
  if (as.character(cc) == "TRUE") {
    cc <- 0.5
  }
  zeroflag <- (x[1] + x[2] == 0)
  infflag <- (x[1] + x[3] == 0)
  doublezero <- (x[1] + x[2] + x[3] == 0)

  if (any(x == 0)) {
    x <- x + ctrl # This avoids failure of the polynomial() function call later
  }
  x11 <- x[1]
  x12 <- x[2]
  x21 <- x[3]
  x22 <- x[4]

  xp1 <- x11 + x21
  x1p <- x11 + x12
  estimate <- x1p / xp1
  N <- sum(x)

  alpha <- 1 - level
  z <- qnorm(1 - alpha / 2)
  if (cctype == "delrocco") {
    # Version proposed by DelRocco et al
    corr <- cc * xp1 / N
  } else if (cctype == "constant") {
    # Alternative constant version
    corr <- cc * 2
  }
  if (doublezero) {
    lower <- 0
    upper <- Inf
  } else if (!doublezero) {
    a_lower <- xp1^4 + z^2 * xp1^3
    b_lower <- (-(2 * xp1^2 + z^2 * xp1) *
      (2 * xp1 * (x1p - corr) + z^2 * (x11 + x12 + x21)))
    c_lower <- (6 * xp1^2 * (x1p - corr)^2 +
      z^4 * (x1p + xp1) * (x11 + x12 + x21) +
      z^2 * xp1 * ((x1p - corr)^2 +
        4 * (x11 + x12 + x21) * (x1p - corr) +
        x1p * xp1))
    d_lower <- (-(2 * (x1p - corr)^2 + z^2 * x1p) *
      (2 * xp1 * (x1p - corr) + z^2 * (x11 + x12 + x21)))
    e_lower <- (x1p - corr)^4 + z^2 * x1p * (x1p - corr)^2

    b_lower <- b_lower / a_lower
    c_lower <- c_lower / a_lower
    d_lower <- d_lower / a_lower
    e_lower <- e_lower / a_lower
    a_lower <- a_lower / a_lower

    a_upper <- xp1^4 + z^2 * xp1^3
    b_upper <- (-(2 * xp1^2 + z^2 * xp1) *
      (2 * xp1 * (x1p + corr) + z^2 * (x11 + x12 + x21)))
    c_upper <- (6 * xp1^2 * (x1p + corr)^2 +
      z^4 * (x1p + xp1) * (x11 + x12 + x21) +
      z^2 * xp1 * ((x1p + corr)^2 +
        4 * (x11 + x12 + x21) * (x1p + corr) +
        x1p * xp1))
    d_upper <- (-(2 * (x1p + corr)^2 + z^2 * x1p) *
      (2 * xp1 * (x1p + corr) + z^2 * (x11 + x12 + x21)))
    e_upper <- (x1p + corr)^4 + z^2 * x1p * (x1p + corr)^2

    b_upper <- b_upper / a_upper
    c_upper <- c_upper / a_upper
    d_upper <- d_upper / a_upper
    e_upper <- e_upper / a_upper
    a_upper <- a_upper / a_upper

    if (!zeroflag) {
      lowertheta <- solve(polynom::polynomial(
        c(e_lower, d_lower, c_lower, b_lower, a_lower)
      ))
      lowertheta <- as.numeric(ifelse(Im(lowertheta) == 0, Re(lowertheta), NA))
      lower <- min(lowertheta, na.rm = TRUE)
    } else if (zeroflag) {
      estimate <- 0
      lower <- 0
    }
    if (!infflag) {
      uppertheta <- solve(polynom::polynomial(
        c(e_upper, d_upper, c_upper, b_upper, a_upper)
      ))
      uppertheta <- as.numeric(ifelse(Im(uppertheta) == 0, Re(uppertheta), NA))
      upper <- max(uppertheta, na.rm = TRUE)
    } else if (infflag) {
      estimate <- Inf
      upper <- Inf
    }
  }

  result <- cbind(Lower = lower, Estimate = estimate, Upper = upper)
  return(result)
}

#' MOVER interval for paired RR or RD (binomial only)
#'
#' Method of Variance Estimates Recovery, applied to paired RD and RR.
#' With various options for the marginal rates, and with optional continuity
#' correction, and Newcombe's correction to the Pearson correlation estimate,
#' applied to both contrasts.
#'
#' @author Pete Laud, \email{p.j.laud@@sheffield.ac.uk}
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
#'
#' @noRd
moverpair <- function(x,
                      level = 0.95,
                      contrast = "RR",
                      method = "SCAS",
                      cc = FALSE,
                      corc = FALSE) {
  if (!is.numeric(c(x))) {
    print("Non-numeric inputs!")
    stop()
  }

  x1p <- x[1] + x[2]
  xp1 <- x[1] + x[3]
  N <- sum(x)
  z <- qnorm(1 - (1 - level) / 2)
  ## First calculate l1, u1, l2, u2 based on an interval for p1p and pp1
  if (method == "SCAS") {
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
  p1phat <- j1[, 2]
  pp1hat <- j2[, 2]

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

  estimates <- cbind(Lower = lower, Estimate = estimate, Upper = upper)
  row.names(estimates) <- NULL
  return(estimates)
}
