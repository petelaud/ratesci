#' Confidence intervals for comparisons of paired binomial rates.
#'
#' Confidence intervals for the rate (or risk) difference ("RD"), rate ratio
#' ("RR") or conditional odds ratio ("OR"), for paired binomial data. (For
#' paired Poisson rates, suggest use the tdasci function with `distrib = "poi"`,
#' and `weighting = "MH"`, with pairs as strata.)
#' This function applies the score-based Tango and Tang methods for RD and
#' RR respectively, with iterative and closed-form versions, and an added
#' skewness correction for improved one-sided coverage.
#' Also includes MOVER options using the Method of Variance Estimates Recovery
#' for paired RD and RR, incorporating Newcombe's correlation correction, and
#' some simpler methods by Bonett & Price for RD and RR.
#' For OR, intervals are produced based on transforming various intervals for
#' the single proportion, including SCASp, mid-p and Jeffreys.
#' All methods have options for continuity adjustment, and the magnitude of
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
#' @param level Number specifying confidence level (between 0 and 1, default
#'   0.95).
#' @param contrast Character string indicating the contrast of interest: \cr
#'   "RD" = rate difference (default); \cr
#'   "RR" = rate ratio; \cr
#'   "OR" = conditional odds ratio.
#' @param method Character string indicating the confidence interval method
#'   to be used. The following are available for `contrast = "RD"` or `"RR"`: \cr
#'   "Score" = (default) asymptotic score class of methods including Tango
#'             (for RD) / Tang (for RR), by iterative calculations, with
#'             optional skewness correction; \cr
#'   "Score_closed" = closed form solution for Tango/Tang intervals
#'                    (without skewness correction); \cr
#'   "MOVER" = hybrid MOVER method (as per "method 8" in Newcombe, but with
#'             a choice of input methods - see moverbase); \cr
#'   "MOVER_newc" = hybrid MOVER methods with correction to correlation
#'                  estimate (Newcombe's "method 10"); \cr
#'   "TDAS" = t-distribution asymptotic score (experimental method, now
#'            deprecated); \cr
#'   "BP" = Wald with Bonett-Price adjustment for RD, or Hybrid Bonett-Price
#'          method for RR. \cr
#'   For `contrast = "OR"`, one of the following methods may be selected,
#'   all of which are based on transformation of an interval for a single
#'   proportion `b/(b+c)`: \cr
#'   "SCASp" = transformed skewness-corrected score (default); \cr
#'   "jeff" = transformed Jeffreys; \cr
#'   "midp" = transformed mid-p; \cr
#'   "wilson" = transformed Wilson score - included for reference only, not
#'   recommended.
#' @param moverbase Character string indicating the base method used as input
#'   for the MOVER methods for RD or RR (when method = "MOVER" or "MOVER_newc"),
#'   and for the Hybrid BP method for RR:
#'   "jeff" = Jeffreys equal-tailed interval (default),
#'   "SCASp" = skewness-corrected score,
#'   "midp" = mid-p,
#'   "wilson" = Wilson score (not recommended, known to be skewed).
#' @param bcf Logical (default FALSE) indicating whether to apply variance bias
#'   correction in the score denominator. (Under evaluation, manuscript under
#'   review.)
#' @param skew Logical (default TRUE) indicating whether to apply skewness
#'   correction or not. (Under evaluation, manuscript under review.)
#'   - Only applies for the iterative `method = "Score"`.
#' @param cc Number or logical (default FALSE) specifying (amount of) continuity
#'   adjustment. When a score-based method is used, cc = 0.5 corresponds to the
#'   continuity-corrected McNemar test.
#' @param cctype (deprecated: new equivariant cc method implemented instead.)
#' @param theta0 Number to be used in a one-sided significance test (e.g.
#'   non-inferiority margin). 1-sided p-value will be < 0.025 iff 2-sided 95\% CI
#'   excludes theta0. NB: can also be used for a superiority test by setting
#'   theta0 = 0.
#' @param method_RD (deprecated: parameter renamed to method)
#' @param method_RR (deprecated: parameter renamed to method)
#' @param method_OR (deprecated: parameter renamed to method)
#' @param precis Number (default 6) specifying precision (i.e. number of decimal
#'   places) to be used in optimisation subroutine for the confidence interval.
#' @param warn Logical (default TRUE) giving the option to suppress warnings.
#' @param ... Other arguments.
#' @importFrom stats uniroot pbinom ppois dpois
#' @return A list containing the following components: \describe{
#'   \item{data}{the input data in 2x2 matrix form.}
#'   \item{estimates}{the requested contrast, with its confidence interval and
#'   the specified confidence level, along with estimates of the marginal
#'   probabilities and the correlation coefficient (uncorrected and
#'   corrected).}
#'   \item{pval}{the corresponding 2-sided significance test
#'   against the null hypothesis that p_1 = p_2, and one-sided
#'   significance tests against the null hypothesis that theta >= or <= theta0
#'   as specified.}
#'   \item{call}{details of the function call.}}
#' @examples
#' # Example from Fagerland et al 2014
#' # SCAS method for RD
#' pairbinci(x = c(1, 1, 7, 12), contrast = "RD", method = "Score")
#' # Tango method
#' pairbinci(x = c(1, 1, 7, 12), contrast = "RD", method = "Score", skew = FALSE, bcf = FALSE)
#' # MOVER-NJ method
#' pairbinci(x = c(1, 1, 7, 12), contrast = "RD", method = "MOVER_newc", moverbase = "jeff")
#' # SCAS for RR
#' pairbinci(x = c(1, 1, 7, 12), contrast = "RR", method = "Score")
#' # Tang method
#' pairbinci(x = c(1, 1, 7, 12), contrast = "RR", method = "Score", skew = FALSE, bcf = FALSE)
#' # MOVER-NJ
#' pairbinci(x = c(1, 1, 7, 12), contrast = "RR", method = "MOVER_newc", moverbase = "jeff")
#' # Transformed SCASp method for OR
#' pairbinci(x = c(1, 1, 7, 12), contrast = "OR", method = "SCASp")
#' # Transformed Wilson method
#' pairbinci(x = c(1, 1, 7, 12), contrast = "OR", method = "wilson")
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
#'   Bonett DG, Price RM. Confidence intervals for a ratio of binomial
#'   proportions based on paired data.
#'   Statistics in Medicine 2006; 25:3039-3047
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
#'   Laud PJ. Comments on "New Confidence Intervals for Relative Risk of Two
#'   Correlated Proportions" (2023).
#'   Statistics in Biosciences 2025; https://doi.org/10.1007/s12561-025-09479-4
#'
#'   Laud PJ. Improved confidence intervals and tests for paired binomial
#'   proportions. (2025, Under review)
#'
#' @export
pairbinci <- function(x,
                      level = 0.95,
                      contrast = "RD",
                      method = ifelse(contrast == "OR", "SCASp", "Score"),
                      moverbase = ifelse(method %in% c("MOVER", "MOVER_newc", "BP"), "jeff", NULL),
                      bcf = TRUE,
                      skew = TRUE,
                      cc = FALSE,
                      theta0 = NULL,
                      precis = 6,
                      warn = TRUE,
                      method_RD = NULL,
                      method_RR = NULL,
                      method_OR = NULL,
                      cctype = NULL,
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
  if (!is.null(cctype)) {
    warning(
      "argument cctype is deprecated; previous options are replaced with a single
       cc method that produces an equivariant interval.",
      call. = FALSE
    )
  }
  if (!is.null(method_RD)) {
    warning(
      "argument method_RD is deprecated; please use method instead.",
      call. = FALSE
    )
    method <- method_RD
  }
  if (!is.null(method_RR)) {
    warning(
      "argument method_RR is deprecated; please use method instead.",
      call. = FALSE
    )
    method <- method_RR
  }
  if (!is.null(method_OR)) {
    warning(
      "argument method_OR is deprecated; please use method instead.",
      call. = FALSE
    )
    method <- method_OR
  }
  if (!(tolower(substr(contrast, 1, 2)) %in% c("rd", "rr", "or"))) {
    print("Contrast must be one of 'RD', 'RR' or 'OR'")
    stop()
  }
  if (contrast %in% c("RD", "RR")) {
    if (tolower(method) == "tdas") {
      print("TDAS method is deprecated in pairbinci(); if needed use tdasci() function,
             but method = 'Score' with skew = TRUE is recommended instead.")
      stop()
    }
    if (!(tolower(substr(method, 1, 4)) %in%
      c("scor", "move", "bp"))) {
      print("Method must be one of 'Score_closed', 'Score', 'BP',
          'MOVER' or 'MOVER_newc' for contrast = 'RD' or 'RR'")
      stop()
    }
    # if (FALSE) {
    if (skew == TRUE &&
      (tolower(substr(method, 1, 12)) == "score_closed")) {
      method <- "Score"
      if (warn == TRUE) {
        print(paste("Closed-form calculation not available with skewness correction -
         method is set to 'Score' instead"))
      }
    }
    # }
    if (method %in% c("MOVER", "MOVER_newc", "BP") &&
      !(tolower(substr(moverbase, 1, 4)) %in%
        c("scas", "wils", "midp", "jeff"))) {
      print("moverbase must be one of 'SCASp', 'wilson', 'midp' or 'jeff'")
      stop()
    }
  }
  if (contrast == "OR") {
    if (!(tolower(substr(method, 1, 4)) %in%
      c("scas", "wils", "midp", "jeff"))) {
      print("Method must be one of 'SCASp', 'wilson', 'midp' or 'jeff' for
              contrast = 'OR'")
      stop()
    }
  }
  if (method == "BP" && contrast == "RD" && cc != FALSE) {
    cc <- FALSE
    if (warn == TRUE) {
      print(paste("Continuity adjustment not available for
         Wald with Bonett-Price adjustment for RD -
         cc is set to FALSE"))
    }
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
  xi <- table(Test_1 = factor(x1i, levels = c("Success", "Failure")),
              Test_2 = factor(x2i, levels = c("Success", "Failure")))
  x1 <- x[1] + x[2]
  x2 <- x[1] + x[3]
  N <- sum(x)
  p1hat <- x1 / N
  p2hat <- x2 / N
  phi_hat <- (x[1] * x[4] - x[2] * x[3]) / sqrt(x1 * (N - x1) * x2 * (N - x2))
  #if (is.na(phi_hat) | is.infinite(phi_hat)) {
  if (is.na(phi_hat) ) {
      phi_hat <- 0
  }
  phi_c <- phi_hat
  if (x[1] * x[4] - x[2] * x[3] > 0) {
    phi_c <- (max(x[1] * x[4] - x[2] * x[3] - N / 2, 0) /
      sqrt(x1 * (N - x1) * x2 * (N - x2)))
  }
  psi_hat <- x[1] * x[4] / (x[2] * x[3])
  if (is.na(psi_hat)) psi_hat <- 0


  if (contrast == "OR") {
    # special case for OR, use conditional method based on transforming the
    # SCAS interval for a single proportion
    # Transformed SCAS method with bias correction factor based on
    # N/(N-1) i.e. full sample size, not the number of discordant pairs
    # This gives a p-value matching that for other bias-corrected methods
    x12 <- x[2]
    x21 <- x[3]
    if (method == "SCASp") {
      trans_th0 <- NULL
      if (is.null(theta0)) theta0 <- 1
      trans_th0 <- theta0 / (1 + theta0)
      OR_ci <- scaspci(
        x = x12, n = x12 + x21, distrib = "bin",
        level = level, cc = cc, bcf = bcf, bign = N
      )$estimates[, c(1:3), drop = FALSE]
      estimates <- OR_ci / (1 - OR_ci)
      scorezero <- scoretheta(
        theta = 0.5, x1 = x12, n1 = x12 + x21, n2 = x[1] + x[4],
        contrast = "p", distrib = "bin", bcf = bcf,
        skew = TRUE, cc = cc, level = level
      )
      chisq_zero <- scorezero$score^2
      pval2sided <- pchisq(chisq_zero, 1, lower.tail = FALSE)
      scoreth0 <- scoretheta(
        theta = trans_th0, x1 = x12, n1 = x12 + x21, n2 = x[1] + x[4],
        contrast = "p", distrib = "bin", bcf = bcf,
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
    } else if (method == "midp") {
      trans_ci <- exactci(x = x12, n = x12 + x21, midp = 0.5 - cc, level = level)
      estimates <- (trans_ci / (1 - trans_ci))
      outlist <- list(xi, estimates = estimates)
    } else if (method == "wilson") {
      trans_ci <- wilsonci(x = x12, n = x12 + x21, cc = cc, level = level)
      estimates <- (trans_ci / (1 - trans_ci))
      outlist <- list(xi, estimates = estimates)
    } else if (method == "jeff") {
      trans_ci <- jeffreysci(x = x12, n = x12 + x21, cc = cc, level = level)$estimates[, c(1:3), drop = FALSE]
      estimates <- (trans_ci / (1 - trans_ci))
      outlist <- list(xi, estimates = estimates)
    }
  } else if (contrast != "OR") {
    # Iterative Score methods by Tango (for RD) & Tang (for RR):
    if (method == "Score") {
      myfun <- function(theta) {
        scorepair(
          theta = theta, x = x, contrast = contrast, cc = cc,
          skew = skew, bcf = bcf
        )$score
      }
      myfun0 <- function(theta) {
        scorepair(
          theta = theta, x = x, contrast = contrast, cc = 0,
          skew = skew, bcf = bcf
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
      at_MLE <- scorepair(
        theta = MLE, x = x, contrast = contrast, cc = cc,
        skew = skew, bcf = bcf
      )
      # fix some extreme cases with zero counts
      if (contrast %in% c("RR")) {
        upper[x2 == 0] <- Inf
        MLE[x2 == 0 & skew == FALSE] <- Inf
      }
      estimates <- cbind(
        lower = lower, est = MLE, upper = upper,
        level = level, p1hat = p1hat, p2hat = p2hat,
        p1mle = at_MLE$p1d, p2mle = at_MLE$p2d,
        phi_hat = phi_hat, phi_c = phi_c, psi_hat = psi_hat
      )
      # Closed-form versions of Score methods by Chang (for RD Tango method)
      #  & DelRocco (for RR Tang method):
    } else if (contrast == "RD" && method == "Score_closed") {
      estimates <- tangoci(x = x, level = level, cc = cc, bcf = bcf)
    } else if (contrast == "RR" && method == "Score_closed") {
      estimates <- tangci(
        x = x, level = level, cc = cc, bcf = bcf
      )
    }
    # MOVER methods for RD and RR
    if (method == "MOVER") {
      estimates <- moverpair(
        x = x, contrast = contrast, level = level,
        method = moverbase, cc = cc, corc = FALSE
      )
      outlist <- list(data = xi, estimates = estimates)
    }
    if (method == "MOVER_newc") {
      estimates <- moverpair(
        x = x, contrast = contrast, level = level,
        method = moverbase, cc = cc, corc = TRUE
      )
      outlist <- list(data = xi, estimates = estimates)
    }
    if (contrast == "RR" && tolower(method) == "bp") {
      estimates <- bpci(
        x = x, contrast = contrast, level = level, cc = cc, method = moverbase
      )
      outlist <- list(data = xi, estimates = estimates)
    }
    if (contrast == "RD" && tolower(method) == "bp") {
      estimates <- bpci(
        x = x, contrast = contrast, level = level
      )
      outlist <- list(data = xi, estimates = estimates)
    }

    if ((method == "Score") || (method == "Score_closed")) {
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
        cc = cc, skew = skew, bcf = bcf
      )
      scorenull <- scorepair(
        theta = theta0, x = x, contrast = contrast,
        cc = cc, skew = skew, bcf = bcf
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
  outlist <- list(data = xi, estimates = round(estimates, precis))
  # MOVER methods don't produce p-values
  if (!(method %in% c("MOVER", "MOVER_newc", "midp", "wilson", "jeff", "BP"))) {
    outlist <- append(outlist, list(pval = pval))
  }
  # Set unused arguments to null to omit them from call
  if (!(method %in% c("MOVER", "MOVER_newc") ||
        (method == "BP" && contrast == "RR"))) {
    moverbase <- NULL
  }
  call <- c(
    contrast = contrast, method = method, moverbase = moverbase,
    level = level, bcf = bcf, skew = skew, cc = cc
  )
  outlist <- append(outlist, list(call = call))
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
                      skew = TRUE,
                      bcf = TRUE,
                      cc = FALSE,
                      ...) {
  N <- sum(x)
  lambda <- switch(as.character(bcf),
    "TRUE" = N / (N - 1),
    "FALSE" = 1
  )
  if (as.character(cc) == "TRUE") cc <- 0.5

  if (contrast == "RD") {
    # notation per Tango 1999 letter, divided by N
    # (to get numerator on the right scale for skewness correction)
    # Variance is divided by N^2 accordingly
    # and continuity adjustment term also divided by N
    Stheta <- ((x[2] - x[3]) - N * theta) / N
    A <- 2 * N
    B <- -x[2] - x[3] + (2 * N - x[2] + x[3]) * theta
    C_ <- -x[3] * theta * (1 - theta)
    num <- (-B + Re(sqrt(as.complex(B^2 - 4 * A * C_))))
    p21 <- ifelse(num == 0, 0, num / (2 * A))
    corr <- 2 * cc * sign(Stheta) / N
    p12 <- p21 + theta
    # From Tango
    p11 <- ifelse(x[1] == 0, 0, x[1] / (x[1] + x[4]) * (1 - p12 - p21))
    #    p22 <- ifelse(x[4] == 0, 0, x[4] / (x[1] + x[4]) * (1 - p12 - p21))
    p22 <- 1 - p11 - p12 - p21
    p2d <- pmin(1, pmax(0, p21 + p11))
    #    p1d <- pmin(1, pmax(0, p12 + p11))
    p1d <- p2d + theta

    # Tango variance
    #    V <- pmax(0, N * (2 * p21 + theta * (1 - theta))) * lambda / (N^2)
    # Equivalent
    V <- pmax(0, (p1d * (1 - p1d) + p2d * (1 - p2d) -
      2 * (p11 * p22 - p12 * p21)) / N) * lambda
    mu3 <- (p1d * (1 - p1d) * (1 - 2 * p1d) +
      ((-1)^3) * p2d * (1 - p2d) * (1 - 2 * p2d) +
      3 * (-1) * (p11 * (1 - p1d)^2 + p21 * p1d^2 - p1d * p2d * (1 - p1d)) +
      3 * ((-1)^2) * (p11 * (1 - p2d)^2 + p12 * p2d^2 - p1d * p2d * (1 - p2d))
    ) / (N^2)
  }
  if (contrast == "RR") {
    # per Tang 2003, but divided by N
    Stheta <- ((x[2] + x[1]) - (x[3] + x[1]) * theta) / N
    A <- N * (1 + theta)
    B <- (x[1] + x[3]) * theta^2 - (x[1] + x[2] + 2 * x[3])
    C_ <- x[3] * (1 - theta) * (x[1] + x[2] + x[3]) / N
    num <- (-B + Re(sqrt(as.complex(B^2 - 4 * A * C_))))
    q21 <- ifelse(num == 0, 0, num / (2 * A))

    # Equivariant continuity adjustment for RR, aligned with McNemar cc.
    # replaces previous 'delrocco' and 'constant' cctype options
      corr <- cc * (1 + theta) * sign(Stheta) / N

    q12num <- (q21 + (theta - 1) * (1 - x[4] / N))
    q12 <- ifelse(q12num < 1E-10, 0, q12num / theta)
    # Below from Tang 2003
    q11num <- (1 - x[4] / N - (1 + theta) * q21)
    q11 <- ifelse(q11num < 1E-10, 0, q11num / theta)
    q22 <- 1 - q11 - q12 - q21
    p2d <- q21 + q11
    #    p1d <- q12 + q11
    p1d <- p2d * theta

    # Tang variance
    #    V <- pmax(0, N * (1 + theta) * q21 + (x[1] + x[2] + x[3]) * (theta - 1)) * lambda / (N^2)
    # Nam-Blackwelder variance
    #    V <- pmax(0, theta*(q12 + q21) / N) * lambda
    # Equivalent consistent with scoreci notation
    V <- pmax(0, (p1d * (1 - p1d) + theta^2 * p2d * (1 - p2d) -
      2 * theta * (q11 * q22 - q12 * q21)) / N) * lambda
    mu3 <- (p1d * (1 - p1d) * (1 - 2 * p1d) +
      ((-theta)^3) * p2d * (1 - p2d) * (1 - 2 * p2d) +
      3 * (-theta) * (q11 * (1 - p1d)^2 + q21 * p1d^2 - p1d * p2d * (1 - p1d)) +
      3 * ((-theta)^2) * (q11 * (1 - p2d)^2 + q12 * p2d^2 - p1d * p2d * (1 - p2d))) / (N^2)
  }
  scterm <- mu3 / (6 * V^(3 / 2))
  scterm[abs(mu3) < 1E-10 | abs(V) < 1E-10] <- 0 # Avoids issues with e.g. x = c(1, 3, 0, 6)
  score1 <- ifelse(Stheta == 0, 0, (Stheta - corr) / sqrt(V))
  A <- scterm
  B <- 1
  C_ <- -(score1 + scterm)
  num <- (-B + sqrt(pmax(0, B^2 - 4 * A * C_)))
  score <- ifelse((skew == FALSE | scterm == 0),
    score1, num / (2 * A)
  )
  score[abs(Stheta) < abs(corr)] <- 0

  pval <- pnorm(score)
  outlist <- list(
    score = score, p1d = p1d, p2d = p2d,
    mu3 = mu3, pval = pval
  )
  return(outlist)
}


#' Closed form Tango asymptotic score confidence intervals for a paired
#' difference of proportions (RD).
#'
#' R code to calculate Tango's score-based CI using a non-iterative method.
#' For contrast = "RD" only.
#' Code originates from Appendix B of Yang 2013,
#' with updates to include continuity adjustment from Chang 2024.
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
                    bcf = TRUE,
                    cc = FALSE) {
  if (as.character(cc) == "TRUE") {
    cc <- 0.5 # Default correction for paired RD aligned with cc'd McNemar test
  }

  # output 2x2 table for validation
  x1i <- rep(c(1, 1, 0, 0), x)
  x2i <- rep(c(1, 0, 1, 0), x)
  xi <- table(x1i, x2i)

  N <- sum(x)
  lambda <- switch(as.character(bcf),
    "TRUE" = N / (N - 1),
    "FALSE" = 1
  )

  x12 <- x[2]
  x21 <- x[3]
  alpha <- 1 - level
  g <- qnorm(1 - alpha / 2)^2 * lambda

  # Updated to include continuity adjustment from Chang 2024
  keep1 <- keep2 <- NULL
  for (uplow in c(-1, 1)) {
    corr <- 2 * cc * uplow
    G <- N^2 + g * N
    H <- 0.5 * g * (2 * N - x12 + x21) - 2 * N * (x12 - x21 + corr) - g * N
    I <- (x12 - x21 + corr)^2 - 0.5 * g * (x12 + x21)
    m0 <- G^2
    m1 <- 2 * G * H
    m2 <- H^2 + 2 * G * I - 0.25 * g^2 * ((2 * N - x12 + x21)^2 - 8 * N * x21)
    m3 <- 2 * H * I - 0.25 * g^2 * (8 * N * x21 - 2 * (x12 + x21) * (2 * N - x12 + x21))
    m4 <- I^2 - 0.25 * g^2 * (x12 + x21)^2

    u1 <- m1 / m0
    u2 <- m2 / m0
    u3 <- m3 / m0
    u4 <- m4 / m0

    if (u1 == 0 & u3 == 0) {
      # Special case with x12 = x21 = N/2 reduces to a quadratic equation
      #    x^4 + u2 x^2 + u4 = 0 reduces to
      #    x^2 = sqrt((-u2 + sqrt(u2^2 - 4*u4))/2)
      keep2 <- sqrt((-u2 + sqrt(u2^2 - 4 * u4)) / 2)
      keep1 <- -keep2
    } else if (x12 == x21 & cc == 0) {
      u2 <- -(((x12 + x21) / g + 1) / (N / g + 1)^2)
      root1 <- -sqrt(-u2)
      root2 <- sqrt(-u2)
      root <- c(root1, root2)
    } else if (cc > 0 & x12 == x21 & x12 == N / 2) {
      # Rare special case fails to find solution to the quartic
      # and we have to resort to iterative method
      root <- pairbinci(x = x, contrast = "RD", method = "Score", level = level, cc = cc, bcf = bcf, skew = FALSE)$estimates[c(1, 3)]
      # Alternative method using polynomial()
      #    lowertheta <- solve(polynom::polynomial(c(u4, u3, u2, u1, 1)))
      #    root1 <- min(as.numeric(ifelse(Im(lowertheta) == 0, Re(lowertheta), NA)), na.rm = TRUE)
      #    root2 <- max(as.numeric(ifelse(Im(lowertheta) == 0, Re(lowertheta), NA)), na.rm = TRUE)
    } else if (x12 != x21 | (cc > 0 & !(x12 == x21 & x12 == N / 2))) {
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

      if (x12 == N & x21 == 0) {
        root <- c(lower, 1)
      } else if (x12 == 0 & x21 == N) {
        root <- c(-1, upper)
      } else {
        root <- c(lower, upper)
      }
    }
    if (x12 == x21 & cc == 0) { # Only works without cc
      root1 <- -sqrt(-u2)
      root2 <- sqrt(-u2)
      root <- c(root1, root2)
    }
    if (uplow == -1) {
      keep1 <- root[1]
    } else if (uplow == 1) keep2 <- root[2]
  }

  estimates <- cbind(
    lower = keep1, est = (x12 - x21) / N, upper = keep2,
    level = level
  )

  return(estimates)
}

#' Closed form Tang asymptotic score confidence intervals for a paired
#' ratio of proportions (RR).
#'
#' R code to calculate Tang's score-based CI using a non-iterative method.
#' without dependence on polynom package and without needing ctrl argument.
#' For contrast = "RR" only.
#' This could be combined with tangoci to avoid code repetition.
#' Adapted from code kindly provided by Guogen Shan for the closed-form
#' ASCC method proposed in DelRocco et al. 2023.
#' with modified form of continuity adjustment (Laud 2025)
#' for consistency with McNemar test, and unified code with/without cc.
#'
#' @author Pete Laud, \email{p.j.laud@@sheffield.ac.uk}
#' @references
#'   DelRocco N et al. New Confidence Intervals for Relative Risk of Two
#'   Correlated Proportions.
#'   Statistics in Biosciences 2023; 15:1–30
#'
#'   Laud PJ. Comments on "New Confidence Intervals for Relative Risk of Two
#'   Correlated Proportions" (2023).
#'   Statistics in Biosciences 2025; https://doi.org/10.1007/s12561-025-09479-4
#'
#' @inheritParams pairbinci
#'
#' @noRd
tangci <- function(x,
                   level = 0.95,
                   bcf = TRUE,
                   cc = FALSE) {
  if (as.character(cc) == "TRUE") {
    cc <- 0.5
  }
  zeroflag <- (x[1] + x[2] == 0)
  infflag <- (x[1] + x[3] == 0)
  doublezero <- (x[1] + x[2] + x[3] == 0)

  N <- sum(x)
  x11 <- x[1]
  x12 <- x[2]
  x21 <- x[3]
  x22 <- x[4]

  xp1 <- x11 + x21
  x1p <- x11 + x12
  estimate <- x1p / xp1
  lambda <- switch(as.character(bcf),
    "TRUE" = N / (N - 1),
    "FALSE" = 1
  )

  alpha <- 1 - level
  z <- qnorm(1 - alpha / 2) * sqrt(lambda)
  corr2 <- 0
  # Previous cctype options replaced with an equivariant method
  corr1 <- corr2 <- cc
  if (doublezero) {
    lower <- 0
    upper <- Inf
  } else if (!doublezero) {
    for (uplow in c(-1, 1)) {
      cor1 <- corr1 * uplow
      cor2 <- corr2 * uplow

      m0 <- (xp1 - cor2)^4 + z^2 * xp1 * (xp1 - cor2)^2
      m1 <- (-(2 * (xp1 - cor2)^2 + z^2 * xp1) *
        (2 * (xp1 - cor2) * (x1p + cor1) + z^2 * (x11 + x12 + x21)))
      m2 <- (6 * (xp1 - cor2)^2 * (x1p + cor1)^2 +
        z^4 * (x1p + xp1) * (x11 + x12 + x21) +
        z^2 * (xp1 * (x1p + cor1)^2 +
          4 * (x11 + x12 + x21) * (x1p + cor1) * (xp1 - cor2) +
          x1p * (xp1 - cor2)^2))
      m3 <- (-(2 * (x1p + cor1)^2 + z^2 * x1p) *
        (2 * (xp1 - cor2) * (x1p + cor1) + z^2 * (x11 + x12 + x21)))
      m4 <- (x1p + cor1)^4 + z^2 * x1p * (x1p + cor1)^2

      y1A <- h1 <- h2 <- h3 <- h4 <- root1 <- root2 <- root3 <- root4 <- c()
      if (m0 == 0 & m1 == 0) {
        # Solve quadratic equation instead
        root1 <- quadroot(m2, m3, m4)
      } else {
        u1 <- m1 / m0
        u2 <- m2 / m0
        u3 <- m3 / m0
        u4 <- m4 / m0

        nu1 <- -u2
        nu2 <- u1 * u3 - 4 * u4
        nu3 <- -(u3^2 + u1^2 * u4 - 4 * u2 * u4)

        d1 <- nu2 - nu1^2 / 3
        d2 <- 2 * nu1^3 / 27 - nu1 * nu2 / 3 + nu3
        CritQ <- d2^2 / 4 + d1^3 / 27

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
      }
      if (uplow == -1) {
        lower <- ifelse(zeroflag, 0,
          max(0, min(root1, root2, root3, root4, na.rm = TRUE))
        )
      } else if (uplow == 1) {
        upper <- ifelse(infflag, Inf,
          max(root1, root2, root3, root4, na.rm = TRUE)
        )
      }
    }
  }

  result <- cbind(lower = lower, est = estimate, upper = upper)
  return(result)
}

#' MOVER interval for paired RR or RD (binomial only)
#'
#' Method of Variance Estimates Recovery, applied to paired RD and RR.
#' With various options for the marginal rates, and with optional continuity
#' adjustment, and Newcombe's correction to the Pearson correlation estimate,
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


#' Bonett-Price confidence intervals for a paired difference (RD) or ratio (RR)
#'
#' R code to calculate Bonett-Price's hybrid score-based CI for contrast = "RR",
#' and adjusted Wald-based method by the same authors for contrast = "RD".
#'
#' Code is adapted from Appendix A of Bonett & Price 2006
#'
#' @author Pete Laud, \email{p.j.laud@@sheffield.ac.uk}
#' @references
#'   Bonett DG, Price RM. Confidence intervals for a ratio of binomial
#'   proportions based on paired data.
#'   Statistics in Medicine 2006; 25:3039-3047
#'
#'   Bonett DG, Price RM. Adjusted Wald Confidence Interval for a Difference
#'   of Binomial Proportions Based on Paired Data.
#'   Journal of Educational and Behavioral Statistics 2012; 37(4):479-488
#'
#' @inheritParams pairbinci
#'
#' @noRd
bpci <- function(x,
                 level = 0.95,
                 contrast = "RD",
                 method = "wilson",
                 cc = FALSE) {
  x11 <- x[1]
  x10 <- x[2]
  x01 <- x[3]
  z0 <- qnorm(1 - (1 - level) / 2)
  x1 <- x11 + x10
  x0 <- x11 + x01

  if (contrast == "RD") {
    n <- sum(x)
    p12 <- (x10 + 1) / (n + 2)
    p21 <- (x01 + 1) / (n + 2)
    estimate <- p12 - p21
    v <- (p12 + p21 - (p12 - p21)^2) / (n + 2)
    estimates <- cbind(
      lower = pmax(-1, estimate - z0 * sqrt(v)),
      Estimate = estimate,
      upper = pmin(1, estimate + z0 * sqrt(v))
    )
  } else if (contrast == "RR") {
    estimate <- x1 / x0
    n <- x11 + x10 + x01
    s1 <- sqrt((1 - (x1 + 1) / (n + 2)) / (x1 + 1))
    s0 <- sqrt((1 - (x0 + 1) / (n + 2)) / (x0 + 1))
    s2 <- sqrt((x10 + x01 + 2) / ((x1 + 1) * (x0 + 1)))
    k <- s2 / (s1 + s0)
    z <- k * z0
    adjlevel <- 1 - 2 * (1 - pnorm(z))
    if (method == "wilson") {
      j1 <- wilsonci(x = x1, n = n, distrib = "bin", level = adjlevel, cc = cc)
      j2 <- wilsonci(x = x0, n = n, distrib = "bin", level = adjlevel, cc = cc)
    } else if (method == "jeff") {
      j1 <- rateci(x = x1, n = n, distrib = "bin", level = adjlevel, cc = cc)$jeff
      j2 <- rateci(x = x0, n = n, distrib = "bin", level = adjlevel, cc = cc)$jeff
    }
    ul <- j1[, 3] / j2[, 1]
    ll <- j1[, 1] / j2[, 3]
    estimates <- cbind(lower = ll, est = estimate, upper = ul)
    row.names(estimates) <- NULL
  }
  estimates
}
