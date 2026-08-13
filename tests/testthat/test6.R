# Tests relating to Vincent Jaquet bug report July 2026
# whereby evaluation of the stratified score can produce a spurious zero result
# resolved by using an alternative 'citardauq' formula for quadratic solution
# as per https://math.stackexchange.com/questions/866331/numerically-stable-algorithm-for-solving-the-quadratic-equation-when-a-is-very

test_that("point estimate not too far from crude difference", {

  x1 <- c(5, 9)
  n1 <- c(10, 11)
  x2 <- c(0, 1)
  n2 <- c(6, 5)
  res <- ratesci::scoreci(x1 = x1,
                  n1 = n1,
                  x2 = x2,
                  n2 = n2,
                  contrast = "RD",
                  stratified = TRUE,
                  cc = FALSE,
                  skew = TRUE,
                  weighting = "MH")

  mle <- res$estimates[2]
  crude <- sum(res$stratdata[, "theta_j"] * res$stratdata[, "wtpct_fixed"]) /100
  expect_true(
    !(abs(mle - crude) > 0.01 & abs(mle) < 0.0001)
  )

  if(FALSE) {
    # One-off check that unstratified intervals are not affected by the issue
    maxn1 <- 25
    maxn2 <- 25
    ns <- expand.grid(10:maxn1, 10:maxn2)
    n1s <- ns[ , 1]
    n2s <- ns[ , 2]

    x1 <- x2 <- n1 <- n2 <- c()

    for (i in 1:length(n1s)) {
      n1i <- n1s[i]
      n2i <- n2s[i]
      xs <- expand.grid(0:n1i, 0:n2i)
      x1i <- xs[ , 1]
      x2i <- xs[ , 2]

      x1 <- c(x1, x1i)
      x2 <- c(x2, x2i)
      n1 <- c(n1, rep(n1i, length(x1i)))
      n2 <- c(n2, rep(n2i, length(x2i)))

    }


    res <- ratesci::scoreci(x1 = x1,
                            n1 = n1,
                            x2 = x2,
                            n2 = n2,
                            contrast = "RD",
                            stratified = FALSE,
                            cc = FALSE,
                            skew = TRUE)

    mle <- res$estimates[, "est"]
    crude <- x1 / n1 - x2 / n2
    myindex <- (1:length(x1))[abs(mle - crude) > 0.02] # & abs(mle) < 0.0001]
    res$estimates[myindex, ]

  }


})
