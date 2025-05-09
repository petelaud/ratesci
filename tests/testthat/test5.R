# Tests of iterative vs closed-form calculations
# context("Consistency")

# Exploring the extent of a problem with the non-iterative SCAS method when x=1 or 2
# (or n-1 or n-2) which had been apparent during development with level > 0.99.

# n <- round(runif(1000, 15, 9000), 0)
n <- round(seq(15, 9000, length.out = 1000))
xs <- rep_len(0:10, 1000)

rounded <- 10
for (level in c(0.9, 0.95, 0.99, 0.999)) {
  test_that("noniterative scas matches iterative version", {
    expect_equal(
      unname(scoreci(x1 = xs, n1 = n, contrast = "p", bcf = FALSE, level = level,
               precis = rounded + 0)$estimates[, c(1:3)]),
      unname(scaspci(x = xs, n = n, level = level)$estimates[, c(1:3)]) # Env test bug 9Nov2021 relates to this line
    )
    expect_equal(
      unname(scoreci(x1 = xs, n1 = n, contrast = "p", bcf = FALSE, cc = T,
               level = level, precis = rounded + 0)$estimates[, c(1:3)]),
      unname(scaspci(x = xs, n = n, cc = T, level = level)$estimates[, c(1:3)]) # Or is it this one? Fixed by removing Rmpfr code
    )
    expect_equal(
      unname(scoreci(x1 = xs, n1 = n, contrast = "p", bcf = FALSE, level = level,
               distrib = "poi", precis = rounded + 0)$estimates[, c(1:3)]),
      unname(scaspci(x = xs, n = n, level = level, distrib = "poi")$estimates[, c(1:3)])
    )
    expect_equal(
      unname(scoreci(x1 = xs, n1 = n, contrast = "p", bcf = FALSE, cc = T,
               level = level, distrib = "poi",
               precis = rounded + 0)$estimates[, c(1:3)]),
      unname(scaspci(x = xs, n = n, cc = T, level = level, distrib = "poi")$estimates[, c(1:3)])
    )
  })
}
