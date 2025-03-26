# Tests for consistency of confidence intervals vs p-values
# context("p-values")

n1 <- 40
n2 <- 20
xs <- expand.grid(0:n1, 0:n2)
x1 <- xs[, 1]
x2 <- xs[, 2]

for (skew in c(TRUE, FALSE)) {
  test_that("p-values consistent with confidence interval", {
    expect_equal(
      scoreci(
        x1 = x1, n1 = n1, x2 = x2, n2 = n2, skew = skew,
        contrast = "RD", precis = 10
      )$estimates[, "lower"] > 0 |
        scoreci(
          x1 = x1, n1 = n1, x2 = x2, n2 = n2, skew = skew,
          contrast = "RD", precis = 10
        )$estimates[, "upper"] < 0,
      scoreci(
        x1 = x1, n1 = n1, x2 = x2, n2 = n2, skew = skew,
        contrast = "RD", precis = 10
      )$pval[, "pval2sided"] < 0.05
    )
    expect_equal(
      scoreci(
        x1 = x1, n1 = n1, x2 = x2, n2 = n2, skew = skew,
        contrast = "RD", precis = 10
      )$estimates[, "lower"] > 0,
      scoreci(
        x1 = x1, n1 = n1, x2 = x2, n2 = n2, skew = skew,
        contrast = "RD", precis = 10
      )$pval[, "pval_right"] < 0.025
    )
    expect_equal(
      scoreci(
        x1 = x1, n1 = n1, x2 = x2, n2 = n2, skew = skew,
        contrast = "RD", precis = 10
      )$estimates[, "upper"] < 0,
      scoreci(
        x1 = x1, n1 = n1, x2 = x2, n2 = n2, skew = skew,
        contrast = "RD", precis = 10
      )$pval[, "pval_left"] < 0.025
    )
    expect_equal(
      scoreci(
        x1 = x1, n1 = n1, x2 = x2, n2 = n2, skew = skew,
        contrast = "RR", distrib = "poi", precis = 10
      )$estimates[, "lower"] > 1 |
        scoreci(
          x1 = x1, n1 = n1, x2 = x2, n2 = n2, skew = skew,
          contrast = "RR", distrib = "poi", precis = 10
        )$estimates[, "upper"] < 1,
      scoreci(
        x1 = x1, n1 = n1, x2 = x2, n2 = n2, skew = skew,
        contrast = "RR", distrib = "poi", precis = 10
      )$pval[, "pval2sided"] < 0.05
    )
    expect_equal(
      scoreci(
        x1 = x1, n1 = n1, x2 = x2, n2 = n2, skew = skew,
        contrast = "RR", distrib = "poi", precis = 10, theta0 = 2
      )$estimates[, "lower"] > 2,
      scoreci(
        x1 = x1, n1 = n1, x2 = x2, n2 = n2, skew = skew,
        contrast = "RR", distrib = "poi", precis = 10, theta0 = 2
      )$pval[, "pval_right"] < 0.025
    )
    expect_equal(
      scoreci(
        x1 = x1, n1 = n1, x2 = x2, n2 = n2, skew = skew,
        contrast = "RR", distrib = "poi", precis = 10, theta0 = 2
      )$estimates[, "upper"] < 2,
      scoreci(
        x1 = x1, n1 = n1, x2 = x2, n2 = n2, skew = skew,
        contrast = "RR", distrib = "poi", precis = 10, theta0 = 2
      )$pval[, "pval_left"] < 0.025
    )
    expect_equal(
      scoreci(
        x1 = x1, n1 = n1, x2 = x2, n2 = n2, skew = skew,
        contrast = "OR", precis = 10
      )$estimates[, "lower"] > 1 |
        scoreci(
          x1 = x1, n1 = n1, x2 = x2, n2 = n2, skew = skew,
          contrast = "OR", precis = 10
        )$estimates[, "upper"] < 1,
      scoreci(
        x1 = x1, n1 = n1, x2 = x2, n2 = n2, skew = skew,
        contrast = "OR", precis = 10
      )$pval[, "pval2sided"] < 0.05
    )
    expect_equal(
      scoreci(
        x1 = x1, n1 = n1, x2 = x2, n2 = n2, skew = skew,
        contrast = "OR", precis = 10
      )$estimates[, "lower"] > 1,
      scoreci(
        x1 = x1, n1 = n1, x2 = x2, n2 = n2, skew = skew,
        contrast = "OR", precis = 10
      )$pval[, "pval_right"] < 0.025
    )
    expect_equal(
      scoreci(
        x1 = x1, n1 = n1, x2 = x2, n2 = n2, skew = skew,
        contrast = "OR", precis = 10
      )$estimates[, "upper"] < 1,
      scoreci(
        x1 = x1, n1 = n1, x2 = x2, n2 = n2, skew = skew,
        contrast = "OR", precis = 10
      )$pval[, "pval_left"] < 0.025
    )
  })
}


n <- 10
combos <- NULL
for (a in 0:n) {
  for (b in 0:(n - a)) {
    for (c in 0:(n - a - b)) {
      for (d in 0:(n - a - b - c)) {
        combos <- rbind(combos, c(a = a, b = b, c = c, d = d))
      }
    }
  }
}
combos <- combos[rowSums(combos) == n, ]

test_that("RD p-values consistent with paired confidence interval", {
  expect_equal(
    (sapply(
      1:dim(combos)[1],
      function(i) pairbinci(x = combos[i, ], contrast = "RD",
                            method = "Score")$estimates[1] > 0
    )),
    (sapply(
      1:dim(combos)[1],
      function(i) pairbinci(x = combos[i, ], contrast = "RD",
                            method = "Score")$pval[, "pval_right"] < 0.025
    ))
  )
})
test_that("RR p-values consistent with paired confidence interval", {
  expect_equal(
    (sapply(
      1:dim(combos)[1],
      function(i) pairbinci(x = combos[i, ], contrast = "RR",
                            method = "Score")$estimates[1] > 1
    )),
    (sapply(
      1:dim(combos)[1],
      function(i) pairbinci(x = combos[i, ], contrast = "RR",
                            method = "Score")$pval[, "pval_right"] < 0.025
    ))
  )
})
test_that("OR p-values consistent with paired confidence interval", {
  expect_equal(
    unname(sapply(
      1:dim(combos)[1],
      function(i) {
        pairbinci(
          x = combos[i, ], contrast = "OR",
          method = "SCASp"
        )$estimates[, "lower", drop = FALSE] > 1
      }
    )),
    unname(sapply(
      1:dim(combos)[1],
      function(i) {
        pairbinci(
          x = combos[i, ], contrast = "OR",
          method = "SCASp"
        )$pval[, "pval_right", drop = FALSE] < 0.025
      }
    ))
  )
})
test_that("paired confidence interval (with bcf=FALSE) consistent with McNemar test", {
  expect_equal(
    (sapply(
      1:dim(combos)[1],
      function(i) {
        pairbinci(x = combos[i, ], contrast = "RD", method = "Score",
                  bcf = FALSE)$estimates[1] > 0 |
          pairbinci(x = combos[i, ], contrast = "RD", method = "Score",
                    bcf = FALSE)$estimates[3] < 0
      }
    )),
    (sapply(
      1:dim(combos)[1],
      function(i) {
        mctest <- mcnemar.test(x = matrix(combos[i, ], nrow = 2),
                               correct = FALSE)$p.value
        ifelse(is.na(mctest), FALSE, mctest < 0.05)
      }
    ))
  )
})
