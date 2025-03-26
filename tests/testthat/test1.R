# Tests for equivariant confidence intervals
# context("Symmetry")

n1 <- 40
n2 <- 20
xs <- expand.grid(0:n1, 0:n2)
x1 <- xs[ , 1]
x2 <- xs[ , 2]

for (cc in c(FALSE, TRUE)) {
  for (skew in c(TRUE, FALSE)) {
    test_that("Transposed inputs produce inverted intervals", {
      expect_equal(
        unname(scoreci(
          x1 = x1, n1 = n1, x2 = x2, n2 = n2, skew = skew,
          cc = cc, contrast = "RD", precis = 10
        )$estimates[, 1]),
        -unname(scoreci(
          x1 = x2, n1 = n2, x2 = x1, n2 = n1, skew = skew,
          cc = cc, contrast = "RD", precis = 10
        )$estimates[, 3])
      )
      expect_equal(
        unname(scoreci(
          x1 = x1, n1 = n1, x2 = x2, n2 = n2, skew = skew,
          cc = cc, contrast = "RD", distrib = "poi", precis = 10
        )$estimates[, 1]),
        -unname(scoreci(
          x1 = x2, n1 = n2, x2 = x1, n2 = n1, skew = skew,
          cc = cc, contrast = "RD", distrib = "poi", precis = 10
        )$estimates[, 3])
      )
      expect_equal(
        (unname(scoreci(
          x1 = x1, n1 = n1, x2 = x2, n2 = n2, skew = skew, cc = cc,
          contrast = "RR", rr_tang = F, precis = 10
        )$estimates)[, 1]),
        signif(1 / unname(scoreci(
          x1 = x2, n1 = n2, x2 = x1, n2 = n1, skew = skew, cc = cc,
          contrast = "RR", rr_tang = F, precis = 10
        )$estimates)[, 3]),
        tolerance = 1E-5
      )
      expect_equal(
        (unname(scoreci(
          x1 = x1, n1 = n1, x2 = x2, n2 = n2, skew = skew, cc = cc,
          contrast = "RR", rr_tang = T, precis = 10
        )$estimates)[, 1]),
        (1 / unname(scoreci(
          x1 = x2, n1 = n2, x2 = x1, n2 = n1, skew = skew, cc = cc,
          contrast = "RR", rr_tang = T, precis = 10
        )$estimates)[, 3]),
        tolerance = 1E-5
      )
      expect_equal(
        (unname(scoreci(
          x1 = x1, n1 = n1, x2 = x2, n2 = n2, skew = skew, cc = cc,
          contrast = "RR", distrib = "poi", precis = 10
        )$estimates)[, 1]),
        (1 / unname(scoreci(
          x1 = x2, n1 = n2, x2 = x1, n2 = n1, skew = skew, cc = cc,
          contrast = "RR", distrib = "poi", precis = 10
        )$estimates)[, 3]),
        tolerance = 1E-5
      )
      expect_equal(
        (unname(scoreci(
          x1 = x1, n1 = n1, x2 = x2, n2 = n2, skew = skew,
          cc = cc, contrast = "OR", precis = 10
        )$estimates)[, 1]),
        (unname(1 / scoreci(
          x1 = x2, n1 = n2, x2 = x1, n2 = n1, skew = skew,
          cc = cc, contrast = "OR", precis = 10
        )$estimates)[, 3]),
        tolerance = 1E-5
      )
    })
  }
}


for (midp in c(FALSE, TRUE)) {
    test_that("Transposed inputs produce inverted intervals", {
      expect_equal(
        unname(exactci(x = 0:40, n = 40,
                       midp = midp, distrib = "bin"
                       )[, 1]),
        1 - unname(exactci(x = 40:0, n = 40,
                           midp = midp, distrib = "bin"
                           )[, 3])
        )
    })
}
for (distrib in c("bin", "poi")) {
  for (midp in c(FALSE, TRUE)) {
  test_that("Midpoint equals UCL and LCL for level = 0", {
    expect_equal(
      unname(exactci(x = 0:40, n = 40,
                     midp = midp, distrib = "poi"
      )[, 2]),
      unname(exactci(x = 0:40, n = 40,
                         midp = TRUE, distrib = "poi", level = 0
      )[, 1]),
      tolerance = 1E-5
    )
    expect_equal(
      unname(exactci(x = 0:40, n = 40,
                     midp = midp, distrib = "poi"
      )[, 2]),
      unname(exactci(x = 0:40, n = 40,
                     midp = TRUE, distrib = "poi", level = 0
      )[, 3]),
      tolerance = 1E-5
    )
  })
}}


test_that("Tang score matches IVS", {
  expect_equal(
    signif(unname(scoreci(
      x1 = x1, n1 = n1, x2 = x2, n2 = n2, skew = skew,
      cc = cc, contrast = "RR", rr_tang = F, precis = 100
    )$estimates)[, 1], digits = 5),
    signif(unname(scoreci(
      x1 = x1, n1 = n1, x2 = x2, n2 = n2, skew = skew,
      cc = cc, contrast = "RR", rr_tang = T, precis = 100
    )$estimates)[, 1], digits = 5)
  )
})

if (FALSE) {
  # Check symmetry of results for 2-stratum examples
  # (Warning: long runtime for this test - not included in CRAN submission)
  n1 <- n2 <- 10
  xs <- expand.grid(0:n1, 0:n1, 0:n2, 0:n2)
  x1 <- xs[, c(1, 2)]
  x2 <- xs[, c(3, 4)]

  cc <- FALSE
  skew <- FALSE
  i <- 3
  # for (i in 1:dim(x1)[1]) {
  for (i in 1:1000) {
    test_that("Transposed stratified inputs produce inverted intervals", {
      expect_equal(
        unname(scoreci(
          x1 = unlist(x1[i, ]), n1 = c(n1, n1), x2 = unlist(x2[i, ]),
          n2 = c(n2, n2), stratified = T, skew = skew, cc = cc,
          contrast = "RD", precis = 10
        )$estimates[, 1]),
        -unname(scoreci(
          x1 = unlist(x2[i, ]), n1 = c(n2, n2), x2 = unlist(x1[i, ]),
          n2 = c(n1, n1), stratified = T, skew = skew, cc = cc,
          contrast = "RD", precis = 10
        )$estimates[, 3])
      )
    })
  }
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

# Exclude cases with zero marginal totals, to avoid comparing results that are
# = Inf using closed-form code, vs very large using iterative method
zero1 <- rowSums(combos[, c(1, 2)]) == 0
zero2 <- rowSums(combos[, c(1, 3)]) == 0
combos <- combos[(!zero1 & !zero2), ]

test_that("Transposed inputs produce inverted paired intervals", {
  expect_equal(
    t(unname(round(sapply(
      1:dim(combos)[1],
      function(i) pairbinci(x = combos[i, ],
                            contrast = "RD")$estimates[, c(1, 3)]
    ), 5))),
    t(unname(round(-sapply(
      1:dim(combos)[1],
      function(i) pairbinci(x = combos[i, c(1, 3, 2, 4)],
                            contrast = "RD")$estimates[, c(3, 1)]
    ), 5)))
  )
  expect_equal(
    t(unname(round(sapply(
      1:dim(combos)[1],
      function(i) pairbinci(x = combos[i, ],
                            contrast = "RD", cc = T)$estimates[, c(1, 3)]
    ), 5))),
    t(unname(round(-sapply(
      1:dim(combos)[1],
      function(i) pairbinci(x = combos[i, c(1, 3, 2, 4)],
                            contrast = "RD", cc = T)$estimates[, c(3, 1)]
    ), 5)))
  )
  expect_equal(
    t(unname(round(sapply(
      1:dim(combos)[1],
      function(i) pairbinci(x = combos[i, ], contrast = "RD",
                            method = "Score", skew = T)$estimates[, c(1, 3)]
    ), 5))),
    t(unname(round(-sapply(
      1:dim(combos)[1],
      function(i) pairbinci(x = combos[i, c(1, 3, 2, 4)], contrast = "RD",
                            method = "Score", skew = T)$estimates[, c(3, 1)]
    ), 5)))
  )
  expect_equal(
    t(unname(round(1 / round(1 / sapply(
      1:dim(combos)[1],
      function(i) pairbinci(x = combos[i, ], contrast = "RR",
                            precis = 8)$estimates[, c(1, 3)]
    ), 5), 5))),
    t(unname(round(1 / round(sapply(
      1:dim(combos)[1],
      function(i) pairbinci(x = combos[i, c(1, 3, 2, 4)], contrast = "RR",
                            precis = 8)$estimates[, c(3, 1)]
    ), 5), 5)))
  )
  expect_equal(
    t(unname(round(1 / round(1 / sapply(
      1:dim(combos)[1],
      function(i) pairbinci(x = combos[i, ], contrast = "RR", cc = T,
                            precis = 8)$estimates[, c(1, 3)]
    ), 5), 5))),
    t(unname(round(1 / round(sapply(
      1:dim(combos)[1],
      function(i) pairbinci(x = combos[i, c(1, 3, 2, 4)], contrast = "RR",
                            cc = T, precis = 8)$estimates[, c(3, 1)]
    ), 5), 5)))
  )
  expect_equal(
    t(unname(round(1 / round(1 / sapply(
      1:dim(combos)[1],
      function(i) pairbinci(x = combos[i, ], contrast = "RR", method = "Score",
                            skew = T, precis = 8)$estimates[, c(1, 3)]
    ), 5), 5))),
    t(unname(round(1 / round(sapply(
      1:dim(combos)[1],
      function(i) pairbinci(x = combos[i, c(1, 3, 2, 4)],
                            contrast = "RR", method = "Score", skew = T,
                            precis = 8)$estimates[, c(3, 1)]
    ), 5), 5)))
  )
})

test_that("Closed form methods match iterative results for paired intervals", {
  expect_equal(
    t(unname(
      round(sapply(
        1:dim(combos)[1],
        function(i) {
          pairbinci(x = combos[i, ], contrast = "RD", method = "Score",
                    skew = FALSE, precis = 8)$estimates[, c(1, 3)]
        }
      ), 5)
    )),
    t(unname(
      round(sapply(
        1:dim(combos)[1],
        function(i) {
          pairbinci(x = combos[i, ], contrast = "RD", method = "Score_closed",
                    skew = FALSE, precis = 8)$estimates[, c(1, 3)]
        }
      ), 5)
    ))
  )

  expect_equal(
    t(unname(
      round(sapply(
        1:dim(combos)[1],
        function(i) {
          pairbinci(x = combos[i, ], contrast = "RD", method = "Score",
                    skew = FALSE, bcf = FALSE, cc = TRUE,
                    precis = 8)$estimates[, c(1, 3)]
        }
      ), 5)
    )),
    t(unname(
      round(sapply(
        1:dim(combos)[1],
        function(i) {
          pairbinci(x = combos[i, ], contrast = "RD", method = "Score_closed",
                    skew = FALSE, bcf = FALSE, cc = TRUE,
                    precis = 8)$estimates[, c(1, 3)]
        }
      ), 5)
    ))
  )

  expect_equal(
    t(unname(
      round(sapply(
        1:dim(combos)[1],
        function(i) {
          pairbinci(x = combos[i, ], contrast = "RR", method = "Score",
                    skew = FALSE, precis = 8)$estimates[, c(1, 3)]
        }
      ), 3)
    )),
    t(unname(
      round(sapply(
        1:dim(combos)[1],
        function(i) {
          pairbinci(x = combos[i, ], contrast = "RR", method = "Score_closed",
                    skew = FALSE, precis = 8)$estimates[, c(1, 3)]
        }
      ), 3)
    ))
  )

  # These ones match to 2dps
  expect_equal(
    t(unname(
      round(sapply(
        1:dim(combos)[1],
        function(i) {
          pairbinci(x = combos[i, ], contrast = "RR", method = "Score",
                    skew = FALSE, cc = 0.5, precis = 8)$estimates[, c(1, 3)]
        }
      ), 2)
    )),
    t(unname(
      round(sapply(
        1:dim(combos)[1],
        function(i) {
#          for (i in 1:dim(combos)[1]) {
          pairbinci(x = combos[i, ], contrast = "RR", method = "Score_closed",
                    skew = FALSE, cc = 0.5)$estimates[, c(1, 3)]
        }
#          }
      ), 2)
    ))
  )
})


if (FALSE) {
  scoreci(x1 = 0, x2 = 0, n1 = 20, n2 = 10, contrast = "RR", plot = T)
  scoreci(x1 = 100, x2 = 0, n1 = 200, n2 = 5, plot = T)
  scoreci(x1 = 0, x2 = 100, n1 = 5, n2 = 200, plot = T)

  scoreci(x1 = 0, n1 = 80000, contrast = "p", plot = T)

  x1 <- 0
  n1 <- 800
  x2 <- NULL
  n2 <- NULL
  contrast <- "p"
  distrib <- "bin"
  bcf <- TRUE
  skew <- TRUE
  myfun <- function(theta, ...) {
    #  scoretheta(theta = theta, x1 = x1, x2 = x2, n1 = n1, n2 = n2, bcf = bcf,
    #             contrast = contrast, distrib = distrib, skew = skew)$score
    scoretheta(theta = theta, x1 = x1, n1 = n1, contrast = "p", skew = T)$score
  }
  myfun(theta = 0.00064)
  myfun(theta = 0.00065)

  myseq <- seq(0, 0.001, length.out = 400)
  dim1 <- 1
  sc <- array(sapply(myseq, function(x) myfun(x)), dim = c(dim1, length(myseq)))
  myfun(theta = 0.00064)
  myfun(theta = 0.00065)
  scoretheta(theta = 0.00064, x1 = 0, n1 = 800, distrib = "bin", contrast = "p", skew = T)$score
  scoretheta(theta = 0.00065, x1 = 0, n1 = 800, distrib = "bin", contrast = "p", skew = T)$score
}
