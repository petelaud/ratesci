library(contingencytables)
library(exact2x2)

?chap4

options(digits = 3)

# Attempt to reproduce results in Chapter 4, Table 4.12

x1 <- 0
n1 <- 5
x2 <- 0
n2 <- 5

x1 <- 6
n1 <- 10
x2 <- 6
n2 <- 20


x1 <- 0
n1 <- 16
x2 <- 15
n2 <- 72

x1 <- 7
n1 <- 34
x2 <- 1
n2 <- 34
ct <- matrix(c(x1, n1 - x1, x2, n2 - x2), nrow = 2, byrow = TRUE)

x <- matrix(c(x11, x12, x21, x22), nrow = 2, byrow = TRUE)
T0 <- -Fisher_exact_test_2x2(ct, "hypergeometric")$P




# Suissa-Schuster method (statistic = "Pearson") falls over when x1 = x2 = 0
x11 <- x12 <- 5
x21 <- x22 <- 5
N <- x11 + x12 + x21 + x22
T0 <- (N * (x11 * x22 - x12 * x21)^2) / ((x11 + x12) * (x21 + x22) * (x11 + x21) * (x12 + x22))


Pearson_chi_squared_test_2x2(ct)

ct <- matrix(c(0, 5, 0, 5), nrow = 2, byrow = TRUE)
Exact_unconditional_test_2x2(ct, statistic = "Pearson")$Pvalue



ct[1, 1] + ct[1, 2]

MiettinenNurminen_asymptotic_score_CI_difference_2x2(ct)
scoreci(x1, n1, x2, n2, skew=F)

Fisher_exact_test_2x2(ct)$Pvalue
Fisher_midP_test_2x2(ct)$Pvalue
Exact_unconditional_test_2x2(ct)$Pvalue
Exact_unconditional_test_2x2(ct, statistic = "LR")$Pvalue
Exact_unconditional_test_2x2(ct, statistic = "unpooled")$Pvalue

# Fisher-Boschloo result doesn't match vs Table 4.13 or 4.12,
# though it does match for Table 4.15
Exact_unconditional_test_2x2(ct, statistic = "Fisher")$Pvalue
Exact_unconditional_test_2x2(ct, statistic = "Fisher", gamma = 0)$Pvalue



names(MiettinenNurminen_asymptotic_score_CI_difference_2x2(t(perondi_2004)))

Exact_unconditional_test_2x2(perondi_2004, statistic = "Fisher", gamma = 0)$Pvalue
Exact_unconditional_test_2x2(perondi_2004, statistic = "Fisher", gamma = 1e-04)$Pvalue

# Internally calls
-Fisher_exact_test_2x2(ct, "hypergeometric")$P
# instead of
-Fisher_exact_test_2x2(ct, "Pearson")$P
# This ought to be documented, since all the other packages use the Pearson statistic

# This function does match T4.12, but doesn't have the Berger-Boos option
#
exact2x2::boschloo(x1, n1, x2, n2, midp = FALSE)$p.value
# Internally calls
exact2x2(ct)$p.value
ct <- matrix(c(5,29,3,31), 2, byrow=TRUE)
ct <- matrix(c(6,28,2,32), 2, byrow=TRUE)
ct <- matrix(c(5,29,3,31), 2, byrow=TRUE)
exact2x2(ct)$p.value

Fisher_exact_test_2x2(ct, "hypergeometric")$P
# exact2x2 function
exact2x2::fisher.exact(ct)$p.value
# base R function
stats::fisher.test(ct)$p.value


Fisher_exact_test_2x2(ct, "hypergeometric")$P
stats::fisher.test(ct)$p.value
exact2x2::fisher.exact(ct)$p.value

