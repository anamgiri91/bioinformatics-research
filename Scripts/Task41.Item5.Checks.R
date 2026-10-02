# Independent mathematical checks, including exhaustive predictor tie orders.
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- dirname(dirname(normalizePath(sub("^--file=", "", file_arg[1L]))))
source(file.path(root, "Scripts", "Task41.RRBS.Item5.Feasibility.Sep28.2026.R"))

stopifnot(identical(as.integer(levels_beta(c(0, .099, .1, .9, 1))), c(1L, 1L, 2L, 10L, 10L)))
stopifnot(levels_beta(.1 - 1e-12) == 1, levels_beta(.1 + 1e-12) == 2,
          is.na(levels_beta(NA_real_)))

all_orders <- function(v) {
  if (length(v) == 1L) return(matrix(v, 1))
  do.call(rbind, lapply(seq_along(v), function(i) cbind(v[i], all_orders(v[-i]))))
}
xi_ordered <- function(y) {
  n <- length(y); r <- vapply(y, function(z) sum(y <= z), 0L)
  l <- vapply(y, function(z) sum(y >= z), 0L)
  den <- 2 * sum(l * (n - l))
  if (!den) return(NA_real_)
  1 - n * sum(abs(diff(r))) / den
}
cases <- list(list(c(1, 1, 2, 2, 2), c(4, 1, 1, 3, 2)),
              list(rep(1, 5), c(1, 2, 2, 4, 5)),
              list(c(1, 2, 3, 4, 5), c(2, 1, 1, 3, 2)))
for (z in cases) {
  orders <- all_orders(seq_along(z[[1]]))
  orders <- orders[apply(orders, 1, function(i) !is.unsorted(z[[1]][i])), , drop = FALSE]
  brute <- mean(apply(orders, 1, function(i) xi_ordered(z[[2]][i])))
  stopifnot(abs(brute - xi_tie_average(z[[1]], z[[2]])) < 1e-12)
}
stopifnot(is.na(xi_tie_average(1:5, rep(1, 5))),
          abs(xi_tie_average(rep(1, 5), 1:5)) < 1e-12,
          abs(xi_tie_average(1:18, 1:18) - 16/19) < 1e-12)

a <- rbind(rep(.05, 18), rep(.05, 18), seq(.01, .18, length.out = 18),
           seq(.02, .98, length.out = 18))
b <- rbind(rep(.05, 18), rep(.95, 18), seq(.71, .88, length.out = 18),
           seq(.98, .02, length.out = 18))
z <- pair_metrics(a, b)
stopifnot(z$agreement_exact[1] == 1, z$ordinal_similarity[1] == 1,
          is.na(z$pearson_raw[1]), is.na(z$kappa_quadratic[1]),
          z$kappa_quadratic[2] == 0, z$ordinal_similarity[2] == 0,
          abs(z$pearson_raw[3] - 1) < 1e-12, abs(z$beta_MAE[3] - .7) < 1e-12,
          abs(z$pearson_raw[4] + 1) < 1e-12)
# Missingness must use the identical set of patients for both margins.
a <- matrix(c(.05, .2, NA, .7, .8), 1); b <- matrix(c(.1, NA, .9, .6, .7), 1)
z <- pair_metrics(a, b)
stopifnot(z$n_common == 3, abs(z$beta_MAE - mean(c(.05, .1, .1))) < 1e-12,
          abs(z$pearson_raw - cor(c(.05, .7, .8), c(.1, .6, .7))) < 1e-12)
cat("PASS: bin boundaries, exhaustive xi tie expectations, constants, identity limit,\n",
    "agreement versus dependence, paired missingness and exact MSD decomposition.\n")
