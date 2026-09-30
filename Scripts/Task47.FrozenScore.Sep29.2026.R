# ================================================================
# Task 47 - The frozen similarity score F+ for neighbouring CpGs
# Date: Sep 29, 2026
#
# FROZEN on 2026-09-29, before any new data were scored. The test
# plan that uses it is Results/Task47_FrozenProtocol.md. Do not edit
# this file: any change is a deviation and must be reported as one.
#
# THE SCORE (from Task 46: Codex's F plus the sign guard)
#   F+ = S * [w + (1 - w) * E+]
#   S   agreement: 1 - mean |x_k - y_k| over the reference samples,
#       on raw beta
#   E   co-methylation: squared distance correlation q, minus its
#       exact permutation mean q0, divided by (1 - q0), floored at 0.
#       q0 = tr(A) tr(B) / ((n - 1) sqrt(sum A^2 sum B^2)) for the
#       double-centred distance matrices A and B. E = 0 when a profile
#       is constant or 1 - q0 <= 1e-12.
#   E+  the sign guard: E when Spearman(x, y) > 0, else 0. Distance
#       correlation cannot see the sign, so an inverse pair would
#       otherwise get full credit. In Task 46 the guard changed the
#       test curves by at most 0.00001.
#   w   distance weight: exp(-gap / 200 bp). 200 bp is the decay
#       length of neighbour Spearman fitted in Task 46 (197 bp on
#       all 18 RRBS patients, 160 to 213 across folds).
#
# The reliability flags are reported beside F+, never multiplied in
# (Task 46: multiplying by mappability did not help). They are
# defined per data type in the protocol.
#
# USE
#   source("Scripts/Task47.FrozenScore.Sep29.2026.R")
#   score_F_plus(X, Y, gap)   # X, Y: pairs x samples, no NA
#   Running this file directly runs the self-test.
# ================================================================

suppressPackageStartupMessages(library(matrixStats))

FROZEN <- list(ell = 200, sign_guard = TRUE, chunk = 3000)

# E for many pairs at once (rows = pairs). Exact: matches the
# one-pair version of Task 46 to about 1e-15.
dcor_excess_vec <- function(X, Y) {
  n <- ncol(X); ia <- rep(seq_len(n), times = n); ib <- rep(seq_len(n), each = n)
  G <- matrix(0, n * n, n); G[cbind(seq_len(n * n), ia)] <- 1 / n
  dg <- (seq_len(n) - 1) * n + seq_len(n)
  cen <- function(Z) {
    D <- abs(Z[, ia, drop = FALSE] - Z[, ib, drop = FALSE])
    rm <- D %*% G
    D - rm[, ia, drop = FALSE] - rm[, ib, drop = FALSE] + rowMeans(rm)
  }
  A <- cen(X); B <- cen(Y)
  den <- sqrt(rowSums(A * A) * rowSums(B * B))
  q <- rowSums(A * B) / den
  q0 <- rowSums(A[, dg, drop = FALSE]) * rowSums(B[, dg, drop = FALSE]) / ((n - 1) * den)
  e <- (q - q0) / (1 - q0)
  e[!is.finite(e) | !(den > 0) | 1 - q0 <= 1e-12] <- 0
  pmax(e, 0)
}

spearman_rows <- function(X, Y) {
  a <- rowRanks(X, ties.method = "average"); b <- rowRanks(Y, ties.method = "average")
  a <- a - rowMeans(a); b <- b - rowMeans(b)
  den <- sqrt(rowSums(a^2) * rowSums(b^2))
  ifelse(den > 0, rowSums(a * b) / den, NA_real_)
}

# The frozen score. Returns the parts too, so tests can rebuild
# variants of the formula without recomputing E.
score_F_plus <- function(X, Y, gap, ell = FROZEN$ell, sign_guard = FROZEN$sign_guard,
                         chunk = FROZEN$chunk) {
  stopifnot(identical(dim(X), dim(Y)), nrow(X) == length(gap), !anyNA(X), !anyNA(Y), all(gap >= 0))
  S <- 1 - rowMeans(abs(X - Y))
  E <- numeric(nrow(X))
  for (s in seq(1, nrow(X), by = chunk)) {
    ix <- s:min(nrow(X), s + chunk - 1)
    E[ix] <- dcor_excess_vec(X[ix, , drop = FALSE], Y[ix, , drop = FALSE])
  }
  rho <- spearman_rows(X, Y)
  Ep <- if (sign_guard) E * (!is.na(rho) & rho > 0) else E
  w <- exp(-gap / ell)
  list(S = S, E = E, E_plus = Ep, spearman = rho, w = w, F_plus = S * (w + (1 - w) * Ep))
}

# ----------------------------------------------------------------
# self-test: the Task 46 constructed cases, 50 bp apart
# ----------------------------------------------------------------
if (sys.nframe() == 0L) {
  base <- seq(0.1, 0.9, length.out = 18)
  cases <- list(identical_constant = list(rep(0.05, 18), rep(0.05, 18)),
                different_constants = list(rep(0.05, 18), rep(0.95, 18)),
                identical_varying = list(base, base),
                offset_same_order = list(base * 0.5, base * 0.5 + 0.3),
                inverse = list(base, 1 - base),
                one_shared_rare = list(c(0.95, rep(0.05, 17)), c(0.95, rep(0.05, 17))),
                tiny_shared_wiggle = list(seq(0.49, 0.51, length.out = 18), seq(0.49, 0.51, length.out = 18)))
  X <- do.call(rbind, lapply(cases, `[[`, 1)); Y <- do.call(rbind, lapply(cases, `[[`, 2))
  # Codex F as run in Task 46 (100 bp, no sign guard): Results/Task46_Examples.csv
  old <- score_F_plus(X, Y, rep(50, 7), ell = 100, sign_guard = FALSE)$F_plus
  stopifnot(isTRUE(all.equal(unname(round(old, 4)), c(0.6065, 0.0607, 1, 0.7, 0.5765, 1, 1))))
  new <- score_F_plus(X, Y, rep(50, 7))
  out <- data.frame(case = names(cases), S = round(new$S, 4), E = round(new$E, 4),
                    spearman = round(new$spearman, 4), F_plus = round(new$F_plus, 4))
  print(out, row.names = FALSE)
  stopifnot(out$F_plus[out$case == "inverse"] < out$S[out$case == "inverse"])   # the guard acts
  cat("self-test passed\n")
}
