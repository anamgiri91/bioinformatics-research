# Section 6 comparisons of plan_shared_read_noise.md, on the three-way split.
# Extends SharedReadNoise.Benchmark.R, which stays unchanged because it is
# hashed in the Task 58 lock. Requires SharedReadNoise.Core.R,
# SharedReadNoise.Benchmark.R and Task47.FrozenScore.Sep29.2026.R.
#
# Correlation recovery (all eligible pairs, from sums):
#   Pearson        s_A / sqrt(s2_i s2_j)
#   variance-only  s_A / sqrt(V_i V_j),    V = s2 - mean(vhat)
#   full           Chat_A / sqrt(V_i V_j)
#   reference      R / sqrt(Ri Rj), Ri = s(b_i^B, b_i^C): B and C share no reads,
#                  so Ri carries no read noise
#   An estimate is valid when its variances are positive and it lies in [-1, 1].
# Pair ranking (subsample pairs, per-donor A values): scores from A only, judged
# by the mean B/C reference squared disagreement Qref of the top pairs.
suppressPackageStartupMessages({ library(data.table); library(matrixStats) })

S6_EXTRA <- c("sx2", "sy2", "siBiC", "sjBjC")
s6_accumulator <- function(size) {
  z <- srb_accumulator(size)
  for (v in S6_EXTRA) set(z, j = v, value = numeric(size))
  z
}
s6_update <- function(acc, m, rows) {
  srb_update(acc, m, rows)
  add <- function(v, val) set(acc, rows, v, acc[[v]][rows] + val)
  add("sx2", m$xA^2); add("sy2", m$yA^2); add("siBiC", m$biB * m$biC); add("sjBjC", m$bjB * m$bjC)
  invisible(acc)
}
s6_finish <- function(a, min_donors = 20L) {
  a <- a[n >= min_donors]; f <- srb_finish(a, min_donors); n <- a$n
  cv <- function(sxy, sx, sy) (sxy - sx * sy / n) / (n - 1)
  vx <- cv(a$sx2, a$sx, a$sx); vy <- cv(a$sy2, a$sy, a$sy)
  Vx <- vx - a$svx / n; Vy <- vy - a$svy / n
  Ri <- cv(a$siBiC, a$siB, a$siC); Rj <- cv(a$sjBjC, a$sjB, a$sjC)
  ratio <- function(num, d1, d2) { r <- num / sqrt(d1 * d2); r[!(d1 > 0 & d2 > 0)] <- NA_real_; r }
  f[, `:=`(var_x_A = vx, var_y_A = vy, vhat_x_A = Vx, vhat_y_A = Vy, ref_var_i = Ri, ref_var_j = Rj,
           r_pearson = ratio(raw_cov_A, vx, vy), r_vo = ratio(raw_cov_A, Vx, Vy),
           r_full = ratio(cor_cov_A, Vx, Vy), rho_ref = ratio(ref_cov, Ri, Rj))]
  for (v in c("r_vo", "r_full", "rho_ref")) set(f, which(abs(f[[v]]) > 1), v, NA_real_)
  f[]
}

# Rows of a pairs x donors matrix X, keeping only eligible donors (E == 1),
# packed into one matrix per eligible-donor count.
s6_pack <- function(E, X, rows) {
  ix <- which(E[rows, , drop = FALSE] == 1, arr.ind = TRUE); ix <- ix[order(ix[, 1], ix[, 2]), , drop = FALSE]
  n <- rowSums(E[rows, , drop = FALSE]); stopifnot(length(unique(n)) == 1L)
  matrix(X[cbind(rows[ix[, 1]], ix[, 2])], ncol = n[1], byrow = TRUE)
}

# Chatterjee's xi, rows = pairs: x -> y, ties in x broken at random (seeded).
s6_xi_rows <- function(X, Y, seed) {
  n <- ncol(X); set.seed(seed)
  o <- t(apply(X + matrix(runif(length(X)), nrow(X)) * 1e-9, 1, order))
  Ys <- matrix(Y[cbind(rep(seq_len(nrow(Y)), n), as.vector(o))], nrow(Y))
  r <- rowRanks(Ys, ties.method = "max"); l <- n + 1 - rowRanks(Ys, ties.method = "min")
  num <- n * rowSums(abs(r[, -1, drop = FALSE] - r[, -n, drop = FALSE])); den <- 2 * rowSums(l * (n - l))
  ifelse(den > 0, 1 - num / den, NA_real_)
}

# Per pair: mean and standard error of the per-donor corrected disagreement
# d_k = (x_k - y_k)^2 - vx_k - vy_k + 2 c_k, over eligible donors.
s6_d2_se <- function(subm) {
  E <- subm$E; n <- rowSums(E)
  d <- E * ((subm$xA - subm$yA)^2 - subm$vxA - subm$vyA + 2 * subm$cA)
  m <- rowSums(d) / n; v <- (rowSums(d^2) - n * m^2) / (n - 1)
  data.table(d2_mean = m, d2_se = sqrt(pmax(v, 0) / n))
}

# Ranking scores for subsample pairs. Higher = more similar. NA = undefined.
# The three noise-aware variants of the corrected disagreement were added in
# Task 63 development, after the plain one ranked noisy pairs first:
#   pos     max(Dhat2, 0)                 (an estimate of a square cannot be < 0)
#   ub      Dhat2 + 1.645 SE              (one-sided 95% upper bound)
#   pos_ub  max(Dhat2, 0) + 1.645 SE
s6_scores <- function(subm, gap, raw_d2, cor_d2, seed) {
  E <- subm$E; N <- rowSums(E); out <- data.table(S = NA_real_, spearman = NA_real_, dcor_E = NA_real_,
                                                  xi = NA_real_, F_plus = NA_real_)[rep(1L, nrow(E))]
  for (k in sort(unique(N[N >= 2]))) {
    rows <- which(N == k); X <- s6_pack(E, subm$xA, rows); Y <- s6_pack(E, subm$yA, rows)
    sc <- score_F_plus(X, Y, gap[rows])
    constant <- rowSds(X) == 0 | rowSds(Y) == 0
    xi <- pmax(s6_xi_rows(X, Y, seed + k), s6_xi_rows(Y, X, seed + 1000L + k)); xi[constant] <- NA_real_
    set(out, rows, "S", sc$S); set(out, rows, "spearman", sc$spearman); set(out, rows, "dcor_E", sc$E)
    set(out, rows, "xi", xi); set(out, rows, "F_plus", sc$F_plus)
  }
  se <- s6_d2_se(subm)$d2_se
  out[, `:=`(neg_raw_d2 = -raw_d2, neg_cor_d2 = -cor_d2, neg_pos_cor_d2 = -pmax(cor_d2, 0),
             neg_ub_cor_d2 = -(cor_d2 + 1.645 * se), neg_pos_ub_cor_d2 = -(pmax(cor_d2, 0) + 1.645 * se))][]
}

# Mean outcome y over the top fraction of pairs by score (NA ranked last),
# with fractional weights for pairs tied at the cut.
s6_top_mean <- function(score, y, frac) {
  s <- ifelse(is.na(score), -Inf, score); K <- frac * length(s)
  v <- sort(s, decreasing = TRUE)[ceiling(K)]
  above <- s > v; tie <- s == v
  w <- numeric(length(s)); w[above] <- 1; w[tie] <- (K - sum(above)) / sum(tie)
  sum(w * y) / K
}
S6_FRACTIONS <- c(0.01, 0.05, 0.10, 0.25)
S6_SCORES <- c(neg_cor_d2 = "corrected squared disagreement (lower first)",
               neg_pos_cor_d2 = "corrected, positive part (lower first)",
               neg_ub_cor_d2 = "corrected, upper 95% bound (lower first)",
               neg_pos_ub_cor_d2 = "corrected, positive part + 1.645 SE (lower first)",
               neg_raw_d2 = "observed squared disagreement (lower first)",
               S = "agreement S = 1 - mean |x - y| (Manhattan)",
               F_plus = "frozen score F+ (Task 47)", spearman = "Spearman", dcor_E = "distance correlation excess E",
               xi = "Chatterjee's xi (symmetric)")
# set = "all pairs" or "pairs that vary" (both sites vary across donors in A):
# among constant pairs every score agrees, so the second set shows where the
# choice of score matters.
s6_rank_table <- function(scores, qref, varies) {
  one <- function(keep, set) rbindlist(lapply(names(S6_SCORES), function(v) rbindlist(lapply(S6_FRACTIONS, function(f)
    data.table(set = set, score = v, label = S6_SCORES[[v]], top_fraction = f, pairs = sum(keep),
               valid_share = mean(!is.na(scores[[v]][keep])),
               mean_qref_top = s6_top_mean(scores[[v]][keep], qref[keep], f), mean_qref_all = mean(qref[keep]))))))
  rbind(one(rep(TRUE, length(qref)), "all pairs"), one(varies, "pairs that vary"))
}

# Correlation recovery over all eligible pairs. Failures are kept and counted;
# errors are against the reference where it is defined, on the pairs where all
# three estimates are valid (common) and on each estimator's own valid pairs.
s6_corr_table <- function(est) {
  d <- est[!is.na(rho_ref)]; all3 <- d[!is.na(r_pearson) & !is.na(r_vo) & !is.na(r_full)]
  rbindlist(lapply(c(r_pearson = "Pearson", r_vo = "variance-only", r_full = "full correction"), function(l) {
    v <- names(which(c(r_pearson = "Pearson", r_vo = "variance-only", r_full = "full correction") == l))
    own <- d[!is.na(d[[v]])]
    data.table(estimator = l, pairs = nrow(est), valid_share = mean(!is.na(est[[v]])),
               reference_defined_share = mean(!is.na(est$rho_ref)), common_pairs = nrow(all3),
               mse_common = mean((all3[[v]] - all3$rho_ref)^2), mean_error_common = mean(all3[[v]] - all3$rho_ref),
               own_pairs = nrow(own), mse_own = mean((own[[v]] - own$rho_ref)^2))
  }))
}
