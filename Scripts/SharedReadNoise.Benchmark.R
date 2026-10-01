# Three-way physical-fragment benchmark (plan_shared_read_noise.md sections 9-10).
# Each fragment of a donor goes to A, B or C with probability 1/3; all its calls
# stay together. Estimates use A only. The reference uses B and C, which share no
# fragments, so it carries no shared-read noise:
#   R    = (s(b_i^B, b_j^C) + s(b_i^C, b_j^B)) / 2         (covariance across donors)
#   Qref = mean_k (b_ik^B - b_jk^B)(b_ik^C - b_jk^C)        (squared disagreement)
# Per pair: dC = (Chat_A - R)^2 - (s_A - R)^2, negative favours the correction.
# Requires SharedReadNoise.Core.R (srn_pat_counts).
suppressPackageStartupMessages(library(data.table))

# Split each .pat row's multiplicity n into A, B, C (multinomial, p = 1/3 each).
srb_split_pat <- function(p, seed) {
  stopifnot(all(c("idx", "pat", "n") %in% names(p)), all(p$n >= 1))
  set.seed(seed)
  nA <- rbinom(nrow(p), p$n, 1/3); nB <- rbinom(nrow(p), p$n - nA, 1/2); nC <- p$n - nA - nB
  part <- function(m) { k <- which(m > 0); data.table(idx = p$idx[k], pat = p$pat[k], n = m[k]) }
  list(A = part(nA), B = part(nB), C = part(nC))
}

# One donor, one chromosome: eligible donor-pair rows with the A, B, C quantities.
srb_donor_chr <- function(pk, offset, ncpg, positions, seed, max_gap = 200L) {
  sp <- srb_split_pat(pk, seed)
  pc <- lapply(sp, function(q) if (nrow(q)) srn_pat_counts(q, offset, ncpg, positions, max_gap = max_gap) else NULL)
  if (any(vapply(pc, function(x) is.null(x) || !nrow(x), TRUE))) return(NULL)
  m <- merge(pc$A[, .(idx, NiA = Ni, MiA = Mi, NjA = Nj, MjA = Mj, n00, n01, n10, n11)],
             pc$B[, .(idx, NiB = Ni, MiB = Mi, NjB = Nj, MjB = Mj)], by = "idx")
  m <- merge(m, pc$C[, .(idx, NiC = Ni, MiC = Mi, NjC = Nj, MjC = Mj)], by = "idx")
  m[, rA := n00 + n01 + n10 + n11]
  m <- m[NiA >= 2 & NjA >= 2 & NiB >= 2 & NjB >= 2 & NiC >= 2 & NjC >= 2 & (rA == 0 | rA >= 2)]
  if (!nrow(m)) return(NULL)
  m[, `:=`(xA = MiA / NiA, yA = MjA / NjA, biB = MiB / NiB, bjB = MjB / NjB, biC = MiC / NiC, bjC = MjC / NjC)]
  m[, `:=`(vxA = xA * (1 - xA) / (NiA - 1), vyA = yA * (1 - yA) / (NjA - 1),
           cA = fifelse(rA >= 2, (rA * n11 - (n10 + n11) * (n01 + n11)) / (NiA * NjA * pmax(rA - 1, 1)), 0))]
  m[, .(idx, xA, yA, vxA, vyA, cA, biB, bjB, biC, bjC, NiA, NjA, rA)]
}

SRB_SUMS <- c("sx", "sy", "sxy", "sc", "sd2", "svx", "svy", "siB", "sjC", "siBjC", "siC", "sjB", "siCjB",
              "sq", "sNiA", "sNjA", "srA")
srb_accumulator <- function(size) {
  z <- data.table(n = integer(size))
  for (v in SRB_SUMS) set(z, j = v, value = numeric(size))
  z
}
srb_update <- function(acc, m, rows) {
  stopifnot(length(rows) == nrow(m), !anyNA(rows), !anyDuplicated(rows))
  add <- function(v, val) set(acc, rows, v, acc[[v]][rows] + val)
  set(acc, rows, "n", acc$n[rows] + 1L)
  add("sx", m$xA); add("sy", m$yA); add("sxy", m$xA * m$yA); add("sc", m$cA); add("sd2", (m$xA - m$yA)^2)
  add("svx", m$vxA); add("svy", m$vyA)
  add("siB", m$biB); add("sjC", m$bjC); add("siBjC", m$biB * m$bjC)
  add("siC", m$biC); add("sjB", m$bjB); add("siCjB", m$biC * m$bjB)
  add("sq", (m$biB - m$bjB) * (m$biC - m$bjC))
  add("sNiA", m$NiA); add("sNjA", m$NjA); add("srA", m$rA)
  invisible(acc)
}
# Per-pair estimates and losses from the sums (pairs with at least min_donors donors).
srb_finish <- function(a, min_donors = 20L) {
  a <- a[n >= min_donors]; n <- a$n
  cv <- function(sxy, sx, sy) (sxy - sx * sy / n) / (n - 1)
  W <- cv(a$sxy, a$sx, a$sy); mc <- a$sc / n; U <- W - mc
  R <- (cv(a$siBjC, a$siB, a$sjC) + cv(a$siCjB, a$siC, a$sjB)) / 2
  Draw <- a$sd2 / n; Dcor <- Draw - a$svx / n - a$svy / n + 2 * mc; Q <- a$sq / n
  data.table(n, raw_cov_A = W, cor_cov_A = U, ref_cov = R, mean_noise_A = mc,
             loss_raw = (W - R)^2, loss_cor = (U - R)^2, dC = (U - R)^2 - (W - R)^2,
             raw_d2_A = Draw, cor_d2_A = Dcor, ref_d2 = Q, dD = (Dcor - Q)^2 - (Draw - Q)^2,
             depth_A = (a$sNiA + a$sNjA) / (2 * n), overlap_A = (a$srA / n) / pmin(a$sNiA / n, a$sNjA / n))
}
