# Method-development response to item 5: define, implement and stress-test scores.
# Proposed use is a two-score descriptor, not a claim of a novel correlation.
# C is a single-scale, attainable-range normalization of linearly weighted kappa.
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- dirname(dirname(normalizePath(sub("^--file=", "", file_arg[1L]))))
source(file.path(root, "Scripts", "Task41.RRBS.Item5.Feasibility.Sep28.2026.R"))

ordinal_candidate <- function(x, y) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  n <- length(x)
  if (!n) stop("No shared observations")
  d <- mean(abs(x - y))
  d0 <- mean(abs(outer(x, y, "-")))
  dmin <- mean(abs(sort(x) - sort(y)))
  dmax <- mean(abs(sort(x) - sort(y, decreasing = TRUE)))
  scale <- max(d0 - dmin, dmax - d0)
  defined <- scale > 1e-12
  c(D = d, D0 = d0, Dmin = dmin, Dmax = dmax, alignment_scale = scale,
    absolute_similarity = 1 - d / 9,
    linear_kappa = if (d0 > 1e-12) 1 - d / d0 else NA_real_,
    centered_alignment = if (defined) (d0 - d) / scale else NA_real_,
    piecewise_alignment_rejected = if (defined) {
      if (d <= d0) (d0 - d) / (d0 - dmin) else (d0 - d) / (dmax - d0)
    } else NA_real_)
}

# Independent checks against exhaustive permutations and edge cases.
permute_all <- function(v) {
  if (length(v) == 1L) return(matrix(v, 1))
  do.call(rbind, lapply(seq_along(v), function(i) cbind(v[i], permute_all(v[-i]))))
}
xcheck <- c(1, 1, 2, 3, 8); ycheck <- c(1, 2, 2, 3, 9)
all_y <- permute_all(ycheck)
all_scores <- t(apply(all_y, 1, function(y) ordinal_candidate(xcheck, y)))
stopifnot(abs(mean(all_scores[, "centered_alignment"])) < 1e-12,
          max(abs(all_scores[, "centered_alignment"])) <= 1 + 1e-12,
          abs(min(all_scores[, "D"]) - all_scores[1, "Dmin"]) < 1e-12,
          abs(max(all_scores[, "D"]) - all_scores[1, "Dmax"]) < 1e-12)
# Rare binary margins have exactly 18 possible locations of the rare sample.
rare_x <- c(10, rep(1, 17))
rare_scores <- t(vapply(1:18, function(j) {
  y <- rep(1, 18); y[j] <- 10; ordinal_candidate(rare_x, y)
}, numeric(9)))
stopifnot(abs(mean(rare_scores[, "centered_alignment"])) < 1e-12,
          abs(mean(rare_scores[, "piecewise_alignment_rejected"]) + 8/9) < 1e-12,
          abs(rare_scores[2, "centered_alignment"] + 1/17) < 1e-12,
          is.na(ordinal_candidate(rep(1:3, 6), rep(8:10, 6))["centered_alignment"]))
cat("PASS: exhaustive null mean and bounds; rare-event counterexample; disjoint-support degeneracy\n")

d <- fread(file.path(root, "Data", "NN.hg38.18P.forw.chr22.w.header.txt"))
m <- as.matrix(d[, 3:20]); obs <- is.finite(m)
idx <- which(rowSums(obs[-nrow(m), ] & obs[-1, ]) == 18)
x <- levels_beta(m[idx, , drop = FALSE]); y <- levels_beta(m[idx + 1L, , drop = FALSE])
score <- as.data.table(t(vapply(seq_along(idx), function(i) ordinal_candidate(x[i, ], y[i, ]), numeric(9))))
score <- cbind(data.table(row1 = idx, pos1 = d$pos[idx], pos2 = d$pos[idx + 1L]), score)
score[, defined := is.finite(centered_alignment)]
stopifnot(max(abs(score$centered_alignment), na.rm = TRUE) <= 1 + 1e-12,
          max(abs(score[(defined), centered_alignment - linear_kappa * D0 / alignment_scale])) < 1e-12)
fwrite(score, file.path(root, "Results", "Task42_CandidatePairs.csv.gz"))

summary <- data.table(pairs = nrow(score), defined = sum(score$defined),
  defined_share = mean(score$defined), median_absolute_similarity = median(score$absolute_similarity),
  median_centered_alignment = median(score$centered_alignment, na.rm = TRUE),
  median_linear_kappa = median(score$linear_kappa, na.rm = TRUE),
  share_positive_given_defined = mean(score$centered_alignment[score$defined] > 1e-12),
  share_at_plus_one_given_defined = mean(score$centered_alignment[score$defined] >= 1 - 1e-12),
  share_at_minus_one_given_defined = mean(score$centered_alignment[score$defined] <= -1 + 1e-12),
  max_weighted_kappa_equivalence_error = max(abs(score[(defined),
    centered_alignment - linear_kappa * D0 / alignment_scale])))
fwrite(summary, file.path(root, "Results", "Task42_CandidateSummary.csv"))

# Same random primary-pair selection as Task41, no tuning by candidate scores.
set.seed(20260928L)
sel <- sort(sample.int(length(idx), min(512L, length(idx))))
B <- 199L
null <- matrix(NA_real_, B, 3, dimnames = list(NULL, c("centered", "piecewise", "linear_kappa")))
exceed_d <- exceed_c <- integer(length(sel))
for (b in seq_len(B)) {
  order <- sample.int(18L)
  z <- t(vapply(sel, function(i) ordinal_candidate(x[i, ], y[i, order]), numeric(9)))
  null[b, ] <- colMeans(z[, c("centered_alignment", "piecewise_alignment_rejected", "linear_kappa")], na.rm = TRUE)
  exceed_d <- exceed_d + as.integer(z[, "D"] <= score$D[sel] + 1e-12)
  ok <- is.finite(z[, "centered_alignment"])
  exceed_c[ok] <- exceed_c[ok] + as.integer(z[ok, "centered_alignment"] >=
                                          score$centered_alignment[sel][ok] - 1e-12)
}
defined_sel <- score$defined[sel]
stopifnot(all(exceed_d[defined_sel] == exceed_c[defined_sel]))
diagnostic <- data.table(measure = colnames(null), null_mean = colMeans(null),
  null_mean_q025 = apply(null, 2, quantile, .025),
  null_mean_q975 = apply(null, 2, quantile, .975),
  sampled_pairs = length(sel), permutations = B,
  candidate_vs_L1_permutation_p_mismatches = sum(exceed_d[defined_sel] != exceed_c[defined_sel]))
fwrite(diagnostic, file.path(root, "Results", "Task42_CandidateNullDiagnostic.csv"))

# Representative constructed cases test semantics, not real-data biological truth.
examples <- list(
  identical_constant = list(rep(1, 18), rep(1, 18)),
  different_constants = list(rep(1, 18), rep(10, 18)),
  identical_varying = list(rep(1:9, 2), rep(1:9, 2)),
  inverse_varying = list(rep(1:9, 2), rep(9:1, 2)),
  shifted_perfect_correlation = list(rep(1:3, 6), rep(8:10, 6)),
  shared_single_rare_patient = list(rare_x, rare_x),
  separate_rare_patients = list(rare_x, c(1, 10, rep(1, 16))))
examples_table <- rbindlist(lapply(names(examples), function(name) {
  z <- examples[[name]]
  cbind(data.table(example = name), as.data.table(as.list(ordinal_candidate(z[[1]], z[[2]]))))
}))
fwrite(examples_table, file.path(root, "Results", "Task42_CandidateExamples.csv"))

# Quantify one-sample influence without selecting a favorable patient subset.
loo_changes <- matrix(NA_real_, length(sel), 18)
loo_undefined <- matrix(FALSE, length(sel), 18)
for (j in 1:18) {
  z <- vapply(sel, function(i) ordinal_candidate(x[i, -j], y[i, -j])["centered_alignment"], 0.0)
  loo_changes[, j] <- abs(z - score$centered_alignment[sel])
  loo_undefined[, j] <- !is.finite(z)
}
influence <- data.table(sampled_pairs = length(sel), full_defined = sum(defined_sel),
  full_defined_becomes_undefined_when_one_removed = sum(defined_sel & rowSums(loo_undefined) > 0),
  median_absolute_change = median(loo_changes, na.rm = TRUE),
  q95_absolute_change = unname(quantile(loo_changes, .95, na.rm = TRUE)),
  max_absolute_change = max(loo_changes, na.rm = TRUE))
fwrite(influence, file.path(root, "Results", "Task42_CandidateInfluence.csv"))

jsonlite::write_json(list(input_md5 = unname(tools::md5sum(file.path(root, "Data", "NN.hg38.18P.forw.chr22.w.header.txt"))),
  seed = 20260928L, complete_consecutive_pairs = length(idx), bins = 10L,
  permutation_pairs = length(sel), permutations = B,
  definition = "C=(D0-D)/max(D0-Dmin,Dmax-D0); undefined when denominator<=1e-12",
  novelty = "A proposed descriptor derived from weighted kappa; not established as a new correlation",
  script_md5 = unname(tools::md5sum(file.path(root, "Scripts", "Task42.OrdinalSimilarityCandidate.Sep28.2026.R"))),
  R_version = R.version.string), file.path(root, "Results", "Task42_RunManifest.json"), pretty = TRUE, auto_unbox = TRUE)
cat("PASS: candidate equals scaled weighted kappa; same per-pair permutation p-values as L1\n")
print(summary); print(diagnostic); print(influence)
