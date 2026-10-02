# Item 5: reproducible RRBS adjacent-CpG similarity analysis.
# Run from the repository: Rscript Scripts/Task41.RRBS.Item5.Feasibility.Sep28.2026.R
# Existing Task39 outputs and the separate TCGA outlier project are not modified.
# Dependencies: data.table, matrixStats, ggplot2, jsonlite (already installed).
# Agreement, dependence and sequencing uncertainty are different estimands.
suppressPackageStartupMessages({library(data.table); library(matrixStats)})

levels_beta <- function(x, bins = 10L) {
  if (any(is.finite(x) & (x < 0 | x > 1))) stop("Beta outside [0,1]")
  # findInterval implements left-closed bins; it does not move values just
  # below a boundary upward by rounding. Preserve missingness and dimensions.
  ans <- rep(NA_real_, length(x)); ok <- is.finite(x)
  ans[ok] <- pmin(findInterval(x[ok], (0:bins) / bins), bins)
  dim(ans) <- dim(x)
  ans
}

xi_tie_average <- function(x, y) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  n <- length(y)
  if (n < 2L || length(unique(y)) == 1L) return(NA_real_)
  r <- rank(y, ties.method = "max")
  l <- n - rank(y, ties.method = "min") + 1
  den <- 2 * sum(l * (n - l))
  # Exact expectation over all uniform random orders within tied X groups.
  # Within a group of m, each unordered pair is adjacent with probability 2/m.
  # Between groups, the boundary endpoints are uniform in their own groups.
  groups <- split(seq_len(n), factor(x, levels = sort(unique(x))))
  distance <- abs(outer(r, r, "-"))
  jump <- sum(vapply(groups, function(g) sum(distance[g, g, drop = FALSE]) /
                      length(g), 0.0))
  if (length(groups) > 1L) for (k in seq_len(length(groups) - 1L)) {
    jump <- jump + mean(distance[groups[[k]], groups[[k + 1L]], drop = FALSE])
  }
  1 - n * jump / den
}

row_correlation <- function(a, b) {
  constant <- rowMins(a, na.rm = TRUE) == rowMaxs(a, na.rm = TRUE) |
    rowMins(b, na.rm = TRUE) == rowMaxs(b, na.rm = TRUE)
  a <- a - rowMeans(a, na.rm = TRUE); b <- b - rowMeans(b, na.rm = TRUE)
  den <- sqrt(rowSums(a^2, na.rm = TRUE) * rowSums(b^2, na.rm = TRUE))
  ans <- rowSums(a * b, na.rm = TRUE) / den
  ans[constant | den == 0 | !is.finite(ans)] <- NA_real_
  pmax(-1, pmin(1, ans))
}

pair_metrics <- function(a, b, bins = 10L, include_xi = TRUE) {
  stopifnot(identical(dim(a), dim(b)))
  ok <- is.finite(a) & is.finite(b); a[!ok] <- NA; b[!ok] <- NA
  n <- rowSums(ok); x <- levels_beta(a, bins); y <- levels_beta(b, bins)
  am <- rowMeans(a, na.rm = TRUE); bm <- rowMeans(b, na.rm = TRUE)
  va <- rowMeans((a - am)^2, na.rm = TRUE)
  vb <- rowMeans((b - bm)^2, na.rm = TRUE)
  covab <- rowMeans((a - am) * (b - bm), na.rm = TRUE)
  ca <- rowMins(a, na.rm = TRUE) == rowMaxs(a, na.rm = TRUE)
  cb <- rowMins(b, na.rm = TRUE) == rowMaxs(b, na.rm = TRUE)
  va[ca] <- 0; vb[cb] <- 0; covab[ca | cb] <- 0
  ordinal_obs <- rowMeans((x - y)^2, na.rm = TRUE)
  ordinal_exp <- rowMeans(x^2, na.rm = TRUE) + rowMeans(y^2, na.rm = TRUE) -
    2 * rowMeans(x, na.rm = TRUE) * rowMeans(y, na.rm = TRUE)
  ccc_den <- va + vb + (am - bm)^2
  out <- data.table(n_common = n,
    agreement_exact = rowMeans(x == y, na.rm = TRUE),
    agreement_within_one = rowMeans(abs(x - y) <= 1, na.rm = TRUE),
    ordinal_L1 = rowMeans(abs(x - y), na.rm = TRUE),
    ordinal_similarity = 1 - rowMeans(abs(x - y), na.rm = TRUE) / (bins - 1),
    beta_MAE = rowMeans(abs(a - b), na.rm = TRUE),
    beta_RMSE = sqrt(rowMeans((a - b)^2, na.rm = TRUE)),
    kappa_quadratic = fifelse(ordinal_exp > 0, 1 - ordinal_obs / ordinal_exp, NA_real_),
    concordance_raw = fifelse(ccc_den > 0, 2 * covab / ccc_den, NA_real_),
    pearson_raw = row_correlation(a, b),
    spearman_raw = row_correlation(rowRanks(a, ties.method = "average", preserveShape = TRUE),
                                  rowRanks(b, ties.method = "average", preserveShape = TRUE)),
    spearman_levels = row_correlation(rowRanks(x, ties.method = "average", preserveShape = TRUE),
                                     rowRanks(y, ties.method = "average", preserveShape = TRUE)),
    constant_raw_any = va == 0 | vb == 0,
    constant_levels_any = rowSds(x, na.rm = TRUE) == 0 | rowSds(y, na.rm = TRUE) == 0,
    squared_level_difference = (am - bm)^2,
    squared_spread_difference = (sqrt(va) - sqrt(vb))^2,
    pattern_MSD = 2 * (sqrt(va * vb) - covab))
  # Constant-safe exact identity: RMSE^2 = level + spread + pattern.
  stopifnot(max(abs(out$beta_RMSE^2 - out$squared_level_difference -
                       out$squared_spread_difference - out$pattern_MSD), na.rm = TRUE) < 1e-10)
  if (include_xi) {
    out[, xi_xy := vapply(seq_len(.N), function(i) xi_tie_average(x[i, ], y[i, ]), 0.0)]
    out[, xi_yx := vapply(seq_len(.N), function(i) xi_tie_average(y[i, ], x[i, ]), 0.0)]
    out[, xi_mean := (xi_xy + xi_yx) / 2]
    out[, xi_max := pmax(xi_xy, xi_yx)]
  }
  out
}

metric_names <- c("agreement_exact", "agreement_within_one", "ordinal_similarity",
                  "beta_MAE", "beta_RMSE", "kappa_quadratic", "concordance_raw",
                  "pearson_raw", "spearman_raw", "spearman_levels", "xi_xy", "xi_yx",
                  "xi_mean", "xi_max")

metric_summary <- function(tab, grouping) {
  long <- melt(tab, id.vars = grouping, measure.vars = metric_names,
               variable.name = "measure", value.name = "value")
  long[, {
    z <- value[is.finite(value)]
    list(pairs = .N, defined = length(z), defined_share = length(z) / .N,
         mean = if (length(z)) mean(z) else NA_real_,
         median = if (length(z)) median(z) else NA_real_,
         q25 = if (length(z)) unname(quantile(z, .25)) else NA_real_,
         q75 = if (length(z)) unname(quantile(z, .75)) else NA_real_)
  }, by = c(grouping, "measure")]
}

run_item5 <- function(root) {
  started <- proc.time()[3]
  input <- file.path(root, "Data", "NN.hg38.18P.forw.chr22.w.header.txt")
  output <- file.path(root, "Results")
  write_result <- function(x, name) fwrite(x, file.path(output, paste0("Task41_", name)))
  d <- fread(input, na.strings = "NA")
  samples <- paste0("P", 1:18, "NF")
  stopifnot(identical(names(d), c("chr", "pos", samples)), !anyDuplicated(d[, .(chr, pos)]),
            uniqueN(d$chr) == 1L, all(diff(d$pos) > 0))
  m <- as.matrix(d[, ..samples]); observed <- is.finite(m)
  stopifnot(all(m[observed] >= 0 & m[observed] <= 1))
  n_site <- nrow(m); complete <- which(rowSums(observed) == 18L)
  common <- rowSums(observed[-n_site, ] & observed[-1, ])
  cat("Input:", n_site, "sites;", length(complete), "complete;",
      sum(common == 18), "true consecutive complete pairs\n")
  write_result(data.table(n_common = 0:18,
                         true_consecutive_pairs = tabulate(common + 1L, nbins = 19)),
               "Coverage.csv")
  input_audit <- list(input = input, input_md5 = unname(tools::md5sum(input)),
    sites = n_site, samples = samples, observed_cells = sum(observed),
    missing_cells = sum(!observed), missing_share = mean(!observed),
    complete_sites = length(complete), consecutive_file_pairs = length(common),
    consecutive_file_complete_pairs = sum(common == 18),
    consecutive_complete_pairs = length(complete) - 1,
    complete_bridges = sum(diff(complete) > 1),
    coordinate_basis = "Input hg38 coordinates preserved; 0/1-based origin not independently verified",
    provenance_limit = "No methylated counts, total coverage, covariates or patient metadata supplied")
  jsonlite::write_json(input_audit, file.path(output, "Task41_InputAudit.json"),
                       pretty = TRUE, auto_unbox = TRUE)

  # All pairs with >=6 jointly observed values form the sensitivity table.
  # This is exploratory, not a validated coverage threshold. Primary = all 18.
  idx <- which(common >= 6L)
  cat("Computing", length(idx), "true consecutive pairs with >=6 shared observations\n")
  meta <- data.table(row1 = idx, row2 = idx + 1L, pos1 = d$pos[idx], pos2 = d$pos[idx + 1L],
                     skipped_rows = 0L, pair_set = "true_consecutive")
  calc <- pair_metrics(m[idx, , drop = FALSE], m[idx + 1L, , drop = FALSE], include_xi = FALSE)
  # Xi below 12 joint observations is deliberately not interpreted or exported.
  # Computed scores at >=12 remain descriptive, not certified population estimates.
  xi_idx <- which(calc$n_common >= 12L)
  lx <- levels_beta(m[idx[xi_idx], , drop = FALSE]); ly <- levels_beta(m[idx[xi_idx] + 1L, , drop = FALSE])
  calc[, c("xi_xy", "xi_yx", "xi_mean", "xi_max") := NA_real_]
  xy <- vapply(seq_along(xi_idx), function(k) xi_tie_average(lx[k, ], ly[k, ]), 0.0)
  yx <- vapply(seq_along(xi_idx), function(k) xi_tie_average(ly[k, ], lx[k, ]), 0.0)
  calc[xi_idx, `:=`(xi_xy = xy, xi_yx = yx, xi_mean = (xy + yx) / 2, xi_max = pmax(xy, yx))]
  pairs <- cbind(meta, calc)

  # Reconstruct the complete-site bridges explicitly instead of calling them adjacent.
  bridge_idx <- which(diff(complete) > 1L)
  bi <- complete[bridge_idx]; bj <- complete[bridge_idx + 1L]
  bridges <- cbind(data.table(row1 = bi, row2 = bj, pos1 = d$pos[bi], pos2 = d$pos[bj],
                             skipped_rows = bj - bi - 1L, pair_set = "complete_site_bridge"),
                    pair_metrics(m[bi, , drop = FALSE], m[bj, , drop = FALSE]))
  pairs <- rbind(pairs, bridges)
  pairs[, gap_bp := pos2 - pos1]
  pairs[, gap_bin := cut(gap_bp, c(0, 10, 50, 100, 200, 500, 1000, 10000, Inf),
                         right = TRUE, dig.lab = 6)]
  setorder(pairs, pair_set, pos1)
  write_result(pairs, "PairMetrics.csv.gz")
  primary <- pairs[pair_set == "true_consecutive" & n_common == 18]
  complete_pairs <- pairs[n_common == 18]
  write_result(metric_summary(complete_pairs, c("pair_set", "gap_bin")), "MetricsByGap.csv")
  write_result(metric_summary(complete_pairs, "pair_set"), "MetricsOverall.csv")
  # States from Task 39 (the TSG-paper rule), so the report can show where each
  # score breaks down. Primary pairs only: both sites observed in all 18.
  state_path <- file.path(output, "Task39_RRBS_CompleteSiteStates.csv.gz")
  if (file.exists(state_path)) {
    st <- fread(state_path, select = c("pos", "state"))
    sp <- copy(primary)
    sp[, s1 := st$state[match(pos1, st$pos)]][, s2 := st$state[match(pos2, st$pos)]]
    stopifnot(!anyNA(sp$s1), !anyNA(sp$s2))
    sp[, state_pair := fifelse(s1 == s2, paste0(s1, "-", s2), "different")]
    write_result(metric_summary(sp, "state_pair"), "MetricsByStatePair.csv")
  }
  sensitivity <- rbindlist(lapply(c(6, 9, 12, 15, 18), function(nmin) {
    z <- copy(pairs[pair_set == "true_consecutive" & n_common >= nmin]); z[, min_common := nmin]
    metric_summary(z, "min_common")
  }))
  write_result(sensitivity, "CoverageSensitivity.csv")

  # Verify all deterministic Task39 measures against its saved rounded export.
  old_path <- file.path(output, "Task39_RRBS_PairSimilarity.csv.gz")
  if (file.exists(old_path)) {
    old <- fread(old_path)
    joined <- merge(old, complete_pairs, by = c("pos1", "pos2"))
    stopifnot(nrow(joined) == length(complete) - 1L)
    mapping <- c(agree = "agreement_exact", agree1 = "agreement_within_one", L1 = "ordinal_L1",
                 L2beta = "beta_RMSE", kappa_w = "kappa_quadratic", pearson_beta = "pearson_raw",
                 spearman = "spearman_levels")
    # Task 39 renamed pearson -> pearson_beta on 2026-09-28; a missing column
    # must stop the check, not compare nothing and pass.
    stopifnot(all(names(mapping) %in% names(joined)), all(mapping %in% names(joined)))
    check <- rbindlist(lapply(names(mapping), function(old_name) {
      new_name <- mapping[[old_name]]; a <- joined[[old_name]]; b <- joined[[new_name]]
      data.table(Task39 = old_name, Task41 = new_name, pairs = length(a),
        missing_disagreements = sum(is.na(a) != is.na(b)),
        max_absolute_difference = max(abs(a - b), na.rm = TRUE))
    }))
    write_result(check, "Task39_Reproduction.csv")
    stopifnot(all(check$missing_disagreements == 0), all(check$max_absolute_difference <= 5.01e-5))
  }

  # Binning sensitivity uses precisely the same primary pairs for every resolution.
  a <- m[primary$row1, , drop = FALSE]; b <- m[primary$row2, , drop = FALSE]
  bin_results <- rbindlist(lapply(c(5L, 10L, 20L), function(bins) {
    z <- pair_metrics(a, b, bins = bins, include_xi = FALSE)
    data.table(bins = bins, pairs = nrow(z), median_exact = median(z$agreement_exact),
      median_ordinal_similarity = median(z$ordinal_similarity),
      share_constant_levels_any = mean(z$constant_levels_any),
      spearman_similarity_vs_raw_MAE = cor(z$ordinal_similarity, -z$beta_MAE, method = "spearman"),
      max_MAE_quantization_error = max(abs(z$ordinal_L1 / bins - z$beta_MAE)))
  }))
  write_result(bin_results, "BinSensitivity.csv")

  # Leave-one-sample-out influence is a sensitivity diagnostic, not a confidence
  # interval and not independent-patient validation. No threshold is tuned on it.
  loo <- rbindlist(lapply(seq_len(18L), function(j) {
    z <- pair_metrics(a[, -j, drop = FALSE], b[, -j, drop = FALSE], include_xi = FALSE)
    data.table(omitted_sample = samples[j],
      median_ordinal_change = median(abs(z$ordinal_similarity - primary$ordinal_similarity)),
      q95_ordinal_change = unname(quantile(abs(z$ordinal_similarity - primary$ordinal_similarity), .95)),
      max_ordinal_change = max(abs(z$ordinal_similarity - primary$ordinal_similarity)),
      pearson_defined_both = sum(is.finite(z$pearson_raw) & is.finite(primary$pearson_raw)),
      median_pearson_change = median(abs(z$pearson_raw - primary$pearson_raw), na.rm = TRUE),
      q95_pearson_change = unname(quantile(abs(z$pearson_raw - primary$pearson_raw), .95, na.rm = TRUE)))
  }))
  write_result(loo, "LeaveOneSampleOut.csv")

  # Fixed, score-independent sample of primary pairs; NOT all-pairs discovery.
  # Each permutation applies one common patient-column shuffle to every second
  # profile, preserving inter-site structure while removing original alignment.
  set.seed(20260928L)
  selected <- sort(sample.int(nrow(primary), min(512L, nrow(primary))))
  pa <- a[selected, , drop = FALSE]; pb <- b[selected, , drop = FALSE]
  observed_stats <- pair_metrics(pa, pb)
  permutation_measures <- c("ordinal_similarity", "kappa_quadratic", "pearson_raw",
                             "spearman_levels", "xi_mean", "xi_max")
  permutations <- 199L
  null_means <- matrix(NA_real_, permutations, length(permutation_measures),
                        dimnames = list(NULL, permutation_measures))
  exceed <- matrix(0L, nrow(pa), length(permutation_measures), dimnames = list(NULL, permutation_measures))
  cat("Running", permutations, "patient permutations on", nrow(pa), "random primary pairs\n")
  for (r in seq_len(permutations)) {
    z <- pair_metrics(pa, pb[, sample.int(18L), drop = FALSE])
    for (name in permutation_measures) {
      obs <- observed_stats[[name]]; v <- z[[name]]
      # Pearson/Spearman are tested two-sided; the other statistics one-sided.
      if (name %in% c("pearson_raw", "spearman_levels")) {obs <- abs(obs); v <- abs(v)}
      null_means[r, name] <- mean(v, na.rm = TRUE)
      exceed[, name] <- exceed[, name] + as.integer(is.finite(obs) & v >= obs - 1e-12)
    }
    if (r %% 50L == 0L) cat("Permutation", r, "of", permutations, "\n")
  }
  perm_pair <- primary[selected, .(pos1, pos2, n_common, gap_bp)]
  perm_summary <- rbindlist(lapply(permutation_measures, function(name) {
    obs <- observed_stats[[name]]
    if (name %in% c("pearson_raw", "spearman_levels")) obs <- abs(obs)
    pp <- (exceed[, name] + 1) / (permutations + 1)
    pp[!is.finite(obs)] <- NA_real_
    perm_pair[, (paste0(name, "_permutation_p")) := pp]
    data.table(measure = name, sampled_pairs = length(obs), defined = sum(is.finite(obs)),
      observed_mean = mean(obs, na.rm = TRUE), null_mean = mean(null_means[, name]),
      null_mean_q025 = unname(quantile(null_means[, name], .025)),
      null_mean_q975 = unname(quantile(null_means[, name], .975)),
      aggregate_permutation_p = (1 + sum(null_means[, name] >= mean(obs, na.rm = TRUE) - 1e-12)) /
        (permutations + 1),
      pairs_unadjusted_p_le_05 = sum(pp <= .05, na.rm = TRUE),
      min_attainable_p = 1 / (permutations + 1))
  }))
  write_result(perm_summary, "PermutationDiagnostic.csv")
  write_result(perm_pair, "PermutationPairs.csv")
  write_result(as.data.table(null_means), "PermutationNullMeans.csv")

  # Finite-sample illustrative cases with known relationships, no truth labels
  # inferred from methylation states or existing outlier flags.
  base <- seq(.02, .98, length.out = 18)
  low <- seq(.01, .18, length.out = 18)
  cases <- list(
    identical_constant = list(rep(.05, 18), rep(.05, 18)),
    different_constants = list(rep(.05, 18), rep(.95, 18)),
    identical_varying = list(base, base),
    same_order_offset = list(low, low + .7),
    opposite_order = list(base, 1 - base),
    nonlinear_u = list(base, (2 * base - 1)^2),
    one_patient_change = list(rep(.05, 18), c(rep(.05, 17), .95)),
    close_across_bin_boundary = list(rep(.099, 18), rep(.101, 18)),
    # the supervisor's example, levels 1 1 1 1 2 1 1 and 1 2 1 1 2 1 1, as
    # betas inside those levels (0.05 = level 1, 0.15 = level 2)
    email_cg1_cg2 = list(c(.05, .05, .05, .05, .15, .05, .05),
                         c(.05, .15, .05, .05, .15, .05, .05)),
    email_cg1_constant = list(rep(.05, 7), c(.05, .15, .05, .05, .15, .05, .05)))
  examples <- rbindlist(lapply(names(cases), function(name) {
    x <- cases[[name]]
    cbind(data.table(example = name), pair_metrics(matrix(x[[1]], 1), matrix(x[[2]], 1)))
  }))
  write_result(examples, "IllustrativeCases.csv")

  # Standalone figure: observed data only, common denominators and clear panels.
  suppressPackageStartupMessages(library(ggplot2))
  plot_data <- rbindlist(lapply(c("agreement_exact", "ordinal_similarity", "pearson_raw", "xi_mean"),
    function(name) {
      z <- primary[, .(value = median(get(name), na.rm = TRUE)), by = gap_bin]
      z[, measure := name]; z
    }))
  plot_data[, measure := factor(measure,
    levels = c("agreement_exact", "ordinal_similarity", "pearson_raw", "xi_mean"),
    labels = c("Exact ten-level agreement", "Ordinal similarity (1 - L1/9)",
               "Pearson correlation (raw beta)", "Xi (mean, ties averaged)"))]
  g <- ggplot(plot_data, aes(gap_bin, value, group = measure, color = measure)) +
    geom_line(linewidth = .8) + geom_point(size = 2.8) +
    facet_wrap(~ measure, ncol = 2, scales = "fixed") +
    scale_y_continuous(limits = c(-.05, 1)) +
    labs(title = "Item 5: agreement and dependence of consecutive CpGs",
      subtitle = "32,664 true consecutive pairs; all 18 samples observed in both sites",
      x = "Genomic gap (bp)", y = "Median score",
      caption = "Scores have different meanings. Xi is not an agreement percentage.\nLong-gap complete-site bridges are excluded; large-gap bins have few pairs.") +
    theme_bw(base_size = 12) + theme(legend.position = "none",
      axis.text.x = element_text(angle = 30, hjust = 1),
      plot.title = element_text(face = "bold"), strip.background = element_rect(fill = "#eef3f8"))
  ggsave(file.path(output, "Task41_Item5_Similarity.png"), g, width = 10, height = 7, dpi = 180)

  manifest <- list(seed = 20260928L, main_n_common = 18L,
    coverage_sensitivity = c(6, 9, 12, 15, 18), xi_min_n_common = 12L,
    bin_sensitivity = c(5, 10, 20), permutation_pairs = nrow(pa), permutations = permutations,
    permutation_scope = "Diagnostic only; no multiplicity-adjusted discoveries or causal claims",
    xi_definition = "Exact average over predictor-tie rearrangements in each direction; mean/max of directional averages",
    seconds = unname(proc.time()[3] - started), R_version = R.version.string,
    packages = sapply(c("data.table", "matrixStats", "ggplot2", "jsonlite"), function(p) as.character(packageVersion(p))),
    script_md5 = unname(tools::md5sum(file.path(root, "Scripts", "Task41.RRBS.Item5.Feasibility.Sep28.2026.R"))))
  jsonlite::write_json(manifest, file.path(output, "Task41_RunManifest.json"), pretty = TRUE, auto_unbox = TRUE)
  capture.output(sessionInfo(), file = file.path(output, "Task41_SessionInfo.txt"))
  cat("Completed in", round(proc.time()[3] - started, 1), "seconds\n")
}

if (sys.nframe() == 0L) {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script <- normalizePath(sub("^--file=", "", file_arg[1L]))
  run_item5(dirname(dirname(script)))
}
