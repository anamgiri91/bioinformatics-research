# ================================================================
# Task 53 - Shared-read noise on real data: how much of the
#           co-methylation between close CpGs is read noise?
# Date: Oct 1, 2026
#
# EXPLORATORY DEVELOPMENT WORK on the GTEx colon cohort, which has
# already been examined. This is the core of plan_shared_read_noise.md
# (count rebuild, estimator, development diagnostics). It is not the
# plan's external validation, and nothing here is tuned for one.
#
# THE IDEA
#   When two close CpGs are read on the same DNA fragment, their
#   measured betas in one person share read noise. Across people, that
#   adds a noise term to the observed covariance:
#     observed covariance = true covariance + mean shared-read noise
#   The noise term can be estimated in each person from the reads that
#   call both sites (n00, n01, n10, n11) and subtracted. The estimator
#   is Scripts/SharedReadNoise.Core.R, unchanged (300 checks pass).
#
# DATA
#   GSE233417 colon RRBS, 29 donors, hg19. GEO's processing notes say
#   PCR duplicates were removed with UMIs and the .pat files hold
#   fragment-level calls from the deduplicated reads, so each counted
#   unit should be a distinct molecule. The raw reads are in EGA under
#   controlled access, so this cannot be checked directly here.
#   Pairs: every autosomal pair of adjacent genome CpGs within 200 bp.
#   A donor counts for a pair when both sites have at least 2 reads and
#   the shared-read count is 0 or at least 2; a pair needs 20 donors.
#
# QUESTIONS
#   1. How big is the shared-read noise term, by distance, depth,
#      overlap and spread?
#   2. The open caveat of Task 52 (H9): low-spread pairs with strong
#      read linkage predict a held-out patient. Is their covariance
#      shared biology, or mostly shared read noise?
#
# OUTPUTS
#   Data/gtex_colon_rrbs/srn_pair_estimates.rds   per-pair estimates
#   Results/Task53_Coverage.csv          pairs and donor exclusions
#   Results/Task53_ByStratum.csv         noise term by distance, depth,
#                                        overlap and spread
#   Results/Task53_H9_Link.csv           the Task 52 H9 groups
#   Results/Fig39_SharedReadNoise.png
#
# LIMITS
#   Unbiased only under the plan's model: reads are independent
#   molecules, and shared reads represent each site like the others.
#   Intervals resample 1-Mb blocks of pairs, not donors, so they
#   understate donor uncertainty.
# ================================================================

suppressPackageStartupMessages({ library(data.table); library(jsonlite); library(ggplot2) })

pick_dir <- function(...) { for (d in c(...)) if (dir.exists(d)) return(d)
  stop("no candidate directory exists") }
repo_dir <- pick_dir("/mmfs1/home/wln26/Experiments.Outlier.July31.2026",
                     path.expand("~/Desktop/bioinformatics-research"))
res <- file.path(repo_dir, "Results"); dat <- file.path(repo_dir, "Data", "gtex_colon_rrbs")
source(file.path(repo_dir, "Scripts", "SharedReadNoise.Core.R"))
AUTO <- paste0("chr", 1:22); MAX_GAP <- 200L; MIN_DONORS <- 20L
r4 <- function(x) round(x, 4); s6 <- function(x) signif(x, 4)

cpg <- readRDS(file.path(repo_dir, "Data", "hg19_seq", "cpg_positions_hg19.rds"))
I0 <- c(0L, cumsum(lengths(cpg))[-22])
samples <- fromJSON(file.path(dat, "colon_samples.json"))
sid <- sub("^GSM[0-9]+_(GTEX-[A-Z0-9]+)-.*", "\\1", basename(sapply(samples$files, `[`, 1)))
stopifnot(length(sid) == 29, !anyDuplicated(sid))

# ----------------------------------------------------------------
# 1. rebuild joint counts donor by donor and stream them into the
#    estimator (cached)
# ----------------------------------------------------------------
est_file <- file.path(dat, "srn_pair_estimates.rds")
if (!file.exists(est_file)) {
  # accumulator row for each adjacent pair within 200 bp, by CpG rank of site i
  lookup <- lapply(cpg, function(p) { ok <- which(diff(p) <= MAX_GAP)
    v <- rep(NA_integer_, length(p)); v[ok] <- seq_along(ok); v })
  acc <- lapply(lookup, function(v) srn_accumulator(sum(!is.na(v))))
  t0 <- Sys.time()
  for (s in seq_along(sid)) {
    fp <- file.path(dat, basename(samples$files[[s]][2]))
    p <- fread(cmd = sprintf("gzip -dc '%s'", fp), header = FALSE, sep = "\t",
               col.names = c("chr", "idx", "pat", "n"),
               colClasses = c("character", "integer", "character", "integer"))
    ps <- split(p[chr %in% AUTO, .(idx, pat, n)], p$chr[p$chr %in% AUTO]); rm(p)
    for (k in 1:22) {
      pk <- ps[[AUTO[k]]]
      if (is.null(pk) || !nrow(pk)) next
      pr <- srn_pat_counts(pk, I0[k], length(cpg[[k]]), cpg[[k]], max_gap = MAX_GAP)
      if (!nrow(pr)) next
      ix <- lookup[[k]][pr$idx - I0[k]]
      stopifnot(!anyNA(ix))
      srn_update(acc[[k]], pr, ix)
    }
    rm(ps); invisible(gc())
    cat(sprintf("donor %2d/29 %s done (%.1f min)\n", s, sid[s], as.numeric(Sys.time() - t0, units = "mins")))
  }
  cover <- rbindlist(lapply(1:22, function(k) acc[[k]][, .(chrn = k, pairs_within_200bp = .N,
    pairs_seen = sum(observed > 0), pairs_with_20_donors = sum(n >= MIN_DONORS),
    donor_rows_seen = sum(observed), donor_rows_low_coverage = sum(low_coverage),
    donor_rows_single_shared = sum(single_shared), donor_rows_eligible = sum(n),
    donor_rows_no_shared = sum(zero_shared))]))
  est <- rbindlist(lapply(1:22, function(k) {
    f <- srn_finish(acc[[k]], min_donors = MIN_DONORS)
    r <- which(!is.na(lookup[[k]]))[f$index]                   # CpG rank of site i
    f[, `:=`(chrn = k, idx = I0[k] + r, pos = cpg[[k]][r], gap = cpg[[k]][r + 1L] - cpg[[k]][r])]
  }))
  rm(acc); invisible(gc())
  saveRDS(list(est = est, cover = cover, built = format(Sys.time())), est_file)
}
z <- readRDS(est_file); est <- z$est; cover <- z$cover
cat("pairs with an estimate (20+ donors):", nrow(est), "\n")

cov_tab <- cover[, lapply(.SD, sum), .SDcols = -"chrn"]
cov_tab[, `:=`(share_donor_rows_low_coverage = r4(donor_rows_low_coverage / donor_rows_seen),
               share_donor_rows_single_shared = r4(donor_rows_single_shared / donor_rows_seen),
               share_eligible_rows_no_shared = r4(donor_rows_no_shared / donor_rows_eligible))]
fwrite(cov_tab, file.path(res, "Task53_Coverage.csv"))
cat("\nCoverage:\n"); print(t(cov_tab))

# ----------------------------------------------------------------
# 2. how big is the noise term, and where?
# ----------------------------------------------------------------
est[, `:=`(depth = (mean_Ni + mean_Nj) / 2, overlap = mean_shared / pmin(mean_Ni, mean_Nj),
           spread = (sqrt(pmax(raw_var_i, 0)) + sqrt(pmax(raw_var_j, 0))) / 2)]
summ <- function(d) d[, .(pairs = .N, mean_raw_cov = s6(mean(raw_cov)), mean_noise_cov = s6(mean(mean_noise_cov)),
                          mean_corrected_cov = s6(mean(cov_corrected)),
                          noise_share_of_raw_cov = r4(sum(mean_noise_cov) / sum(raw_cov)),
                          share_corrected_rho_valid = r4(mean(rho_status == "valid")),
                          share_corrected_rho_outside = r4(mean(rho_status == "outside_unit_interval")),
                          share_corrected_variance_nonpositive = r4(mean(rho_status == "nonpositive_corrected_variance")),
                          median_raw_pearson = r4(median(raw_pearson, na.rm = TRUE)),
                          median_corrected_rho_if_valid = r4(median(rho_corrected[rho_status == "valid"])),
                          mean_raw_d2 = s6(mean(raw_d2)), mean_corrected_d2 = s6(mean(d2_corrected)))]
by_one <- function(name, levels_expr) {
  d <- copy(est)[, level := eval(levels_expr)]
  out <- d[, summ(.SD), by = level][order(level)]
  out[, `:=`(stratum = name, level = as.character(level))]
}
strata <- rbindlist(list(
  summ(est)[, `:=`(stratum = "all pairs", level = "all")],
  by_one("distance (bp)", quote(cut(gap, c(0, 10, 20, 40, 60, 100, 150, 200)))),
  by_one("mean depth (reads)", quote(cut(depth, c(2, 10, 20, 50, Inf), right = FALSE))),
  by_one("overlap (shared / smaller depth)", quote(cut(overlap, c(0, 0.25, 0.5, 0.75, 1.0001), right = FALSE))),
  by_one("spread (mean raw SD)", quote(cut(spread, c(0, 0.02, 0.05, 0.1, Inf), right = FALSE)))), use.names = TRUE)
setcolorder(strata, c("stratum", "level"))
fwrite(strata, file.path(res, "Task53_ByStratum.csv"))
cat("\nNoise term by stratum:\n"); print(strata[, .(stratum, level, pairs, mean_raw_cov, mean_noise_cov,
                                                     mean_corrected_cov, noise_share_of_raw_cov,
                                                     share_corrected_rho_valid)], nrows = 60)

# ----------------------------------------------------------------
# 3. the Task 52 H9 groups: shared biology or shared read noise?
# ----------------------------------------------------------------
po <- readRDS(file.path(dat, "pair_outcomes.rds"))
lk <- readRDS(file.path(dat, "cache_cohort.rds"))$linkage[, .(chrn, idx, linkage)]
h9 <- merge(merge(po[low_spread == TRUE & gap <= 200], lk, by = c("chrn", "idx")),
            est[, .(chrn, idx, raw_cov, mean_noise_cov, cov_corrected, n)], by = c("chrn", "idx"))
h9[, above := linkage > median(linkage)]
h9[, block := paste(chrn, pos %/% 1e6L)]
set.seed(20261001L)
block_ci <- function(d, num, den = NULL) {      # 95% interval by resampling 1-Mb blocks (exact for sums)
  S <- rowsum(cbind(d[[num]], if (is.null(den)) 1 else d[[den]]), d$block, reorder = FALSE); nb <- nrow(S)
  v <- replicate(1000, { w <- colSums(S[sample.int(nb, nb, replace = TRUE), , drop = FALSE]); w[1] / w[2] })
  unname(quantile(v, c(0.025, 0.975)))
}
link <- rbindlist(lapply(c(FALSE, TRUE), function(a) {
  d <- h9[above == a]
  ci_c <- block_ci(d, "cov_corrected"); ci_s <- block_ci(d, "mean_noise_cov", "raw_cov")
  data.table(group = if (a) "above-median linkage" else "at or below median linkage", pairs = nrow(d),
             heldout_gain = s6(mean(d$heldout_gain)), mean_raw_cov = s6(mean(d$raw_cov)),
             mean_noise_cov = s6(mean(d$mean_noise_cov)), mean_corrected_cov = s6(mean(d$cov_corrected)),
             corrected_cov_lo95 = s6(ci_c[1]), corrected_cov_hi95 = s6(ci_c[2]),
             noise_share_of_raw_cov = r4(sum(d$mean_noise_cov) / sum(d$raw_cov)),
             noise_share_lo95 = r4(ci_s[1]), noise_share_hi95 = r4(ci_s[2]))
}))
fwrite(link, file.path(res, "Task53_H9_Link.csv"))
cat("\nThe Task 52 H9 groups (low-spread pairs within 200 bp):\n"); print(link)
al <- merge(merge(po[gap <= 200], lk, by = c("chrn", "idx")),
            est[, .(chrn, idx, raw_cov, mean_noise_cov, cov_corrected)], by = c("chrn", "idx"))
sp <- al[, .(pairs = .N, raw_cov = cor(heldout_gain, raw_cov, method = "spearman"),
             noise_cov = cor(heldout_gain, mean_noise_cov, method = "spearman"),
             corrected_cov = cor(heldout_gain, cov_corrected, method = "spearman"))]
fwrite(sp, file.path(res, "Task53_GainCorrelations.csv"))
cat("\nAll cohort pairs within 200 bp: Spearman of held-out gain with each covariance:\n"); print(sp)

# ----------------------------------------------------------------
# 4. figure
# ----------------------------------------------------------------
SURF <- "#fcfcfb"; INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"; GRID <- "#e1e0d9"
COL <- c(`Observed covariance` = "#52514e", `Shared-read noise` = "#eb6834", `Corrected covariance` = "#2a78d6")
dist <- strata[stratum == "distance (bp)", .(level = factor(level, levels = level), `Observed covariance` = mean_raw_cov,
                                             `Shared-read noise` = mean_noise_cov, `Corrected covariance` = mean_corrected_cov)]
dl <- melt(dist, id.vars = "level", variable.name = "measure", value.name = "value")
pA <- ggplot(dl, aes(level, value, colour = measure, group = measure)) +
  geom_hline(yintercept = 0, colour = MUTED, linewidth = 0.4) +
  geom_line(linewidth = 0.7) + geom_point(size = 2.6) +
  scale_colour_manual(values = COL, name = NULL) +
  labs(title = "A. All pairs: mean covariance across donors, by distance",
       x = "Distance between the two CpGs (bp)", y = "Mean covariance across donors")
gl <- melt(link[, .(group, `Observed covariance` = mean_raw_cov, `Shared-read noise` = mean_noise_cov,
                    `Corrected covariance` = mean_corrected_cov)], id.vars = "group", variable.name = "measure")
gl[, group := factor(group, levels = c("at or below median linkage", "above-median linkage"),
                     labels = c("Below-median linkage", "Above-median linkage"))]
pB <- ggplot(gl, aes(value, group, colour = measure)) +
  geom_vline(xintercept = 0, colour = MUTED, linewidth = 0.4) +
  geom_point(size = 3.2, alpha = 0.9) +
  scale_colour_manual(values = COL, name = NULL) +
  labs(title = "B. Low-spread pairs from Task 52 (H9)", x = "Mean covariance across donors", y = NULL)
th <- theme_minimal(base_size = 10.5) +
  theme(plot.background = element_rect(fill = SURF, colour = NA), panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = GRID, linewidth = 0.3), legend.position = "top",
        legend.justification = "left", plot.title = element_text(colour = INK, face = "bold", size = 11.5),
        axis.text = element_text(colour = INK2), axis.title = element_text(colour = INK2))
png(file.path(res, "Fig39_SharedReadNoise.png"), width = 11.5, height = 4.8, units = "in", res = 200, bg = SURF)
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 2, heights = grid::unit(c(0.11, 0.89), "npc"))))
grid::grid.text("How much of the covariance between close CpGs is shared-read noise? GTEx colon RRBS, 29 donors",
                x = 0.01, hjust = 0, gp = grid::gpar(fontface = "bold", fontsize = 13, col = INK),
                vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1:2))
print(pA + th, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
print(pB + th, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 2))
invisible(dev.off())
cat("\ndone\n")
