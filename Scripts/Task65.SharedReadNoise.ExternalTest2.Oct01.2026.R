# ================================================================
# Task 65 - Second external test of the shared-read noise correction,
#           on the GSE149438 plasma controls kept by the Task 62 screen
#           (a different lab, assay and pipeline from GSE233417)
# Date: Oct 1, 2026
#
# Runs Results/Task64_LockedProtocol.md exactly. The script stops if
# the protocol, its JSON, the code it names, the sample list or any
# input file has changed since the lock (SHA-256 in
# Results/Task64_LockedProtocol.json).
#
# STEPS
#   1. Pair universe (counts only, no outcome): adjacent genome CpGs
#      within 200 bp on autosomes where at least 20 donors have at
#      least 6 reads at both sites.
#   2. Three-way split of every donor's molecules (seed split_seed +
#      1000 x donor + chromosome); A-only estimates; B/C reference.
#   3. Primary endpoint: mean dC = (Chat_A - R)^2 - (s_A - R)^2 over
#      pairs with 20+ eligible donors. 95% interval: 2,000 joint
#      bootstrap draws of donors and 1-Mb blocks within chromosomes,
#      bias-corrected percentile (version 2, Task 57), on the bootstrap
#      pair set: all pairs if there are at most 300,000, otherwise
#      every pair whose CpG index is divisible by 25. Primary
#      criterion: upper end below 0.
#   4. Secondary: disagreement endpoint; strata by distance, depth and
#      overlap; the declared sensitivity analysis on pairs at most
#      40 bp apart (Task 61: mates whose CpG ranges do not overlap are
#      two lines in .mhap), with its own interval from the same draws.
#   5. Section 6 comparisons with the Task 63 definitions: correlation
#      recovery, pair ranking (on the bootstrap pair set) and the
#      Ding & Gentleman baseline (2,000 pairs). Descriptive only.
#
# DRY RUN: `R CMD BATCH "--args --dry-run-colon" <this file>` runs
#   steps 1 to 5 on the GTEx colon development data instead, with no
#   lock check and outputs under Data/gtex_colon_rrbs/task65_dry_run/.
#   It uses chromosomes 21 and 22 only, and exists to test the code
#   before the lock.
#
# OUTPUTS
#   Data/gse149438_plasma_mhap/s6_pair_results.rds
#   Results/Task65_Primary.csv, Task65_ByStratum.csv,
#   Task65_CorrelationRecovery.csv, Task65_Ranking.csv,
#   Task65_MEcorBaseline.csv, Fig46_ExternalTest2.png
# ================================================================

suppressPackageStartupMessages({ library(data.table); library(jsonlite); library(ggplot2); library(MeasurementError.cor) })
pick_dir <- function(...) { for (d in c(...)) if (dir.exists(d)) return(d)
  stop("no candidate directory exists") }
repo_dir <- pick_dir("/mmfs1/home/wln26/Experiments.Outlier.July31.2026",
                     path.expand("~/Desktop/bioinformatics-research"))
DRY <- "--dry-run-colon" %in% commandArgs(TRUE)
res <- file.path(repo_dir, "Results")

# ----------------------------------------------------------------
# 0. the lock is intact (or the dry run's stand-in parameters)
# ----------------------------------------------------------------
if (!DRY) {
  lock <- fromJSON(file.path(res, "Task64_LockedProtocol.json"))
  for (f in names(lock$sha256))
    stopifnot(digest::digest(file = file.path(repo_dir, f), algo = "sha256") == lock$sha256[[f]])
  P <- lock$parameters; dat <- file.path(repo_dir, "Data", "gse149438_plasma_mhap"); outp <- res; tag <- "Task65"
} else {
  P <- list(max_gap_bp = 200L, min_reads_per_part = 2L, min_donors = 20L, split_seed = 20261002L, bootstrap_seed = 20261003L,
            bootstrap_draws = 200L, bootstrap_all_pairs_max = 300000L, bootstrap_subsample_every = 25L,
            sensitivity_max_gap_bp = 40L, xi_seed = 20261205L, me_seed = 20261204L, me_pairs = 300L)
  dat <- file.path(repo_dir, "Data", "gtex_colon_rrbs"); outp <- file.path(dat, "task65_dry_run"); tag <- "DryRun"
  dir.create(outp, showWarnings = FALSE)
}
for (f in c("SharedReadNoise.Core.R", "SharedReadNoise.Benchmark.R", "Task47.FrozenScore.Sep29.2026.R",
            "SharedReadNoise.Section6.R")) source(file.path(repo_dir, "Scripts", f))
AUTO <- paste0("chr", 1:22); s6 <- function(x) signif(x, 4); r4 <- function(x) round(x, 4)
CHR <- if (DRY) 21:22 else 1:22                   # the dry run uses two chromosomes, to be quick
cat(if (DRY) "DRY RUN on GTEx colon (development data)\n" else "lock verified\n")

cpg <- readRDS(file.path(repo_dir, "Data", "hg19_seq", "cpg_positions_hg19.rds"))
I0 <- c(0L, cumsum(lengths(cpg))[-22])
if (!DRY) {
  kept <- fread(file.path(res, "Task62_PlasmaSamples.csv"))[kept == TRUE][order(gsm)]
  files <- file.path(dat, paste0(kept$srx, ".mhap.gz")); stopifnot(length(files) == lock$donors)
  read_donor <- function(f) {                     # .mhap -> .pat-like rows by chromosome (Task 62 reader)
    m <- fread(cmd = sprintf("gzip -dc '%s'", f), header = FALSE, sep = "\t",
               col.names = c("chr", "start", "end", "hap", "n", "strand"),
               colClasses = c("character", "integer", "integer", "character", "integer", "character"))[chr %in% AUTO]
    m[, k := match(chr, AUTO)]
    m[, li := { p <- cpg[[k[1]]]; match(start, p) }, by = k]; m[, lj := { p <- cpg[[k[1]]]; match(end, p) }, by = k]
    if (m[is.na(li) | is.na(lj) | (lj - li + 1L) != nchar(hap), .N]) stop(basename(f), ": lines off the hg19 CpG list")
    p <- m[, .(chr, idx = I0[k] + li, pat = chartr("01", "TC", hap), n)]
    split(p[, .(idx, pat, n)], p$chr)
  }
} else {
  samples <- fromJSON(file.path(dat, "colon_samples.json"))
  files <- file.path(dat, basename(sapply(samples$files, `[`, 2)))
  read_donor <- function(f) {
    p <- fread(cmd = sprintf("gzip -dc '%s'", f), header = FALSE, sep = "\t", col.names = c("chr", "idx", "pat", "n"),
               colClasses = c("character", "integer", "character", "integer"))
    split(p[chr %in% AUTO, .(idx, pat, n)], p$chr[p$chr %in% AUTO])
  }
}
D <- length(files); stopifnot(D >= 20, all(file.exists(files)))

# ----------------------------------------------------------------
# 1. pair universe from counts only
# ----------------------------------------------------------------
lookup <- lapply(cpg, function(p) { ok <- which(diff(p) <= P$max_gap_bp); v <- rep(NA_integer_, length(p)); v[ok] <- seq_along(ok); v })
deep <- lapply(lookup, function(v) integer(sum(!is.na(v))))
t0 <- Sys.time()
for (s in seq_len(D)) {
  ps <- read_donor(files[s])
  for (k in CHR) {
    pk <- ps[[AUTO[k]]]; if (is.null(pk) || !nrow(pk)) next
    pr <- srn_pat_counts(pk, I0[k], length(cpg[[k]]), cpg[[k]], max_gap = P$max_gap_bp)
    pr <- pr[Ni >= 3L * P$min_reads_per_part & Nj >= 3L * P$min_reads_per_part]
    if (!nrow(pr)) next
    ix <- lookup[[k]][pr$idx - I0[k]]; deep[[k]][ix] <- deep[[k]][ix] + 1L
  }
  rm(ps); invisible(gc())
}
uni <- rbindlist(lapply(1:22, function(k) { r <- which(!is.na(lookup[[k]]))[deep[[k]] >= P$min_donors]
  data.table(chrn = k, idx = I0[k] + r, pos = cpg[[k]][r], gap = cpg[[k]][r + 1L] - cpg[[k]][r]) }))
setkey(uni, chrn, idx); rm(deep, lookup); invisible(gc())
boot_all <- nrow(uni) <= P$bootstrap_all_pairs_max
cat("pair universe:", nrow(uni), "pairs (", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min ); bootstrap on",
    if (boot_all) "all pairs\n" else paste("every", P$bootstrap_subsample_every, "th pair\n"))

# ----------------------------------------------------------------
# 2. split, count and accumulate
# ----------------------------------------------------------------
acc <- lapply(1:22, function(k) s6_accumulator(uni[chrn == k, .N]))
uidx <- lapply(1:22, function(k) uni[chrn == k, idx])
sub <- if (boot_all) uni[, .(chrn, idx, pos, gap)] else uni[idx %% P$bootstrap_subsample_every == 0L, .(chrn, idx, pos, gap)]
SM <- c("E", "xA", "yA", "vxA", "vyA", "cA", "biB", "bjB", "biC", "bjC")
subm <- setNames(lapply(SM, function(v) matrix(0, nrow(sub), D)), SM); sub_key <- paste(sub$chrn, sub$idx)
for (s in seq_len(D)) {
  ps <- read_donor(files[s])
  for (k in CHR) {
    pk <- ps[[AUTO[k]]]; if (is.null(pk) || !nrow(pk)) next
    m <- srb_donor_chr(pk, I0[k], length(cpg[[k]]), cpg[[k]], seed = P$split_seed + 1000L * s + k, max_gap = P$max_gap_bp)
    if (is.null(m)) next
    rows <- match(m$idx, uidx[[k]]); keep <- !is.na(rows); if (!any(keep)) next
    m <- m[keep]; s6_update(acc[[k]], m, rows[keep])
    at <- match(paste(k, m$idx), sub_key); h <- which(!is.na(at))
    if (length(h)) { subm$E[at[h], s] <- 1; for (v in SM[-1]) subm[[v]][at[h], s] <- m[[v]][h] }
  }
  rm(ps); invisible(gc())
  cat(sprintf("donor %2d/%d done (%.1f min)\n", s, D, as.numeric(Sys.time() - t0, units = "mins")))
}
est <- rbindlist(lapply(1:22, function(k) {
  a <- copy(acc[[k]])[, `:=`(chrn = k, idx = uidx[[k]])]; keep <- a$n >= P$min_donors
  cbind(a[keep, .(chrn, idx)], s6_finish(a[keep], min_donors = P$min_donors))
}))
est <- uni[est, on = c("chrn", "idx")]
saveRDS(list(est = est, sub = sub, subm = subm, built = format(Sys.time())), file.path(dat, if (DRY) "task65_dry_run/s6_pair_results.rds" else "s6_pair_results.rds"))

# ----------------------------------------------------------------
# 3. primary endpoint and its interval; 4. the sensitivity subset
# ----------------------------------------------------------------
ok <- rowSums(subm$E) >= P$min_donors
E <- subm$E[ok, ]; sb <- sub[ok]; near <- sb$gap <= P$sensitivity_max_gap_bp
pre <- list(EX = E * subm$xA[ok, ], EY = E * subm$yA[ok, ], EXY = E * subm$xA[ok, ] * subm$yA[ok, ], EC = E * subm$cA[ok, ],
            EBI = E * subm$biB[ok, ], ECJ = E * subm$bjC[ok, ], EBIC = E * subm$biB[ok, ] * subm$bjC[ok, ],
            ECI = E * subm$biC[ok, ], EBJ = E * subm$bjB[ok, ], ECIB = E * subm$biC[ok, ] * subm$bjB[ok, ])
dC_w <- function(w) {
  n <- drop(E %*% w); g <- function(M) drop(M %*% w); cv <- function(sxy, sx, sy) (sxy - sx * sy / n) / (n - 1)
  W <- cv(g(pre$EXY), g(pre$EX), g(pre$EY)); U <- W - g(pre$EC) / n
  R <- (cv(g(pre$EBIC), g(pre$EBI), g(pre$ECJ)) + cv(g(pre$ECIB), g(pre$ECI), g(pre$EBJ))) / 2
  d <- (U - R)^2 - (W - R)^2; d[n < 2] <- NA_real_; d
}
blk <- paste(sb$chrn, (sb$pos - 1L) %/% 1e6L); blocks <- unique(data.table(chrn = sb$chrn, blk))
set.seed(P$bootstrap_seed)
bt <- vapply(seq_len(P$bootstrap_draws), function(b) {
  w <- tabulate(sample.int(D, D, replace = TRUE), D)
  bm <- blocks[, .(blk = sample(blk, .N, replace = TRUE)), by = chrn][, .N, by = blk]
  bw <- bm$N[match(blk, bm$blk)]; bw[is.na(bw)] <- 0
  d <- dC_w(w); k <- !is.na(d) & bw > 0; kn <- k & near
  c(sum(bw[k] * d[k]) / sum(bw[k]), sum(bw[kn] * d[kn]) / sum(bw[kn]), sum(is.na(d)))
}, numeric(3))
d0 <- dC_w(rep(1, D)); sub_est <- mean(d0, na.rm = TRUE); near_est <- mean(d0[near], na.rm = TRUE)
bias <- mean(bt[1, ]) - sub_est; bias_n <- mean(bt[2, ]) - near_est
prim <- data.table(donors = D, pairs = nrow(est), delta_mse_cov = s6(mean(est$dC)), bootstrap_pairs = nrow(sb),
                   bootstrap_on_all_pairs = boot_all, subsample_estimate = s6(sub_est),
                   lo95 = s6(quantile(bt[1, ], 0.025) - bias), hi95 = s6(quantile(bt[1, ], 0.975) - bias),
                   plain_percentile_lo95 = s6(quantile(bt[1, ], 0.025)), plain_percentile_hi95 = s6(quantile(bt[1, ], 0.975)),
                   relative_change = r4(mean(est$dC) / mean(est$loss_raw)), share_cor_closer = r4(mean(est$loss_cor < est$loss_raw)),
                   delta_mse_d2 = s6(mean(est$dD)), mean_failed_pairs_per_draw = r4(mean(bt[3, ])),
                   near_pairs = sum(est$gap <= P$sensitivity_max_gap_bp),
                   near_delta_mse_cov = s6(mean(est[gap <= P$sensitivity_max_gap_bp, dC])),
                   near_lo95 = s6(quantile(bt[2, ], 0.025) - bias_n), near_hi95 = s6(quantile(bt[2, ], 0.975) - bias_n))
prim[, verdict := fcase(hi95 < 0, "primary criterion met: correction closer to the independent reference",
                        lo95 > 0, "correction has higher error", default = "improvement not established")]
fwrite(prim, file.path(outp, paste0(tag, "_Primary.csv")))
cat("\nPRIMARY (locked):\n"); print(prim)

endpoint <- function(d) d[, .(pairs = .N, delta_mse_cov = s6(mean(dC)), relative_change = r4(mean(dC) / mean(loss_raw)),
                              share_cor_closer = r4(mean(loss_cor < loss_raw)), delta_mse_d2 = s6(mean(dD)),
                              mean_noise_A = s6(mean(mean_noise_A)), mean_raw_cov_A = s6(mean(raw_cov_A)),
                              mean_ref_cov = s6(mean(ref_cov)), mean_cor_cov_A = s6(mean(cor_cov_A)))]
by_one <- function(name, lev) { d <- copy(est)[, level := lev]
  d[, endpoint(.SD), by = level][order(level)][, `:=`(stratum = name, level = as.character(level))] }
strata <- rbindlist(list(endpoint(est)[, `:=`(stratum = "all pairs", level = "all")],
  by_one("distance (bp)", cut(est$gap, c(0, 10, 20, 40, 60, 100, 150, 200))),
  by_one("mean depth in A (reads)", cut(est$depth_A, c(2, 5, 10, 20, Inf), right = FALSE)),
  by_one("overlap in A", cut(est$overlap_A, c(0, 0.25, 0.5, 0.75, 1.0001), right = FALSE))), use.names = TRUE)
setcolorder(strata, c("stratum", "level"))
fwrite(strata, file.path(outp, paste0(tag, "_ByStratum.csv")))
cat("\nSecondary, by stratum:\n"); print(strata[, .(stratum, level, pairs, delta_mse_cov, relative_change, share_cor_closer)])

# ----------------------------------------------------------------
# 5. section 6 comparisons (Task 63 definitions; descriptive)
# ----------------------------------------------------------------
corr <- s6_corr_table(est); fwrite(corr, file.path(outp, paste0(tag, "_CorrelationRecovery.csv")))
cat("\nCorrelation recovery:\n"); print(corr[, lapply(.SD, function(x) if (is.numeric(x)) s6(x) else x)])
sp <- est[sub, on = c("chrn", "idx"), nomatch = NA]; okr <- !is.na(sp$n)
scores <- s6_scores(lapply(subm, function(M) M[okr, , drop = FALSE]), sp$gap[okr], sp$raw_d2_A[okr], sp$cor_d2_A[okr], P$xi_seed)
rk <- s6_rank_table(scores, sp$ref_d2[okr], varies = sp$var_x_A[okr] > 0 & sp$var_y_A[okr] > 0)
fwrite(rk, file.path(outp, paste0(tag, "_Ranking.csv")))
cat("\nPair ranking, mean reference squared disagreement of the top pairs:\n")
print(dcast(rk, set + label ~ top_fraction, value.var = "mean_qref_top", fun.aggregate = function(x) s6(x)))
set.seed(P$me_seed)
cand <- which(okr & !is.na(sp$rho_ref)); pick <- sort(sample(cand, min(P$me_pairs, length(cand))))
me <- rbindlist(lapply(pick, function(i) {
  e <- subm$E[i, ] == 1
  fit <- tryCatch(cor.me.vector(subm$xA[i, e], sqrt(subm$vxA[i, e]), subm$yA[i, e], sqrt(subm$vyA[i, e])), error = function(x) NULL)
  data.table(row = i, me = if (is.null(fit) || fit$convergence != 0) NA_real_ else unname(fit$estimate["corr.true"]))
}))
me[, `:=`(rho_ref = sp$rho_ref[row], r_pearson = sp$r_pearson[row], r_vo = sp$r_vo[row], r_full = sp$r_full[row])]
v <- me[!is.na(me) & !is.na(r_pearson) & !is.na(r_vo) & !is.na(r_full)]
mtab <- data.table(estimator = c("Pearson", "variance-only", "full correction", "Ding & Gentleman (MeasurementError.cor)"),
                   pairs_tried = nrow(me), valid_share = c(mean(!is.na(me$r_pearson)), mean(!is.na(me$r_vo)), mean(!is.na(me$r_full)),
                                                           mean(!is.na(me$me))), common_pairs = nrow(v),
                   mse_common = c(mean((v$r_pearson - v$rho_ref)^2), mean((v$r_vo - v$rho_ref)^2), mean((v$r_full - v$rho_ref)^2),
                                  mean((v$me - v$rho_ref)^2)),
                   mean_error_common = c(mean(v$r_pearson - v$rho_ref), mean(v$r_vo - v$rho_ref), mean(v$r_full - v$rho_ref),
                                         mean(v$me - v$rho_ref)))
fwrite(mtab, file.path(outp, paste0(tag, "_MEcorBaseline.csv")))
cat("\nDing & Gentleman baseline:\n"); print(mtab[, lapply(.SD, function(x) if (is.numeric(x)) s6(x) else x)])

# ----------------------------------------------------------------
# figure
# ----------------------------------------------------------------
SURF <- "#fcfcfb"; INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"; GRID <- "#e1e0d9"; BLUE <- "#2a78d6"
f46 <- strata[stratum != "all pairs"][, `:=`(stratum = factor(stratum, levels = unique(stratum)), level = factor(level, levels = unique(level)))]
p46 <- ggplot(f46, aes(relative_change, level)) +
  geom_vline(xintercept = 0, colour = MUTED, linewidth = 0.4) + geom_point(colour = BLUE, size = 3) +
  facet_grid(stratum ~ ., scales = "free_y", space = "free_y") +
  labs(title = if (DRY) "DRY RUN on GTEx colon (development data)" else
         "Second external test (locked): is the corrected covariance closer to an independent reference?",
       subtitle = paste0(D, if (DRY) " colon donors. " else " plasma cfDNA controls (GSE149438), another lab and pipeline. ",
                         "Primary: ", prim$verdict, ".\nDelta MSE ", prim$delta_mse_cov, ", 95% interval ", prim$lo95, " to ", prim$hi95, "."),
       x = "Relative change in squared error, corrected vs observed (left of 0 = correction helps)", y = NULL,
       caption = paste0("Source: ", tag, "_Primary.csv and ", tag, "_ByStratum.csv (Task 65)")) +
  theme_minimal(base_size = 10) +
  theme(plot.background = element_rect(fill = SURF, colour = NA), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = GRID, linewidth = 0.3),
        strip.text.y = element_text(angle = 0, hjust = 0, face = "bold", colour = INK),
        plot.title = element_text(colour = INK, face = "bold", size = 12), plot.title.position = "plot",
        plot.subtitle = element_text(colour = INK2, size = 9), plot.caption = element_text(colour = MUTED, hjust = 0),
        plot.caption.position = "plot", axis.text = element_text(colour = INK2))
ggsave(file.path(outp, if (DRY) "DryRun_Fig.png" else "Fig46_ExternalTest2.png"), p46, width = 10, height = 6, dpi = 200, bg = SURF)
cat("\ndone\n")
