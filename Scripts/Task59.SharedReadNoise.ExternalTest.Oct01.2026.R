# ================================================================
# Task 59 - The external test of the shared-read noise correction
#           (plan phase G), on the blood samples screened in Task 56
# Date: Oct 1, 2026
#
# Runs Results/Task58_LockedProtocol.md exactly. The script stops if
# the protocol, its JSON or the code it names have changed since the
# lock (SHA-256 in Results/Task58_LockedProtocol.json).
#
# STEPS
#   1. Pair universe (counts only, no outcome): adjacent genome CpGs
#      within 200 bp on autosomes where at least 20 donors have at
#      least 6 reads at both sites (necessary for the split rule).
#   2. Three-way split of every donor's fragments (seed 20261102 +
#      1000 x donor + chromosome); A-only estimates; B/C reference;
#      per-pair endpoint dC = (Chat_A - R)^2 - (s_A - R)^2.
#   3. Primary endpoint: mean dC over pairs with 20+ eligible donors.
#      95% interval: 2,000 joint bootstrap draws of donors and 1-Mb
#      blocks within chromosomes on the fixed 1-in-25 pair subsample
#      (seed 20261103), bias-corrected percentile interval (protocol
#      version 2, chosen in the Task 57 calibration). Primary
#      criterion: upper end below 0.
#   4. Secondary: disagreement endpoint, strata, share closer.
#
# OUTPUTS
#   Data/gtex_wbc_rrbs/srb_pair_results.rds
#   Results/Task59_Primary.csv, Task59_ByStratum.csv, Fig43_ExternalTest.png
# ================================================================

suppressPackageStartupMessages({ library(data.table); library(jsonlite); library(ggplot2) })
pick_dir <- function(...) { for (d in c(...)) if (dir.exists(d)) return(d)
  stop("no candidate directory exists") }
repo_dir <- pick_dir("/mmfs1/home/wln26/Experiments.Outlier.July31.2026",
                     path.expand("~/Desktop/bioinformatics-research"))
res <- file.path(repo_dir, "Results"); dat <- file.path(repo_dir, "Data", "gtex_wbc_rrbs")

# ----------------------------------------------------------------
# 0. the lock is intact
# ----------------------------------------------------------------
lock <- fromJSON(file.path(res, "Task58_LockedProtocol.json"))
for (f in names(lock$sha256))
  stopifnot(digest::digest(file = file.path(repo_dir, f), algo = "sha256") == lock$sha256[[f]])
source(file.path(repo_dir, "Scripts", "SharedReadNoise.Core.R"))
source(file.path(repo_dir, "Scripts", "SharedReadNoise.Benchmark.R"))
P <- lock$parameters
AUTO <- paste0("chr", 1:22); s6 <- function(x) signif(x, 4); r4 <- function(x) round(x, 4)
cat("lock verified\n")

cpg <- readRDS(file.path(repo_dir, "Data", "hg19_seq", "cpg_positions_hg19.rds"))
I0 <- c(0L, cumsum(lengths(cpg))[-22])
kept <- fread(file.path(res, "Task56_ExternalSamples.csv"))[kept == TRUE]
wbc <- as.data.table(fromJSON(file.path(dat, "wbc_samples.json")))
wbc <- wbc[gsm %in% kept$gsm][order(gsm)]
files <- file.path(dat, basename(sapply(wbc$files, `[`, 2)))
D <- length(files); stopifnot(D == lock$donors, D >= 30, all(file.exists(files)))
read_pat <- function(f) {
  p <- fread(cmd = sprintf("gzip -dc '%s'", f), header = FALSE, sep = "\t", col.names = c("chr", "idx", "pat", "n"),
             colClasses = c("character", "integer", "character", "integer"))
  split(p[chr %in% AUTO, .(idx, pat, n)], p$chr[p$chr %in% AUTO])
}

# ----------------------------------------------------------------
# 1. pair universe from counts only
# ----------------------------------------------------------------
lookup <- lapply(cpg, function(p) { ok <- which(diff(p) <= P$max_gap_bp); v <- rep(NA_integer_, length(p)); v[ok] <- seq_along(ok); v })
deep <- lapply(lookup, function(v) integer(sum(!is.na(v))))
t0 <- Sys.time()
for (s in seq_len(D)) {
  ps <- read_pat(files[s])
  for (k in 1:22) {
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
cat("pair universe:", nrow(uni), "pairs (", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min)\n")

# ----------------------------------------------------------------
# 2. split, count and accumulate
# ----------------------------------------------------------------
acc <- lapply(1:22, function(k) srb_accumulator(uni[chrn == k, .N]))
uidx <- lapply(1:22, function(k) uni[chrn == k, idx])
sub <- uni[idx %% P$bootstrap_subsample_every == 0L, .(chrn, idx, pos)]
SM <- c("E", "xA", "yA", "cA", "biB", "bjB", "biC", "bjC")
subm <- setNames(lapply(SM, function(v) matrix(0, nrow(sub), D)), SM); sub_key <- paste(sub$chrn, sub$idx)
for (s in seq_len(D)) {
  ps <- read_pat(files[s])
  for (k in 1:22) {
    pk <- ps[[AUTO[k]]]; if (is.null(pk) || !nrow(pk)) next
    m <- srb_donor_chr(pk, I0[k], length(cpg[[k]]), cpg[[k]], seed = P$split_seed + 1000L * s + k, max_gap = P$max_gap_bp)
    if (is.null(m)) next
    rows <- match(m$idx, uidx[[k]]); keep <- !is.na(rows); if (!any(keep)) next
    m <- m[keep]; srb_update(acc[[k]], m, rows[keep])
    at <- match(paste(k, m$idx), sub_key); h <- which(!is.na(at))
    if (length(h)) { subm$E[at[h], s] <- 1; for (v in SM[-1]) subm[[v]][at[h], s] <- m[[v]][h] }
  }
  rm(ps); invisible(gc())
  cat(sprintf("donor %2d/%d done (%.1f min)\n", s, D, as.numeric(Sys.time() - t0, units = "mins")))
}
est <- rbindlist(lapply(1:22, function(k) {
  a <- copy(acc[[k]])[, `:=`(chrn = k, idx = uidx[[k]])]; keep <- a$n >= P$min_donors
  cbind(a[keep, .(chrn, idx)], srb_finish(a[keep], min_donors = P$min_donors))
}))
est <- uni[est, on = c("chrn", "idx")]
saveRDS(list(est = est, sub = sub, subm = subm, built = format(Sys.time())), file.path(dat, "srb_pair_results.rds"))

# ----------------------------------------------------------------
# 3. primary endpoint and its interval
# ----------------------------------------------------------------
ok <- rowSums(subm$E) >= P$min_donors
E <- subm$E[ok, ]; sb <- sub[ok]
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
  d <- dC_w(w); k <- !is.na(d) & bw > 0
  c(sum(bw[k] * d[k]) / sum(bw[k]), sum(is.na(d)))
}, numeric(2))
sub_est <- mean(dC_w(rep(1, D)), na.rm = TRUE); bias <- mean(bt[1, ]) - sub_est   # locked v2: bias-corrected percentile
prim <- data.table(donors = D, pairs = nrow(est), delta_mse_cov = s6(mean(est$dC)),
                   subsample_pairs = nrow(sb), subsample_estimate = s6(sub_est),
                   lo95 = s6(quantile(bt[1, ], 0.025) - bias), hi95 = s6(quantile(bt[1, ], 0.975) - bias),
                   plain_percentile_lo95 = s6(quantile(bt[1, ], 0.025)), plain_percentile_hi95 = s6(quantile(bt[1, ], 0.975)),
                   relative_change = r4(mean(est$dC) / mean(est$loss_raw)), share_cor_closer = r4(mean(est$loss_cor < est$loss_raw)),
                   delta_mse_d2 = s6(mean(est$dD)), mean_failed_pairs_per_draw = r4(mean(bt[2, ])))
prim[, verdict := fcase(hi95 < 0, "primary criterion met: correction closer to the independent reference",
                        lo95 > 0, "correction has higher error", default = "improvement not established")]
fwrite(prim, file.path(res, "Task59_Primary.csv"))
cat("\nPRIMARY (locked):\n"); print(prim)

# ----------------------------------------------------------------
# 4. secondary: strata
# ----------------------------------------------------------------
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
fwrite(strata, file.path(res, "Task59_ByStratum.csv"))
cat("\nSecondary, by stratum:\n"); print(strata[, .(stratum, level, pairs, delta_mse_cov, relative_change, share_cor_closer)])

SURF <- "#fcfcfb"; INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"; GRID <- "#e1e0d9"; BLUE <- "#2a78d6"
f43 <- strata[stratum != "all pairs"][, `:=`(stratum = factor(stratum, levels = unique(stratum)), level = factor(level, levels = unique(level)))]
p43 <- ggplot(f43, aes(relative_change, level)) +
  geom_vline(xintercept = 0, colour = MUTED, linewidth = 0.4) + geom_point(colour = BLUE, size = 3) +
  facet_grid(stratum ~ ., scales = "free_y", space = "free_y") +
  labs(title = "External test (locked): is the corrected covariance closer to an independent reference?",
       subtitle = paste0(D, " white-blood-cell RRBS samples (GSE233417). Primary: ", prim$verdict, ".\n",
                         "Delta MSE ", prim$delta_mse_cov, ", 95% interval ", prim$lo95, " to ", prim$hi95, "."),
       x = "Relative change in squared error, corrected vs observed (left of 0 = correction helps)", y = NULL,
       caption = "Source: Results/Task59_Primary.csv and Task59_ByStratum.csv (Task 59)") +
  theme_minimal(base_size = 10) +
  theme(plot.background = element_rect(fill = SURF, colour = NA), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = GRID, linewidth = 0.3),
        strip.text.y = element_text(angle = 0, hjust = 0, face = "bold", colour = INK),
        plot.title = element_text(colour = INK, face = "bold", size = 12), plot.title.position = "plot",
        plot.subtitle = element_text(colour = INK2, size = 9), plot.caption = element_text(colour = MUTED, hjust = 0),
        plot.caption.position = "plot", axis.text = element_text(colour = INK2))
ggsave(file.path(res, "Fig43_ExternalTest.png"), p43, width = 10, height = 6, dpi = 200, bg = SURF)
cat("\ndone\n")
