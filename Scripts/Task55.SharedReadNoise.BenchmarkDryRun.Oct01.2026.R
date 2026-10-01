# ================================================================
# Task 55 - Three-way fragment benchmark: DRY RUN on GTEx colon
# Date: Oct 1, 2026
#
# DEVELOPMENT ONLY (plan_shared_read_noise.md sections 9-10). GTEx colon
# has already been examined, so this tests the pipeline and shows what
# the endpoint looks like. It is not the plan's external validation.
#
# DESIGN (as in the plan)
#   Each donor's fragments are split at random into A, B and C (each
#   fragment, with all its calls, goes to one part with probability
#   1/3; seed 20261002 + 1000 x donor + chromosome). Estimates use A
#   only: observed covariance s_A and corrected Chat_A. The reference R
#   is the B-by-C cross covariance, which shares no reads.
#   A donor counts for a pair if both sites have at least 2 reads in A,
#   B and C, and A has 0 or at least 2 shared fragments. A pair needs 20
#   such donors.
#   Endpoint: mean over pairs of (Chat_A - R)^2 - (s_A - R)^2. Negative
#   means the corrected estimate is closer to the reference.
#   Disagreement: the same with corrected and observed squared
#   disagreement against Qref.
#   Uncertainty (not yet calibrated): 2,000 joint bootstrap draws of
#   donors and 1-Mb blocks within chromosomes, on a fixed 1-in-25 pair
#   subsample (CpG index divisible by 25), seed 20261003.
#
# Pair universe: the 4,447,998 pairs of Task 53 (20+ donors overall).
# A pair that fails there cannot pass here, because here every site
# needs at least 6 reads.
#
# OUTPUTS
#   Data/gtex_colon_rrbs/srb_pair_results.rds
#   Results/Task55_Endpoint.csv, Task55_ByStratum.csv,
#   Task55_Bootstrap.csv, Fig41_BenchmarkDryRun.png
# ================================================================

suppressPackageStartupMessages({ library(data.table); library(jsonlite); library(ggplot2) })

pick_dir <- function(...) { for (d in c(...)) if (dir.exists(d)) return(d)
  stop("no candidate directory exists") }
repo_dir <- pick_dir("/mmfs1/home/wln26/Experiments.Outlier.July31.2026",
                     path.expand("~/Desktop/bioinformatics-research"))
res <- file.path(repo_dir, "Results"); dat <- file.path(repo_dir, "Data", "gtex_colon_rrbs")
source(file.path(repo_dir, "Scripts", "SharedReadNoise.Core.R"))
source(file.path(repo_dir, "Scripts", "SharedReadNoise.Benchmark.R"))
AUTO <- paste0("chr", 1:22); SEED_SPLIT <- 20261002L; SEED_BOOT <- 20261003L; SUB_EVERY <- 25L; NBOOT <- 2000L
s6 <- function(x) signif(x, 4); r4 <- function(x) round(x, 4)

cpg <- readRDS(file.path(repo_dir, "Data", "hg19_seq", "cpg_positions_hg19.rds"))
I0 <- c(0L, cumsum(lengths(cpg))[-22])
samples <- fromJSON(file.path(dat, "colon_samples.json"))
sid <- sub("^GSM[0-9]+_(GTEX-[A-Z0-9]+)-.*", "\\1", basename(sapply(samples$files, `[`, 1)))
D <- length(sid); stopifnot(D == 29)
uni <- readRDS(file.path(dat, "srn_pair_estimates.rds"))$est[, .(chrn, idx, pos, gap, raw_var_i, raw_var_j)]
setkey(uni, chrn, idx)

# ----------------------------------------------------------------
# 0. checks on one donor, chr22
# ----------------------------------------------------------------
p1 <- fread(cmd = sprintf("gzip -dc '%s'", file.path(dat, basename(samples$files[[1]][2]))), header = FALSE, sep = "\t",
            col.names = c("chr", "idx", "pat", "n"), colClasses = c("character", "integer", "character", "integer"))
pk <- p1[chr == "chr22", .(idx, pat, n)]
sp <- srb_split_pat(pk, 1L)
whole <- attr(srn_pat_counts(pk, I0[22], length(cpg[[22]]), cpg[[22]]), "site_counts")
parts <- rbindlist(lapply(sp, function(q) attr(srn_pat_counts(q, I0[22], length(cpg[[22]]), cpg[[22]]), "site_counts")))
parts <- parts[, .(N = sum(N), M = sum(M)), keyby = idx]
stopifnot(identical(whole$idx, parts$idx), all(whole$N == parts$N), all(whole$M == parts$M),
          sum(sp$A$n) + sum(sp$B$n) + sum(sp$C$n) == sum(pk$n))
set.seed(1); a <- runif(30); b <- runif(30); bb <- runif(30); cc <- runif(30)
ac <- srb_accumulator(1L)
for (k in 1:30) srb_update(ac, data.table(xA = a[k], yA = b[k], vxA = 0, vyA = 0, cA = 0.01, biB = a[k], bjB = bb[k],
                                          biC = cc[k], bjC = b[k], NiA = 5, NjA = 5, rA = 5), 1L)
f <- srb_finish(ac, min_donors = 2L)
stopifnot(isTRUE(all.equal(f$raw_cov_A, cov(a, b))), isTRUE(all.equal(f$ref_cov, (cov(a, b) + cov(cc, bb)) / 2)),
          isTRUE(all.equal(f$dC, 0.01^2 - 2 * 0.01 * (f$raw_cov_A - f$ref_cov))))
rm(p1, pk, sp, whole, parts); invisible(gc())
cat("checks passed: the split keeps every count; cross covariance and endpoint algebra match\n")

# ----------------------------------------------------------------
# 1. split, count and accumulate, donor by donor
# ----------------------------------------------------------------
out_file <- file.path(dat, "srb_pair_results.rds")
if (!file.exists(out_file)) {
  acc <- lapply(1:22, function(k) srb_accumulator(uni[chrn == k, .N]))
  uidx <- lapply(1:22, function(k) uni[chrn == k, idx])
  sub <- uni[idx %% SUB_EVERY == 0L, .(chrn, idx, pos)]
  P <- nrow(sub); SM <- c("E", "xA", "yA", "cA", "biB", "bjB", "biC", "bjC")
  subm <- setNames(lapply(SM, function(v) matrix(0, P, D)), SM)
  sub_key <- paste(sub$chrn, sub$idx)
  t0 <- Sys.time()
  for (s in seq_len(D)) {
    p <- fread(cmd = sprintf("gzip -dc '%s'", file.path(dat, basename(samples$files[[s]][2]))), header = FALSE, sep = "\t",
               col.names = c("chr", "idx", "pat", "n"), colClasses = c("character", "integer", "character", "integer"))
    ps <- split(p[chr %in% AUTO, .(idx, pat, n)], p$chr[p$chr %in% AUTO]); rm(p)
    for (k in 1:22) {
      pk <- ps[[AUTO[k]]]; if (is.null(pk) || !nrow(pk)) next
      m <- srb_donor_chr(pk, I0[k], length(cpg[[k]]), cpg[[k]], seed = SEED_SPLIT + 1000L * s + k)
      if (is.null(m)) next
      rows <- match(m$idx, uidx[[k]]); keep <- !is.na(rows)
      if (!any(keep)) next
      m <- m[keep]; srb_update(acc[[k]], m, rows[keep])
      at <- match(paste(k, m$idx), sub_key); h <- which(!is.na(at))
      if (length(h)) {
        subm$E[at[h], s] <- 1
        for (v in SM[-1]) subm[[v]][at[h], s] <- m[[v]][h]
      }
    }
    rm(ps); invisible(gc())
    cat(sprintf("donor %2d/%d %s done (%.1f min)\n", s, D, sid[s], as.numeric(Sys.time() - t0, units = "mins")))
  }
  est <- rbindlist(lapply(1:22, function(k) {
    a <- copy(acc[[k]])[, `:=`(chrn = k, idx = uidx[[k]])]
    keep <- a$n >= 20L
    cbind(a[keep, .(chrn, idx)], srb_finish(a[keep], min_donors = 20L))
  }))
  saveRDS(list(est = est, sub = sub, subm = subm, built = format(Sys.time())), out_file)
}
z <- readRDS(out_file); est <- z$est
est <- uni[est, on = c("chrn", "idx")]
cat("\npairs with 20+ donors eligible in A, B and C:", nrow(est), "\n")

# ----------------------------------------------------------------
# 2. the endpoint, overall and by stratum
# ----------------------------------------------------------------
endpoint <- function(d) d[, .(pairs = .N, delta_mse_cov = s6(mean(dC)), mean_loss_raw = s6(mean(loss_raw)),
                              mean_loss_cor = s6(mean(loss_cor)), relative_change = r4(mean(dC) / mean(loss_raw)),
                              share_cor_closer = r4(mean(loss_cor < loss_raw)),
                              delta_mse_d2 = s6(mean(dD)), mean_noise_A = s6(mean(mean_noise_A)),
                              mean_raw_cov_A = s6(mean(raw_cov_A)), mean_ref_cov = s6(mean(ref_cov)),
                              mean_cor_cov_A = s6(mean(cor_cov_A)))]
est[, spread := (sqrt(pmax(raw_var_i, 0)) + sqrt(pmax(raw_var_j, 0))) / 2]
by_one <- function(name, lev) { d <- copy(est)[, level := lev]
  d[, endpoint(.SD), by = level][order(level)][, `:=`(stratum = name, level = as.character(level))] }
strata <- rbindlist(list(
  endpoint(est)[, `:=`(stratum = "all pairs", level = "all")],
  by_one("distance (bp)", cut(est$gap, c(0, 10, 20, 40, 60, 100, 150, 200))),
  by_one("mean depth in A (reads)", cut(est$depth_A, c(2, 5, 10, 20, Inf), right = FALSE)),
  by_one("overlap in A", cut(est$overlap_A, c(0, 0.25, 0.5, 0.75, 1.0001), right = FALSE)),
  by_one("spread (mean raw SD)", cut(est$spread, c(0, 0.02, 0.05, 0.1, Inf), right = FALSE))), use.names = TRUE)
setcolorder(strata, c("stratum", "level"))
fwrite(strata[stratum == "all pairs"], file.path(res, "Task55_Endpoint.csv"))
fwrite(strata, file.path(res, "Task55_ByStratum.csv"))
cat("\nEndpoint (negative = corrected estimate closer to the B/C reference):\n")
print(strata[, .(stratum, level, pairs, delta_mse_cov, relative_change, share_cor_closer, delta_mse_d2)], nrows = 40)

# ----------------------------------------------------------------
# 3. joint donor and block bootstrap on the 1-in-25 subsample
# ----------------------------------------------------------------
S <- z$subm; sub <- z$sub
ok <- rowSums(S$E) >= 20                      # the subsample's eligible pairs (fixed mask)
E <- S$E[ok, ]; sub <- sub[ok]
pre <- list(EX = E * S$xA[ok, ], EY = E * S$yA[ok, ], EXY = E * S$xA[ok, ] * S$yA[ok, ], EC = E * S$cA[ok, ],
            EBI = E * S$biB[ok, ], ECJ = E * S$bjC[ok, ], EBIC = E * S$biB[ok, ] * S$bjC[ok, ],
            ECI = E * S$biC[ok, ], EBJ = E * S$bjB[ok, ], ECIB = E * S$biC[ok, ] * S$bjB[ok, ])
dC_w <- function(w) {
  n <- drop(E %*% w); g <- function(M) drop(M %*% w)
  cv <- function(sxy, sx, sy) (sxy - sx * sy / n) / (n - 1)
  W <- cv(g(pre$EXY), g(pre$EX), g(pre$EY)); U <- W - g(pre$EC) / n
  R <- (cv(g(pre$EBIC), g(pre$EBI), g(pre$ECJ)) + cv(g(pre$ECIB), g(pre$ECI), g(pre$EBJ))) / 2
  d <- (U - R)^2 - (W - R)^2; d[n < 2] <- NA_real_; d
}
blk <- paste(sub$chrn, (sub$pos - 1L) %/% 1e6L)
blocks <- unique(data.table(chrn = sub$chrn, blk))
point_sub <- mean(dC_w(rep(1, D)))
set.seed(SEED_BOOT)
bt <- vapply(seq_len(NBOOT), function(b) {
  w <- tabulate(sample.int(D, D, replace = TRUE), D)
  bm <- blocks[, .(blk = sample(blk, .N, replace = TRUE)), by = chrn][, .N, by = blk]
  bw <- bm$N[match(blk, bm$blk)]; bw[is.na(bw)] <- 0
  d <- dC_w(w); k <- !is.na(d) & bw > 0
  c(sum(bw[k] * d[k]) / sum(bw[k]), sum(is.na(d)))
}, numeric(2))
bias <- mean(bt[1, ]) - point_sub          # Task 57: the plain percentile interval is shifted upward
boot <- data.table(subsample_pairs = nrow(sub), point_estimate_subsample = s6(point_sub),
                   point_estimate_all_pairs = s6(mean(est$dC)),
                   bias_corrected_lo95 = s6(quantile(bt[1, ], 0.025) - bias), bias_corrected_hi95 = s6(quantile(bt[1, ], 0.975) - bias),
                   plain_percentile_lo95 = s6(quantile(bt[1, ], 0.025)), plain_percentile_hi95 = s6(quantile(bt[1, ], 0.975)),
                   draws = NBOOT, mean_failed_pairs_per_draw = r4(mean(bt[2, ])))
fwrite(boot, file.path(res, "Task55_Bootstrap.csv"))
cat("\nJoint donor and block bootstrap (not yet calibrated):\n"); print(boot)

# ----------------------------------------------------------------
# 4. figure
# ----------------------------------------------------------------
SURF <- "#fcfcfb"; INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"; GRID <- "#e1e0d9"; BLUE <- "#2a78d6"
f41 <- strata[stratum != "all pairs"][, `:=`(stratum = factor(stratum, levels = unique(stratum)),
                                             level = factor(level, levels = unique(level)))]
p41 <- ggplot(f41, aes(relative_change, level)) +
  geom_vline(xintercept = 0, colour = MUTED, linewidth = 0.4) +
  geom_point(colour = BLUE, size = 3) +
  facet_grid(stratum ~ ., scales = "free_y", space = "free_y") +
  labs(title = "Dry run on GTEx colon: is the corrected covariance closer to an independent reference?",
       subtitle = paste0("Relative change in squared error against the B/C reference (left of 0 = correction helps). ",
                         "Overall: ", sprintf("%+.1f%%", 100 * strata[stratum == "all pairs", relative_change]),
                         " over ", format(nrow(est), big.mark = ","), " pairs. Development data, not validation."),
       x = "Relative change in squared error, corrected vs observed", y = NULL,
       caption = "Source: Results/Task55_ByStratum.csv (Task 55)") +
  theme_minimal(base_size = 10) +
  theme(plot.background = element_rect(fill = SURF, colour = NA), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = GRID, linewidth = 0.3),
        strip.text.y = element_text(angle = 0, hjust = 0, face = "bold", colour = INK),
        plot.title = element_text(colour = INK, face = "bold", size = 12), plot.title.position = "plot",
        plot.subtitle = element_text(colour = INK2, size = 9), plot.caption = element_text(colour = MUTED, hjust = 0),
        plot.caption.position = "plot", axis.text = element_text(colour = INK2))
ggsave(file.path(res, "Fig41_BenchmarkDryRun.png"), p41, width = 10, height = 6.5, dpi = 200, bg = SURF)
cat("\ndone\n")
