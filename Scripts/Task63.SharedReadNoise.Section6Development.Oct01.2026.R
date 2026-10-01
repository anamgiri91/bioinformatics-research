# ================================================================
# Task 63 - The plan's section 6 comparisons, developed on GTEx colon
#           (development data), before they are locked for the
#           second external cohort.
# Date: Oct 1, 2026
#
# WHAT IS COMPARED (plan_shared_read_noise.md section 6)
#   1. Correlation recovery: Pearson, variance-only and full correction,
#      against a reference correlation built from B and C. Failures are
#      counted, never dropped silently.
#   2. Pair ranking: seven scores from A only. Each is judged by the mean
#      B/C reference squared disagreement (Qref) of its top 1, 5, 10 and
#      25% of pairs; lower means it found pairs that truly agree. This is
#      a secondary utility comparison: the scores do not estimate the
#      same quantity.
#   3. An existing correlated-error method: Ding & Gentleman (2003),
#      Bioconductor MeasurementError.cor. Maximum likelihood under
#      bivariate normality with per-donor standard errors and one error
#      correlation, here with se = sqrt(vhat) per donor.
#
# THE SPLIT
#   Task 55's seeds (20261002 + 1000 x donor + chromosome) and pair
#   universe, so A, B and C are the same fragments as in Task 55. The
#   script stops unless its covariance endpoint equals Task 55's for
#   every pair.
#
# SETTLED HERE, FOR THE LOCK
#   - xi: Chatterjee's xi, the larger of its two directions; ties in x
#     broken at random (seeded); undefined if either profile is constant.
#   - Undefined scores rank last. Pairs tied at the cut get fractional
#     weights.
#   - A correlation estimate is valid when its variances are positive and
#     it lies in [-1, 1]. The reference is defined on the same terms.
#   - The baseline runs on 2,000 random subsample pairs with a defined
#     reference (seed 20261204), because it fits one optimisation per
#     pair.
#
# OUTPUTS
#   Results/Task63_CorrelationRecovery.csv, Task63_Ranking.csv,
#   Task63_MEcorBaseline.csv, Fig45_Section6Development.png
#   Data/gtex_colon_rrbs/s6_pair_results.rds (cache)
# ================================================================

suppressPackageStartupMessages({ library(data.table); library(jsonlite); library(ggplot2); library(MeasurementError.cor) })
pick_dir <- function(...) { for (d in c(...)) if (dir.exists(d)) return(d)
  stop("no candidate directory exists") }
repo_dir <- pick_dir("/mmfs1/home/wln26/Experiments.Outlier.July31.2026",
                     path.expand("~/Desktop/bioinformatics-research"))
res <- file.path(repo_dir, "Results"); dat <- file.path(repo_dir, "Data", "gtex_colon_rrbs")
for (f in c("SharedReadNoise.Core.R", "SharedReadNoise.Benchmark.R", "Task47.FrozenScore.Sep29.2026.R",
            "SharedReadNoise.Section6.R")) source(file.path(repo_dir, "Scripts", f))
AUTO <- paste0("chr", 1:22); SEED_SPLIT <- 20261002L; SUB_EVERY <- 25L; SEED_XI <- 20261205L; SEED_ME <- 20261204L
s4 <- function(x) signif(x, 4)

# ----------------------------------------------------------------
# 0. checks of the new code
# ----------------------------------------------------------------
xi_one <- function(x, y) {                     # textbook formula, x without ties
  n <- length(x); ys <- y[order(x)]; r <- sapply(ys, function(v) sum(y <= v)); l <- sapply(ys, function(v) sum(y >= v))
  1 - n * sum(abs(diff(r))) / (2 * sum(l * (n - l)))
}
set.seed(7); X <- matrix(runif(40 * 25), 40); Y <- matrix(round(runif(40 * 25) * 6) / 6, 40)   # ties in y only
stopifnot(isTRUE(all.equal(s6_xi_rows(X, Y, 1L), sapply(1:40, function(i) xi_one(X[i, ], Y[i, ])))))
stopifnot(isTRUE(all.equal(s6_xi_rows(matrix(1:30, 1), matrix((1:30)^2, 1), 1L), 1 - 3 / 31)))
stopifnot(isTRUE(all.equal(s6_top_mean(c(5, 4, 4, 4, 1), c(1, 2, 3, 4, 5), 0.4), (1 + (2 + 3 + 4) / 3) / 2)),
          isTRUE(all.equal(s6_top_mean(c(NA, 3, NA, 2), c(9, 1, 7, 2), 0.75), (1 + 2 + (9 + 7) / 2) / 3)))
set.seed(8); a <- runif(30); b <- a + rnorm(30, 0, 0.1); bi <- a + rnorm(30, 0, 0.05); ci <- a + rnorm(30, 0, 0.05)
bj <- b + rnorm(30, 0, 0.05); cj <- b + rnorm(30, 0, 0.05)
ac <- s6_accumulator(1L)
for (k in 1:30) s6_update(ac, data.table(xA = a[k], yA = b[k], vxA = 0.001, vyA = 0.002, cA = 0.0005, biB = bi[k], bjB = bj[k],
                                         biC = ci[k], bjC = cj[k], NiA = 5, NjA = 5, rA = 5), 1L)
f <- s6_finish(ac, min_donors = 2L)
stopifnot(isTRUE(all.equal(f$r_pearson, cor(a, b))), isTRUE(all.equal(f$ref_var_i, cov(bi, ci))),
          isTRUE(all.equal(f$r_vo, cov(a, b) / sqrt((var(a) - 0.001) * (var(b) - 0.002)))),
          isTRUE(all.equal(f$r_full, (cov(a, b) - 0.0005) / sqrt((var(a) - 0.001) * (var(b) - 0.002)))),
          isTRUE(all.equal(f$rho_ref, ((cov(bi, cj) + cov(ci, bj)) / 2) / sqrt(cov(bi, ci) * cov(bj, cj)))))
cat("checks passed: xi matches the textbook formula; top-fraction weights; correlation algebra\n")

# ----------------------------------------------------------------
# 1. split (Task 55 seeds), count and accumulate
# ----------------------------------------------------------------
cpg <- readRDS(file.path(repo_dir, "Data", "hg19_seq", "cpg_positions_hg19.rds")); I0 <- c(0L, cumsum(lengths(cpg))[-22])
samples <- fromJSON(file.path(dat, "colon_samples.json")); D <- length(samples$gsm); stopifnot(D == 29)
uni <- readRDS(file.path(dat, "srn_pair_estimates.rds"))$est[, .(chrn, idx, pos, gap)]; setkey(uni, chrn, idx)
out_file <- file.path(dat, "s6_pair_results.rds")
SM <- c("E", "xA", "yA", "vxA", "vyA", "cA", "biB", "bjB", "biC", "bjC")
if (!file.exists(out_file)) {
  acc <- lapply(1:22, function(k) s6_accumulator(uni[chrn == k, .N])); uidx <- lapply(1:22, function(k) uni[chrn == k, idx])
  sub <- uni[idx %% SUB_EVERY == 0L, .(chrn, idx, pos, gap)]; sub_key <- paste(sub$chrn, sub$idx)
  subm <- setNames(lapply(SM, function(v) matrix(0, nrow(sub), D)), SM)
  t0 <- Sys.time()
  for (s in seq_len(D)) {
    p <- fread(cmd = sprintf("gzip -dc '%s'", file.path(dat, basename(samples$files[[s]][2]))), header = FALSE, sep = "\t",
               col.names = c("chr", "idx", "pat", "n"), colClasses = c("character", "integer", "character", "integer"))
    ps <- split(p[chr %in% AUTO, .(idx, pat, n)], p$chr[p$chr %in% AUTO]); rm(p)
    for (k in 1:22) {
      pk <- ps[[AUTO[k]]]; if (is.null(pk) || !nrow(pk)) next
      m <- srb_donor_chr(pk, I0[k], length(cpg[[k]]), cpg[[k]], seed = SEED_SPLIT + 1000L * s + k)
      if (is.null(m)) next
      rows <- match(m$idx, uidx[[k]]); keep <- !is.na(rows); if (!any(keep)) next
      m <- m[keep]; s6_update(acc[[k]], m, rows[keep])
      at <- match(paste(k, m$idx), sub_key); h <- which(!is.na(at))
      if (length(h)) { subm$E[at[h], s] <- 1; for (v in SM[-1]) subm[[v]][at[h], s] <- m[[v]][h] }
    }
    rm(ps); invisible(gc())
    cat(sprintf("donor %2d/%d done (%.1f min)\n", s, D, as.numeric(Sys.time() - t0, units = "mins")))
  }
  est <- rbindlist(lapply(1:22, function(k) { a <- copy(acc[[k]])[, `:=`(chrn = k, idx = uidx[[k]])]; keep <- a$n >= 20L
    cbind(a[keep, .(chrn, idx)], s6_finish(a[keep], min_donors = 20L)) }))
  saveRDS(list(est = est, sub = sub, subm = subm, built = format(Sys.time())), out_file)
}
z <- readRDS(out_file); est <- uni[z$est, on = c("chrn", "idx")]; sub <- z$sub; subm <- z$subm
old <- readRDS(file.path(dat, "srb_pair_results.rds"))$est
chk <- old[est, on = c("chrn", "idx"), nomatch = NULL]
stopifnot(nrow(chk) == nrow(est), nrow(est) == nrow(old), isTRUE(all.equal(chk$dC, chk$i.dC)), isTRUE(all.equal(chk$dD, chk$i.dD)))
cat("same split as Task 55: covariance and disagreement endpoints equal for all", nrow(est), "pairs\n")

# ----------------------------------------------------------------
# 2. correlation recovery
# ----------------------------------------------------------------
corr <- s6_corr_table(est)
fwrite(corr, file.path(res, "Task63_CorrelationRecovery.csv"))
cat("\nCorrelation recovery (reference from B and C):\n"); print(corr[, lapply(.SD, function(x) if (is.numeric(x)) s4(x) else x)])

# ----------------------------------------------------------------
# 3. pair ranking on the subsample
# ----------------------------------------------------------------
sp <- est[sub, on = c("chrn", "idx"), nomatch = NA]
ok <- !is.na(sp$n)                                   # subsample pairs that reached 20 donors
subo <- lapply(subm, function(M) M[ok, , drop = FALSE])
chk2 <- s6_d2_se(subo)
stopifnot(isTRUE(all.equal(chk2$d2_mean, sp$cor_d2_A[ok])))      # per-donor terms rebuild the summed estimate
scores <- s6_scores(subo, sp$gap[ok], sp$raw_d2_A[ok], sp$cor_d2_A[ok], SEED_XI)
rk <- s6_rank_table(scores, sp$ref_d2[ok], varies = sp$var_x_A[ok] > 0 & sp$var_y_A[ok] > 0)
fwrite(rk, file.path(res, "Task63_Ranking.csv"))
cat("\nPair ranking: mean reference squared disagreement of the top pairs (lower is better):\n")
print(dcast(rk, set + label ~ top_fraction, value.var = "mean_qref_top", fun.aggregate = function(x) s4(x)))

# ----------------------------------------------------------------
# 4. existing correlated-error method (Ding & Gentleman 2003)
# ----------------------------------------------------------------
set.seed(SEED_ME)
cand <- which(ok & !is.na(sp$rho_ref)); pick <- sort(sample(cand, min(2000L, length(cand))))
t1 <- Sys.time()
me <- rbindlist(lapply(pick, function(i) {
  e <- subm$E[i, ] == 1
  fit <- tryCatch(cor.me.vector(subm$xA[i, e], sqrt(subm$vxA[i, e]), subm$yA[i, e], sqrt(subm$vyA[i, e])), error = function(x) NULL)
  data.table(row = i, me_corr_true = if (is.null(fit)) NA_real_ else unname(fit$estimate["corr.true"]),
             me_corr_error = if (is.null(fit)) NA_real_ else unname(fit$estimate["corr.me"]),
             converged = !is.null(fit) && fit$convergence == 0)
}))
me[, `:=`(rho_ref = sp$rho_ref[row], r_pearson = sp$r_pearson[row], r_vo = sp$r_vo[row], r_full = sp$r_full[row])]
me[converged == FALSE, me_corr_true := NA_real_]
v <- me[!is.na(me_corr_true) & !is.na(r_pearson) & !is.na(r_vo) & !is.na(r_full)]
mtab <- data.table(estimator = c("Pearson", "variance-only", "full correction", "Ding & Gentleman (MeasurementError.cor)"),
                   pairs_tried = nrow(me), valid_share = c(mean(!is.na(me$r_pearson)), mean(!is.na(me$r_vo)), mean(!is.na(me$r_full)),
                                                           mean(!is.na(me$me_corr_true))),
                   common_pairs = nrow(v),
                   mse_common = c(mean((v$r_pearson - v$rho_ref)^2), mean((v$r_vo - v$rho_ref)^2), mean((v$r_full - v$rho_ref)^2),
                                  mean((v$me_corr_true - v$rho_ref)^2)),
                   mean_error_common = c(mean(v$r_pearson - v$rho_ref), mean(v$r_vo - v$rho_ref), mean(v$r_full - v$rho_ref),
                                         mean(v$me_corr_true - v$rho_ref)),
                   minutes = round(as.numeric(Sys.time() - t1, units = "mins"), 1))
fwrite(mtab, file.path(res, "Task63_MEcorBaseline.csv"))
cat("\nExisting correlated-error method, on", nrow(me), "pairs:\n"); print(mtab[, lapply(.SD, function(x) if (is.numeric(x)) s4(x) else x)])

# ----------------------------------------------------------------
# 5. figure: the ranking comparison
# ----------------------------------------------------------------
SURF <- "#fcfcfb"; INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"; GRID <- "#e1e0d9"
pal <- c("#2a78d6", "#5b93dd", "#8cb2e6", "#bcd2f0", "#d03b3b", "#e39a2d", "#5c9e5a", "#8a63b8", "#3f3f3f", "#898781")
f45 <- copy(rk)[, label := factor(label, levels = unname(S6_SCORES))]
p45 <- ggplot(f45, aes(top_fraction * 100, mean_qref_top, colour = label)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) + facet_wrap(~ set, scales = "free_y") +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1e-7, base = 10), breaks = c(0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2),
                     labels = c("0", "1e-6", "1e-5", "1e-4", "1e-3", "1e-2")) +
  scale_colour_manual(values = setNames(pal, unname(S6_SCORES)), name = NULL) +
  scale_x_continuous(breaks = S6_FRACTIONS * 100) +
  labs(title = "Development (GTEx colon): which score from A finds pairs that truly agree?",
       subtitle = paste0("Mean reference squared disagreement (B with C) of each score's top pairs; lower is better. ",
                         format(sum(ok), big.mark = ","), " subsample pairs."),
       x = "Top % of pairs selected", y = "Mean reference squared disagreement",
       caption = "Source: Results/Task63_Ranking.csv (Task 63)") +
  theme_minimal(base_size = 10) +
  theme(plot.background = element_rect(fill = SURF, colour = NA), panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = GRID, linewidth = 0.3), legend.position = "right",
        plot.title = element_text(colour = INK, face = "bold", size = 12), plot.title.position = "plot",
        plot.subtitle = element_text(colour = INK2, size = 9), plot.caption = element_text(colour = MUTED, hjust = 0),
        plot.caption.position = "plot", axis.text = element_text(colour = INK2))
ggsave(file.path(res, "Fig45_Section6Development.png"), p45, width = 10, height = 5, dpi = 200, bg = SURF)
cat("\ndone\n")
