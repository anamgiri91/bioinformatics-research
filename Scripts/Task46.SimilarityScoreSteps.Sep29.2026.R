# ================================================================
# Task 46 - A better similarity score for neighbouring CG sites,
#           built step by step, with a fair test at every step
# Date: Sep 29, 2026
#
# THE REQUEST (2026-09-29)
#   Improve the Manhattan agreement, use a better co-methylation
#   measure, add genomic (base-pair) distance, and add a factor for
#   known limitations such as repetitive sequence. Build it step by
#   step and check whether each step is an improvement.
#
# PAIRS
#   The 32,664 truly consecutive RRBS pairs: two adjacent rows of the
#   chr22 file, both observed in all 18 samples (as in Task 41).
#
# A FAIR TEST (the Codex test only scored agreement, which Manhattan
# wins by construction)
#   Leave one patient out, 18 times. Every score is computed on the
#   other 17 patients only. Two outcomes are measured on the held-out
#   patient:
#     A  agreement:  |x_k - y_k|, the held-out patient's two values
#                    (lower is better). Rewards "same level".
#     B  prediction gain: how much the neighbour improves the guess of
#                    the held-out value, beyond the site's own mean.
#                    For site x: guess = mean_x + slope * (y_k - mean_y)
#                    from a straight line fitted on the 17 patients,
#                    clipped to [0, 1]; gain = |x_k - mean_x| -
#                    |x_k - guess|. Averaged over both directions
#                    (higher is better). Rewards "moves together".
#   For each score and fold: the mean outcome in the top 10% of pairs
#   by score (unscorable pairs ranked last, ties split evenly), and
#   the Spearman correlation of score with outcome, both over the
#   pairs the score can rate and over the pairs every score can rate
#   (both sites vary in the 17 patients).
#
# THE STEPS (the chosen variant of each step, marked *, was fixed
# before any result was seen; every variant is still reported)
#   Step 0  baselines: Manhattan S (raw beta), ordinal S (10 levels),
#           Pearson, Spearman, and the Codex score
#           F = S [w + (1 - w) E], E = dCor excess over its exact
#           permutation mean, w = exp(-gap / 100 bp).
#   Step 1  agreement:
#           * S_asin     compare asin(sqrt(beta)); this evens out the
#                        noise of read ratios, which is largest near 0.5
#             S_tol      ignore differences below one level (0.1)
#             S_asin_tol both
#   Step 2  co-methylation:
#             D_pos      max(0, Spearman); 0 when a site is constant
#           * D_shrunk   D_pos * m / (m + 3), where m is the smaller of
#                        the two sites' counts of samples more than 0.1
#                        from the site median (little spread, little trust)
#   Step 3  combine and distance:
#             C_fixed    S_asin * [w + (1 - w) D_shrunk], w = exp(-gap/100)
#             C_fitted   the same, with the length scale fitted in each
#                        fold from how Spearman falls with distance in
#                        the 17 training patients
#           * C_rule     C_fitted, but a pair where both sites are
#                        constant scores S_asin (no co-methylation is
#                        possible there, so none is required)
#   Step 4  reliability:
#           * C_rel      C_rule * R_i * R_j, R = Bismap multi-read
#                        mappability for 50-bp reads after bisulfite
#                        conversion (Karimzadeh et al. 2018), 0 to 1
#             C_norep    C_rule, with pairs touching a RepeatMasker
#                        repeat set to 0
#           plus a check of whether pairs in repeats or with low
#           mappability do worse on the two outcomes, within state.
#
# ADDED AFTER THE FIRST RESULTS (so not chosen in advance)
#   The combined scores sit at different points of a trade-off:
#   more weight on agreement gives lower error but less gain. One
#   length scale per score cannot say which formula is better. So
#   each combined formula is traced over many length scales (a
#   curve), and the curves are compared at the same agreement error.
#   The curves are: Codex F; the new C; the new C with raw-beta S;
#   the new C with Codex's dCor excess in place of D_shrunk; and
#   versions of these with one constant weight for every pair (no
#   distance). The first curves showed the distance weight lowering
#   gain in the new C (it raised it in Codex F), so two formulas were
#   then built from what worked:
#     "dCor, no distance":  S_asin * [c + (1 - c) E], c constant
#     "distance as a prior": S_asin * [c + (1 - c) D_prior], where
#        D_prior = (m * D_pos + 3 * r(gap)) / (m + 3) pulls a thin
#        Spearman toward r(gap), the Spearman typical at that gap in
#        the training patients, instead of toward 0.
#   The constructed cases then showed two weak spots of Codex F:
#   dCor ignores the sign (an inverse pair gets full credit) and the
#   scale (a shared 0.49-0.51 wiggle gets full credit). So one more:
#     "F, guarded": S_raw * [w + (1 - w) E_guard], with
#        E_guard = E * [Spearman > 0] * m / (m + 3), the same
#        low-spread shrinkage as D_shrunk.
#     "F, sign guard only": the same with E * [Spearman > 0] alone.
#   Because these were chosen after looking, every curve is also
#   run on each half of chr22 separately: a real difference should
#   show in both halves.
#   The repeat and mappability check is also repeated at matched
#   spread (20 bins of the pair's mean standard deviation), since
#   spread drives both outcomes.
#
# INPUTS
#   Data/NN.hg38.18P.forw.chr22.w.header.txt
#   Results/Task39_RRBS_CompleteSiteStates.csv.gz   (states, for strata)
#   Data/annotation/*   from Task46.FetchAnnotation.Sep29.2026.R
#
# OUTPUTS (Results/)
#   Task46_Steps_Summary.csv       mean over folds, every score
#   Task46_Steps_Folds.csv         every fold, every score
#   Task46_Steps_Comparisons.csv   each step against the one before
#   Task46_Reliability_Strata.csv  outcomes by repeat / mappability
#   Task46_Reliability_SpreadMatched.csv  the same at matched spread
#   Task46_DistanceDecay.csv       median Spearman by gap (all 18 patients)
#   Task46_Frontier.csv            each combined formula over length scales
#   Task46_Frontier_Matched.csv    gain at the same agreement error
#   Task46_Examples.csv            constructed cases through the scores
#   Task46_PairScores.csv.gz       all scores on all 18 patients, per pair
#   Fig33_SimilaritySteps.png      the steps on the two outcomes
#   Fig34_SimilarityFrontier.png   the combined formulas as curves
#
# LIMITS
#   Exploratory, on the one cohort that was already examined. Folds
#   share 16 of 17 patients and neighbouring pairs share sites, so
#   there are no p-values. No read counts, so no read-depth weight.
# ================================================================

suppressPackageStartupMessages({ library(data.table); library(matrixStats); library(ggplot2) })

pick_dir <- function(...) { for (d in c(...)) if (dir.exists(d)) return(d)
  stop("no candidate directory exists") }
repo_dir <- pick_dir("/mmfs1/home/wln26/Experiments.Outlier.July31.2026",
                     path.expand("~/Desktop/bioinformatics-research"))
res <- file.path(repo_dir, "Results")
ann_dir <- file.path(repo_dir, "Data", "annotation")
r4 <- function(x) round(x, 4)
DELTA <- 0.1; KAPPA <- 3; TOP <- 0.10; ELL0 <- 100

# ----------------------------------------------------------------
# 0. pairs
# ----------------------------------------------------------------
d <- fread(file.path(repo_dir, "Data", "NN.hg38.18P.forw.chr22.w.header.txt"), na.strings = "NA",
           showProgress = FALSE)
samp <- grep("^P[0-9]+NF$", names(d), value = TRUE)
m <- as.matrix(d[, ..samp]); pos_all <- d$pos; rm(d)
comp <- rowSums(!is.na(m)) == 18
i1 <- which(comp[-length(comp)] & comp[-1]); i2 <- i1 + 1L
X <- m[i1, ]; Y <- m[i2, ]; rm(m); invisible(gc())
p1 <- pos_all[i1]; p2 <- pos_all[i2]; gap <- p2 - p1
stopifnot(length(i1) == 32664, !anyNA(X), !anyNA(Y))
cs <- fread(file.path(res, "Task39_RRBS_CompleteSiteStates.csv.gz"), select = c("pos", "state"))
st1 <- cs$state[match(p1, cs$pos)]; st2 <- cs$state[match(p2, cs$pos)]
stopifnot(!anyNA(st1), !anyNA(st2))
cat("pairs:", length(i1), " gap range:", range(gap), "bp\n")

# ----------------------------------------------------------------
# 1. per-site reliability annotation (cached)
# ----------------------------------------------------------------
ann_file <- file.path(ann_dir, "rrbs_chr22_pairsite_annotation.tsv.gz")
sites <- sort(unique(c(p1, p2)))
if (!file.exists(ann_file)) {
  at_value <- function(f) {                   # value of a 1-based interval track at each site
    bw <- fread(file.path(ann_dir, f), select = c("start", "end", "score"))
    setorder(bw, start)
    j <- findInterval(sites, bw$start); jj <- pmax(j, 1L)
    v <- fifelse(j > 0 & sites <= bw$end[jj], bw$score[jj], 0)
    rm(bw); invisible(gc()); v
  }
  map50 <- at_value("bismap_k50_multi_chr22.tsv.gz")
  map100 <- at_value("bismap_k100_multi_chr22.tsv.gz")
  un <- fread(file.path(ann_dir, "bismap_k50_C2T_unique_chr22.tsv.gz")); setorder(un, start)
  j <- findInterval(sites, un$start)
  uniq50 <- j > 0 & sites <= un$end[pmax(j, 1L)]
  rmk <- fread(file.path(ann_dir, "rmsk_chr22.tsv.gz"))
  rmk[, `:=`(s1 = start + 1L, e1 = end)]      # UCSC tables are 0-based, half open
  setkey(rmk, s1, e1)
  q <- data.table(pos = sites, s1 = sites, e1 = sites); setkey(q, s1, e1)
  ov <- foverlaps(q, rmk[, .(s1, e1, repClass)], type = "within", nomatch = NA)
  ov <- ov[, .(rep_class = if (all(is.na(repClass))) NA_character_
                           else paste(sort(unique(na.omit(repClass))), collapse = "+")), by = pos]
  ann <- data.table(pos = sites, map50, map100, uniq50)[ov, on = "pos"]
  ann[, in_repeat := !is.na(rep_class)]
  fwrite(ann, ann_file)
}
ann <- fread(ann_file)
stopifnot(identical(ann$pos, sites))
a1 <- match(p1, ann$pos); a2 <- match(p2, ann$pos)
R1 <- ann$map50[a1]; R2 <- ann$map50[a2]
rep_pair <- ann$in_repeat[a1] | ann$in_repeat[a2]
cls1 <- ann$rep_class[a1]; cls2 <- ann$rep_class[a2]
rep_class_pair <- fifelse(!is.na(cls1), cls1, fifelse(!is.na(cls2), cls2, "no repeat"))
cat("sites in a repeat:", sum(ann$in_repeat), "of", nrow(ann),
    " | Bismap k50 mappability: share = 1:", r4(mean(ann$map50 == 1)),
    " share < 0.5:", r4(mean(ann$map50 < 0.5)), " share = 0:", r4(mean(ann$map50 == 0)), "\n")
cat("pairs touching a repeat:", sum(rep_pair), " | pairs with a site below mappability 1:",
    sum(pmin(R1, R2) < 1), "\n")

# ----------------------------------------------------------------
# 2. the scores
# ----------------------------------------------------------------
lev <- function(v) pmin(floor(round(v * 10, 8)) + 1, 10)
rowcor <- function(a, b) {
  a <- a - rowMeans(a); b <- b - rowMeans(b)
  den <- sqrt(rowSums(a^2) * rowSums(b^2))
  ifelse(den > 0, rowSums(a * b) / den, NA_real_)
}
dcor_excess <- function(Xt, Yt) {             # Codex's E: dCor^2 excess over its permutation mean
  n <- ncol(Xt)
  vapply(seq_len(nrow(Xt)), function(p) {
    x <- Xt[p, ]; y <- Yt[p, ]
    if (max(x) == min(x) || max(y) == min(y)) return(0)
    a <- abs(outer(x, x, "-")); b <- abs(outer(y, y, "-"))
    A <- a - rowMeans(a) - rep(colMeans(a), each = n) + mean(a)
    B <- b - rowMeans(b) - rep(colMeans(b), each = n) + mean(b)
    den <- sqrt(sum(A * A) * sum(B * B)); if (den <= 0) return(0)
    q <- sum(A * B) / den
    q0 <- sum(diag(A)) * sum(diag(B)) / ((n - 1) * den)
    if (1 - q0 <= 1e-12) return(0)
    max(0, (q - q0) / (1 - q0))
  }, 0)
}
fit_decay <- function(spear, g) {             # typical Spearman = a * exp(-gap / ell), from gap bins
  b <- cut(g, c(0, 10, 20, 50, 100, Inf))
  t <- data.table(b, g, spear)[!is.na(spear), .(g = as.numeric(median(g)), r = median(spear), n = .N), by = b]
  t <- t[r > 0 & n >= 20]
  if (nrow(t) < 3) return(list(a = 0, ell = Inf))
  cf <- coef(lm(log(r) ~ g, data = t, weights = n))
  if (cf[["g"]] >= 0) return(list(a = weighted.mean(t$r, t$n), ell = Inf))
  list(a = exp(cf[[1]]), ell = -1 / cf[["g"]])
}
score_all <- function(Xt, Yt, dec, g = gap, r1 = R1, r2 = R2, rp = rep_pair) {
  ell <- dec$ell
  cx <- rowMaxs(Xt) == rowMins(Xt); cy <- rowMaxs(Yt) == rowMins(Yt)
  S_raw  <- 1 - rowMeans(abs(Xt - Yt))
  S_ord  <- 1 - rowMeans(abs(lev(Xt) - lev(Yt))) / 9
  da     <- abs(asin(sqrt(Xt)) - asin(sqrt(Yt)))
  S_asin <- 1 - rowMeans(da) / (pi / 2)
  S_tol  <- 1 - rowMeans(pmax(abs(Xt - Yt) - DELTA, 0)) / (1 - DELTA)   # matrix first keeps dim
  S_asin_tol <- 1 - rowMeans(pmax(da - DELTA, 0)) / (pi / 2 - DELTA)
  pear  <- rowcor(Xt, Yt)
  spear <- rowcor(rowRanks(Xt, ties.method = "average"), rowRanks(Yt, ties.method = "average"))
  mx <- rowSums(abs(Xt - rowMedians(Xt)) > DELTA); my <- rowSums(abs(Yt - rowMedians(Yt)) > DELTA)
  msup <- pmin(mx, my)
  D_pos <- pmax(spear, 0); D_pos[is.na(D_pos)] <- 0
  D_shrunk <- D_pos * msup / (msup + KAPPA)
  r0 <- pmin(1, dec$a * exp(-g / ell))        # post hoc: shrink toward the Spearman expected at this gap
  D_prior <- (msup * D_pos + KAPPA * r0) / (msup + KAPPA)
  E <- dcor_excess(Xt, Yt)
  E_sign <- E * fifelse(!is.na(spear) & spear > 0, 1, 0)                            # post hoc
  E_guard <- E_sign * msup / (msup + KAPPA)                                          # post hoc
  w0 <- exp(-g / ELL0); wf <- exp(-g / ell)
  C_fixed  <- S_asin * (w0 + (1 - w0) * D_shrunk)
  C_fitted <- S_asin * (wf + (1 - wf) * D_shrunk)
  C_rule   <- fifelse(cx & cy, S_asin, C_fitted)
  data.table(S_raw, S_ord, Pearson = pear, Spearman = spear,
             F_codex = S_raw * (w0 + (1 - w0) * E), F_guard = S_raw * (w0 + (1 - w0) * E_guard),
             S_asin, S_tol, S_asin_tol, D_pos, D_shrunk, E_dcor = E, D_prior,
             C_fixed, C_fitted, C_rule, C_rel = C_rule * r1 * r2, C_norep = C_rule * !rp,
             both_constant = cx & cy, both_vary = !cx & !cy, support = msup, E_guard, E_sign)
}
SCORES <- c("S_raw", "S_ord", "Pearson", "Spearman", "F_codex", "S_asin", "S_tol", "S_asin_tol",
            "D_pos", "D_shrunk", "E_dcor", "D_prior", "C_fixed", "C_fitted", "C_rule", "F_guard", "C_rel",
            "C_norep")
STEP <- c(S_raw = 0, S_ord = 0, Pearson = 0, Spearman = 0, F_codex = 0, S_asin = 1, S_tol = 1,
          S_asin_tol = 1, D_pos = 2, D_shrunk = 2, E_dcor = 2, D_prior = 2, C_fixed = 3, C_fitted = 3,
          C_rule = 3, F_guard = 3, C_rel = 4, C_norep = 4)
CHOSEN <- c("S_asin", "D_shrunk", "C_rule", "C_rel")
POSTHOC <- c("D_prior", "F_guard")                       # built after the first results

# ----------------------------------------------------------------
# 3. leave one patient out
# ----------------------------------------------------------------
outcomes <- function(Xt, Yt, xk, yk) {
  mx <- rowMeans(Xt); my <- rowMeans(Yt)
  vx <- rowMeans((Xt - mx)^2); vy <- rowMeans((Yt - my)^2); cxy <- rowMeans((Xt - mx) * (Yt - my))
  bxy <- ifelse(vy > 0, cxy / vy, 0); byx <- ifelse(vx > 0, cxy / vx, 0)
  xhat <- pmin(1, pmax(0, mx + bxy * (yk - my))); yhat <- pmin(1, pmax(0, my + byx * (xk - mx)))
  list(A = abs(xk - yk),
       B = ((abs(xk - mx) - abs(xk - xhat)) + (abs(yk - my) - abs(yk - yhat))) / 2)
}
top_mean <- function(s, y, K) {               # mean y over the K highest s; ties split evenly
  t <- -sort(-s, partial = K)[K]
  above <- s > t; tied <- s == t
  (sum(y[above]) + (K - sum(above)) * mean(y[tied])) / K
}
evaluate <- function(sc, out) {
  K <- ceiling(TOP * nrow(sc)); cm <- sc$both_vary
  rbindlist(lapply(SCORES, function(nm) {
    s <- sc[[nm]]; def <- !is.na(s); s2 <- fifelse(def, s, -Inf)
    stopifnot(all(def[cm]))
    data.table(score = nm, defined = mean(def),
               top_A = top_mean(s2, out$A, K), top_B = top_mean(s2, out$B, K),
               rho_A = cor(s[def], -out$A[def], method = "spearman"),
               rho_B = cor(s[def], out$B[def], method = "spearman"),
               rho_A_common = cor(s[cm], -out$A[cm], method = "spearman"),
               rho_B_common = cor(s[cm], out$B[cm], method = "spearman"))
  }))
}
ELL_GRID <- c(0, 3, 5, 10, 20, 30, 50, 75, 100, 150, 200, 300, 500, 1000, 3000, Inf)
MID <- median(p1)
SUBSETS <- list(all = rep(TRUE, length(p1)), first_half = p1 < MID, second_half = p1 >= MID)
frontier <- function(sc, out, g = gap) {      # each combined formula over the length-scale grid
  cc <- sc$both_constant
  rbindlist(lapply(ELL_GRID, function(l) {
    w <- exp(-g / l); cbar <- mean(w)         # l = 0 gives w = 0 (gaps are >= 2); l = Inf gives 1
    rule <- function(S, v) fifelse(cc, S, S * v)
    fam <- list(                              # distance weight w, or one weight cbar for every pair
      `Codex F`                    = sc$S_raw * (w + (1 - w) * sc$E_dcor),
      `Codex F, no distance`       = sc$S_raw * (cbar + (1 - cbar) * sc$E_dcor),
      `F, guarded`                 = sc$S_raw * (w + (1 - w) * sc$E_guard),
      `F, sign guard only`         = sc$S_raw * (w + (1 - w) * sc$E_sign),
      `New C`                      = rule(sc$S_asin, w + (1 - w) * sc$D_shrunk),
      `New C, no distance`         = rule(sc$S_asin, cbar + (1 - cbar) * sc$D_shrunk),
      `New C, raw-beta S`          = rule(sc$S_raw, w + (1 - w) * sc$D_shrunk),
      `New C, dCor excess`         = rule(sc$S_asin, w + (1 - w) * sc$E_dcor),
      `New C, dCor, no distance`   = rule(sc$S_asin, cbar + (1 - cbar) * sc$E_dcor),
      `New C, distance as a prior` = rule(sc$S_asin, cbar + (1 - cbar) * sc$D_prior))
    rbindlist(lapply(names(SUBSETS), function(sn) {
      ix <- SUBSETS[[sn]]; K <- ceiling(TOP * sum(ix))
      rbindlist(lapply(names(fam), function(f) data.table(
        subset = sn, family = f, ell = l, mean_w = cbar,
        top_A = top_mean(fam[[f]][ix], out$A[ix], K), top_B = top_mean(fam[[f]][ix], out$B[ix], K))))
    }))
  }))
}
t0 <- Sys.time()
fold_out <- lapply(seq_len(18), function(k) {
  Xt <- X[, -k]; Yt <- Y[, -k]
  sp <- rowcor(rowRanks(Xt, ties.method = "average"), rowRanks(Yt, ties.method = "average"))
  dec <- fit_decay(sp, gap)
  sc <- score_all(Xt, Yt, dec)
  out <- outcomes(Xt, Yt, X[, k], Y[, k])
  ev <- evaluate(sc, out)[, `:=`(fold = k, ell = dec$ell, decay_a = dec$a)]
  list(ev = ev, fr = frontier(sc, out)[, fold := k])
})
folds <- rbindlist(lapply(fold_out, `[[`, "ev")); fr_folds <- rbindlist(lapply(fold_out, `[[`, "fr"))
cat("18 folds in", round(as.numeric(Sys.time() - t0, units = "secs")), "s; fitted length scale:",
    paste(round(range(folds$ell)), collapse = " to "), "bp\n")
fwrite(folds, file.path(res, "Task46_Steps_Folds.csv"))

summ <- folds[, .(step = STEP[score[1]], chosen = score[1] %in% CHOSEN, defined = r4(mean(defined)),
                  top10_agreement_error = r4(mean(top_A)), top10_prediction_gain = r4(mean(top_B)),
                  rho_agreement = r4(mean(rho_A)), rho_gain = r4(mean(rho_B)),
                  rho_agreement_common = r4(mean(rho_A_common)), rho_gain_common = r4(mean(rho_B_common))),
              by = score]
fwrite(summ, file.path(res, "Task46_Steps_Summary.csv"))
cat("\nMean over 18 folds (agreement error: lower is better; prediction gain: higher is better):\n")
print(summ)

cmp_pairs <- list(
  c("S_asin", "S_raw"), c("S_tol", "S_raw"), c("S_asin_tol", "S_raw"),
  c("D_pos", "Spearman"), c("D_shrunk", "D_pos"), c("D_shrunk", "Spearman"), c("D_shrunk", "Pearson"),
  c("D_shrunk", "E_dcor"), c("D_prior", "D_shrunk"), c("C_fixed", "S_asin"), c("C_fixed", "D_shrunk"), c("C_fitted", "C_fixed"), c("C_rule", "C_fitted"),
  c("C_rel", "C_rule"), c("C_norep", "C_rule"),
  c("C_rule", "F_codex"), c("F_guard", "F_codex"), c("C_rel", "F_codex"), c("C_rel", "S_raw"), c("C_rel", "Spearman"))
cmp <- rbindlist(lapply(cmp_pairs, function(pr) {
  a <- folds[score == pr[1]][order(fold)]; b <- folds[score == pr[2]][order(fold)]
  data.table(candidate = pr[1], against = pr[2],
             folds_lower_error = sum(a$top_A < b$top_A), folds_higher_gain = sum(a$top_B > b$top_B),
             change_error = round(mean(a$top_A - b$top_A), 5), change_gain = round(mean(a$top_B - b$top_B), 5),
             change_rho_error_common = round(mean(a$rho_A_common - b$rho_A_common), 4),
             change_rho_gain_common = round(mean(a$rho_B_common - b$rho_B_common), 4))
}))
fwrite(cmp, file.path(res, "Task46_Steps_Comparisons.csv"))
cat("\nStep by step (folds out of 18 where the candidate is better; change = candidate minus the other):\n")
print(cmp)

# ----------------------------------------------------------------
# 3b. the combined formulas as curves, compared at the same error
#     (added after the first results; see the header)
# ----------------------------------------------------------------
fr <- fr_folds[, .(mean_w = r4(mean_w[1]), top10_agreement_error = round(mean(top_A), 5),
                   top10_prediction_gain = round(mean(top_B), 5),
                   error_lo = round(min(top_A), 5), error_hi = round(max(top_A), 5),
                   gain_lo = round(min(top_B), 5), gain_hi = round(max(top_B), 5)), by = .(subset, family, ell)]
fwrite(fr, file.path(res, "Task46_Frontier.csv"))
cat("\nEach combined formula over length scales, all pairs (mean over 18 folds), gain then error:\n")
print(dcast(fr[subset == "all"], ell + mean_w ~ family, value.var = "top10_prediction_gain"))
print(dcast(fr[subset == "all"], ell + mean_w ~ family, value.var = "top10_agreement_error"))

at_error <- function(A, B, target) {          # gain at a given error, read off one curve
  if (target < min(A) || target > max(A)) return(NA_real_)
  o <- order(A); approx(A[o], B[o], xout = target, ties = mean)$y
}
TARGETS <- c(0.016, 0.018, 0.020, 0.025, 0.030, 0.040)
mt <- fr_folds[, .(target = TARGETS, gain = vapply(TARGETS, function(t) at_error(top_A, top_B, t), 0)),
               by = .(subset, fold, family)]
mw <- dcast(mt, subset + fold + target ~ family, value.var = "gain")
QUESTIONS <- list(                            # question, formula a, formula b (does a beat b?)
  c("1 distance weight vs one weight, in Codex F",     "Codex F", "Codex F, no distance"),
  c("2 distance weight vs one weight, in New C",       "New C", "New C, no distance"),
  c("3 distance weight vs one weight, in C with dCor", "New C, dCor excess", "New C, dCor, no distance"),
  c("4 distance as a prior vs no distance, in New C",  "New C, distance as a prior", "New C, no distance"),
  c("5 asin(sqrt) vs raw beta, in New C",              "New C", "New C, raw-beta S"),
  c("6 asin(sqrt) vs raw beta, in Codex F",            "New C, dCor excess", "Codex F"),
  c("7 shrunk Spearman vs dCor, in New C",             "New C", "New C, dCor excess"),
  c("8 planned New C vs Codex F",                      "New C", "Codex F"),
  c("9 F, guarded vs Codex F",                         "F, guarded", "Codex F"),
  c("10 F, sign guard only vs Codex F",                "F, sign guard only", "Codex F"))
matched <- rbindlist(lapply(QUESTIONS, function(q) {
  mw[, .(question = q[1], a = q[2], b = q[3], folds = sum(!is.na(get(q[2])) & !is.na(get(q[3]))),
         folds_a_higher_gain = sum(get(q[2]) > get(q[3]), na.rm = TRUE),
         gain_a = r4(mean(get(q[2]), na.rm = TRUE)), gain_b = r4(mean(get(q[3]), na.rm = TRUE))),
     by = .(subset, agreement_error = target)]
}))
fwrite(matched, file.path(res, "Task46_Frontier_Matched.csv"))
cat("\nAt the same agreement error, how often formula a gives more gain than b",
    "(fold-by-error cells, summed over the six error levels):\n")
print(dcast(matched[, .(cell = paste0(sum(folds_a_higher_gain), "/", sum(folds))), by = .(question, subset)],
            question ~ subset, value.var = "cell"))
cat("\nAll pairs, by error level:\n"); print(matched[subset == "all", !"subset"])

# ----------------------------------------------------------------
# 4. do repeats and low mappability do worse, within state?
# ----------------------------------------------------------------
full_sp <- rowcor(rowRanks(X, ties.method = "average"), rowRanks(Y, ties.method = "average"))
dec_full <- fit_decay(full_sp, gap)
full <- score_all(X, Y, dec_full)
loo <- rbindlist(lapply(seq_len(18), function(k) {
  o <- outcomes(X[, -k], Y[, -k], X[, k], Y[, k]); data.table(pair = seq_along(o$A), A = o$A, B = o$B)
}))[, .(A = mean(A), B = mean(B)), by = pair][order(pair)]
stopifnot(identical(loo$pair, seq_along(i1)))
grp <- data.table(A = loo$A, B = loo$B, both_constant = full$both_constant,
                  state_pair = fifelse(st1 == st2 & st1 %in% c("L", "H"), paste0(st1, "-", st2), "other"),
                  in_repeat = rep_pair, low_map = pmin(R1, R2) < 1,
                  not_unique = !(ann$uniq50[a1] & ann$uniq50[a2]), rep_class = rep_class_pair)
strat <- function(col) {
  g <- copy(grp)[, group := get(col)]
  rbind(g[, .(state_pair = "all", pairs = .N, agreement_error = r4(mean(A)), prediction_gain = r4(mean(B)),
              both_constant = r4(mean(both_constant))), by = group],
        g[, .(pairs = .N, agreement_error = r4(mean(A)), prediction_gain = r4(mean(B)),
              both_constant = r4(mean(both_constant))), by = .(state_pair, group)])[
    , split := col][order(state_pair, group)]
}
strata <- rbindlist(list(strat("in_repeat"), strat("low_map"), strat("not_unique")), use.names = TRUE)
setcolorder(strata, c("split", "state_pair", "group"))
strata[, group := as.character(group)]
by_class <- grp[, .(split = "rep_class", state_pair = "all", pairs = .N, agreement_error = r4(mean(A)),
                    prediction_gain = r4(mean(B)), both_constant = r4(mean(both_constant))),
                by = .(group = rep_class)][pairs >= 100][order(-pairs)]
strata <- rbind(strata, by_class, use.names = TRUE)
fwrite(strata, file.path(res, "Task46_Reliability_Strata.csv"))
cat("\nHeld-out outcomes by repeat and mappability (mean over the 18 held-out patients):\n"); print(strata)

# the same at matched spread: compare flagged pairs with unflagged pairs of the same spread
spread <- (rowSds(X) + rowSds(Y)) / 2
NB <- 20L
sbin <- pmin(NB, 1L + floor(NB * (frank(spread, ties.method = "first") - 1) / length(spread)))
flags <- list(in_repeat = rep_pair, SINE = grepl("SINE", rep_class_pair),
              simple_or_low_complexity = grepl("Simple_repeat|Low_complexity", rep_class_pair),
              low_map = pmin(R1, R2) < 1, not_unique = !(ann$uniq50[a1] & ann$uniq50[a2]))
spread_matched <- rbindlist(lapply(names(flags), function(nm) {
  f <- flags[[nm]]
  t <- data.table(sbin, f, A = loo$A, B = loo$B)[, .(A = mean(A), B = mean(B), n = .N), by = .(sbin, f)]
  w <- merge(t[f == TRUE], t[f == FALSE], by = "sbin", suffixes = c("_f", "_o"))
  data.table(flag = nm, flagged_pairs = sum(f), matched_pairs = sum(w$n_f),
             error_diff_raw = r4(mean(loo$A[f]) - mean(loo$A[!f])),
             error_diff_matched = r4(sum(w$n_f * (w$A_f - w$A_o)) / sum(w$n_f)),
             gain_diff_raw = r4(mean(loo$B[f]) - mean(loo$B[!f])),
             gain_diff_matched = r4(sum(w$n_f * (w$B_f - w$B_o)) / sum(w$n_f)),
             bins_flagged_worse_error = sum(w$A_f > w$A_o), bins_flagged_lower_gain = sum(w$B_f < w$B_o),
             bins = nrow(w))
}))
fwrite(spread_matched, file.path(res, "Task46_Reliability_SpreadMatched.csv"))
cat("\nFlagged minus unflagged pairs, raw and at matched spread (positive error = worse agreement):\n")
print(spread_matched)

# how co-methylation falls with distance, on all 18 patients
decay <- data.table(bin = cut(gap, c(0, 10, 20, 50, 100, Inf)), gap, sp = full_sp)[
  !is.na(sp), .(pairs = .N, median_gap = as.numeric(median(gap)), median_spearman = r4(median(sp)),
                share_spearman_above_0.5 = r4(mean(sp > 0.5))), by = bin][order(bin)]
fwrite(decay, file.path(res, "Task46_DistanceDecay.csv"))
cat("\nSpearman by gap (the fitted length scale comes from this):\n"); print(decay)

pairs_out <- cbind(data.table(pos1 = p1, pos2 = p2, gap, state1 = st1, state2 = st2,
                              map50_1 = R1, map50_2 = R2, repeat_pair = rep_pair),
                   full[, lapply(.SD, r4), .SDcols = SCORES], full[, .(both_constant, support)])
fwrite(pairs_out, file.path(res, "Task46_PairScores.csv.gz"))
cat("\nfitted on all 18 patients: typical Spearman =", r4(dec_full$a), "x exp(-gap /",
    round(dec_full$ell), "bp)\n")

# ----------------------------------------------------------------
# 5. constructed cases
# ----------------------------------------------------------------
base <- seq(0.1, 0.9, length.out = 18)
cases <- list(
  identical_constant   = list(rep(0.05, 18), rep(0.05, 18)),
  different_constants  = list(rep(0.05, 18), rep(0.95, 18)),
  identical_varying    = list(base, base),
  offset_same_order    = list(base * 0.5, base * 0.5 + 0.3),
  inverse              = list(base, 1 - base),
  one_shared_rare      = list(c(0.95, rep(0.05, 17)), c(0.95, rep(0.05, 17))),
  tiny_shared_wiggle   = list(seq(0.49, 0.51, length.out = 18), seq(0.49, 0.51, length.out = 18)))
ex <- rbindlist(lapply(names(cases), function(nm) {
  z <- cases[[nm]]
  s <- score_all(matrix(z[[1]], 1), matrix(z[[2]], 1), dec_full, g = 50, r1 = 1, r2 = 1, rp = FALSE)
  cbind(data.table(case = nm), s[, .(S_raw, S_asin, Spearman, D_shrunk, E_dcor, D_prior, F_codex, F_guard, C_rule)][, lapply(.SD, r4)])
}))
fwrite(ex, file.path(res, "Task46_Examples.csv"))
cat("\nConstructed cases, 50 bp apart, full mappability:\n"); print(ex)

# ----------------------------------------------------------------
# 6. figure: every score on the two outcomes, mean and range over folds
# ----------------------------------------------------------------
SURF <- "#fcfcfb"; INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"; GRID <- "#e1e0d9"
BLUE <- "#2a78d6"; GRAY <- "#c3c2b7"
NICE <- c(S_raw = "Manhattan S, raw beta", S_ord = "Manhattan S, 10 levels", Pearson = "Pearson",
          Spearman = "Spearman", F_codex = "Codex F", S_asin = "S on asin(sqrt) scale",
          S_tol = "S, ignore < 0.1", S_asin_tol = "S, asin(sqrt) + ignore < 0.1",
          D_pos = "Spearman, negatives to 0", D_shrunk = "Spearman, shrunk for low spread",
          E_dcor = "dCor excess (Codex E)", D_prior = "Spearman, shrunk to distance prior",
          C_fixed = "Combined, 100 bp scale", C_fitted = "Combined, fitted scale",
          C_rule = "Combined, fitted + constant rule", F_guard = "Codex F, guarded", C_rel = "Combined x mappability",
          C_norep = "Combined, repeats set to 0")
STEP_LAB <- c("0 Baselines", "1 Agreement", "2 Co-methylation", "3 Combined with distance", "4 Reliability")
MET <- c(A = "Held-out agreement, 1 - |x - y|", B = "Held-out prediction gain")
fl <- rbind(folds[, .(score, fold, metric = MET[["A"]], v = 1 - top_A)],
            folds[, .(score, fold, metric = MET[["B"]], v = top_B)])
fs <- fl[, .(mean = mean(v), lo = min(v), hi = max(v)), by = .(score, metric)]
fs[, `:=`(name = factor(NICE[score], levels = rev(NICE)),
          step = factor(STEP_LAB[STEP[score] + 1], levels = STEP_LAB),
          metric = factor(metric, levels = MET),
          kind = fcase(score %in% CHOSEN, "Chosen before the run", score %in% POSTHOC, "Added after the first run",
                       default = "Other variant or baseline"))]
ref <- fs[score == "S_raw", .(metric, mean)]
p33 <- ggplot(fs, aes(y = name)) +
  geom_vline(data = ref, aes(xintercept = mean), colour = MUTED, linetype = "22", linewidth = 0.4) +
  geom_segment(aes(x = lo, xend = hi, yend = name), colour = GRAY, linewidth = 0.8, lineend = "round") +
  geom_point(aes(x = mean, fill = kind), shape = 21, colour = SURF, stroke = 0.6, size = 3.2) +
  scale_fill_manual(values = c(`Chosen before the run` = BLUE, `Added after the first run` = "#eb6834",
                               `Other variant or baseline` = MUTED),
                    breaks = c("Chosen before the run", "Added after the first run", "Other variant or baseline")) +
  facet_grid(step ~ metric, scales = "free", space = "free_y", switch = "y") +
  labs(title = "Building the similarity score step by step",
       subtitle = paste0("The top 10% of the 32,664 truly consecutive RRBS pairs by each score, judged on a ",
                         "patient left out of the scoring.\nDot: mean over 18 leave-one-patient-out folds. ",
                         "Line: lowest to highest fold. Right is better in both panels."),
       x = NULL, y = NULL,
       caption = "Dashed line: Manhattan S on raw beta. Source: Results/Task46_Steps_Folds.csv (Task 46)") +
  theme_minimal(base_size = 10.5) +
  theme(plot.background = element_rect(fill = SURF, colour = NA),
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(colour = GRID, linewidth = 0.3),
        panel.spacing.x = unit(14, "pt"), panel.spacing.y = unit(6, "pt"),
        strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1, colour = INK, face = "bold"),
        strip.text.x = element_text(colour = INK, face = "bold", hjust = 0),
        plot.title = element_text(colour = INK, face = "bold", size = 13), plot.title.position = "plot",
        plot.subtitle = element_text(colour = INK2, size = 9.5), plot.caption = element_text(colour = MUTED, hjust = 0),
        plot.caption.position = "plot",
        legend.position = "top", legend.justification = "left", legend.title = element_blank(),
        axis.text = element_text(colour = INK2))
ggsave(file.path(res, "Fig33_SimilaritySteps.png"), p33, width = 11, height = 7.2, dpi = 200, bg = SURF)

# ----------------------------------------------------------------
# 7. figure: the combined formulas as curves
# ----------------------------------------------------------------
# one panel per formula: distance weight (solid) against one weight for every pair (dashed);
# colour = formula, the same in every panel; light gray = the other formulas, for comparison
FAM <- c("Codex F", "Codex F, no distance", "F, guarded", "New C", "New C, no distance", "New C, distance as a prior",
         "New C, dCor excess", "New C, dCor, no distance")
FAM_COL <- c("#52514e", "#52514e", "#eb6834", BLUE, BLUE, "#4a3aa7", "#1baf7a", "#1baf7a")
FAM_LT <- c("solid", "22", "solid", "solid", "22", "solid", "solid", "22")
PANELS <- list(`Codex F (raw-beta S, dCor excess)` = c("Codex F", "Codex F, no distance", "F, guarded"),
               `New C as planned (asin S, shrunk Spearman)` =
                 c("New C", "New C, no distance", "New C, distance as a prior"),
               `New C with dCor excess (asin S, dCor)` = c("New C, dCor excess", "New C, dCor, no distance"))
SOLID <- c("Codex F", "New C", "New C, dCor excess")
f34 <- rbindlist(lapply(names(PANELS), function(p) fr[subset == "all" & family %in% PANELS[[p]]][, panel := p]))
ctx <- rbindlist(lapply(names(PANELS), function(p)
  fr[subset == "all" & family %in% setdiff(SOLID, PANELS[[p]])][, panel := p]))
for (z in list(f34, ctx)) z[, panel := factor(panel, levels = names(PANELS))]
f34[, family := factor(family, levels = FAM)]
setorder(f34, panel, family, -ell); setorder(ctx, panel, family, -ell)
XMAX <- 0.05
yl <- range(f34[top10_agreement_error <= XMAX, top10_prediction_gain])
p34 <- ggplot(f34, aes(top10_agreement_error, top10_prediction_gain)) +
  geom_path(data = ctx, aes(group = family), colour = "#d9d8d0", linewidth = 0.6) +
  geom_path(data = f34[FAM_LT[as.integer(family)] == "solid"], aes(colour = family, linetype = family),
            linewidth = 0.9, lineend = "round") +
  geom_path(data = f34[FAM_LT[as.integer(family)] != "solid"], aes(colour = family, linetype = family),
            linewidth = 0.9, lineend = "round") +          # dashed on top, so overlaps stay visible
  geom_point(data = f34[ell == 100], aes(colour = family), shape = 21, fill = SURF, stroke = 1, size = 2.4,
             show.legend = FALSE) +
  geom_point(data = f34[is.infinite(ell)], aes(colour = family), size = 1.8, show.legend = FALSE) +
  scale_colour_manual(values = setNames(FAM_COL, FAM)) +
  scale_linetype_manual(values = setNames(FAM_LT, FAM)) +
  facet_wrap(~panel, nrow = 1) +
  coord_cartesian(xlim = c(0.014, XMAX), ylim = yl + c(-0.001, 0.002)) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE), linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(title = "The combined formulas as curves: up and to the left is better",
       subtitle = paste0("Each curve moves from agreement alone (filled dot) toward agreement times co-methylation. ",
                         "Open dot: 100 bp length scale.\nSolid: distance weight. Dashed: one weight for every pair. ",
                         "Light gray: the other formulas. Top 10% of 32,664 RRBS pairs,\nmean over 18 ",
                         "leave-one-patient-out folds, zoomed to the range compared in the tables."),
       x = "Held-out agreement error |x - y| (lower is better)",
       y = "Held-out prediction gain (higher is better)",
       caption = "Source: Results/Task46_Frontier.csv and Task46_Frontier_Matched.csv (Task 46)") +
  theme_minimal(base_size = 10.5) +
  theme(plot.background = element_rect(fill = SURF, colour = NA),
        panel.grid.major = element_line(colour = GRID, linewidth = 0.3), panel.grid.minor = element_blank(),
        plot.title = element_text(colour = INK, face = "bold", size = 13), plot.title.position = "plot",
        plot.subtitle = element_text(colour = INK2, size = 9.5), plot.caption = element_text(colour = MUTED, hjust = 0),
        plot.caption.position = "plot",
        legend.position = "top", legend.justification = "left", legend.title = element_blank(),
        legend.key.width = unit(28, "pt"), panel.spacing.x = unit(16, "pt"),
        strip.text = element_text(colour = INK, face = "bold", hjust = 0),
        axis.text = element_text(colour = INK2), axis.title = element_text(colour = INK2))
ggsave(file.path(res, "Fig34_SimilarityFrontier.png"), p34, width = 13, height = 6, dpi = 200, bg = SURF)
cat("\ndone\n")
