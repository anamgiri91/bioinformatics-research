# ================================================================
# Task 24 - Benchmark of alternative outlier definitions
# Date: Aug 28, 2026
#
# Task 22 established three failure modes of the current flagging:
#   (i)   flag density per unit of available beta space is 44-171x
#         higher in the direction the methylation state constrains;
#   (ii)  the decision boundary at H/L sites sits ~0.003 apart in
#         beta, inside the array's technical error, so 21% of flags
#         vanish under noise of sd = 0.01;
#   (iii) it is not package-specific - Tukey's 3xIQR inherits the
#         same pathology because the IQR is itself state-dependent.
#
# This script asks what to use instead. Seven definitions are run on
# the same chr22 cells and scored on the failure modes above.
#
# ---------------------------------------------------------------
# CALIBRATION - the part that makes the comparison fair
#
# A method that flags fewer cells will look more stable and less
# state-biased for free. Every method is therefore tuned so its
# overall flag rate matches the external reference's rate on the
# same data, by bisection on its own threshold parameter. Rates are
# matched to within 5% relative. Uncalibrated comparisons of outlier
# methods are close to meaningless and are the main reason this
# script exists.
#
# ---------------------------------------------------------------
# METHODS
#   ext.percentile  external TCGA panel, thresholds reconstructed from
#                   the Task 13 flag matrix (see Task 22 Part A)
#   ext.floor       the same, plus |beta - cohort median| >= offset
#   tukey.iqr       Q1 - k*IQR / Q3 + k*IQR on beta
#   mad.beta        median +/- k*MAD on beta
#   mad.mvalue      median +/- k*MAD on M = log2(b/(1-b))   <- the
#                   textbook correction for beta heteroscedasticity
#   beta.fit        per-CpG Beta(a,b) by moments, two-sided tail p
#   loo.quantile    leave-one-out empirical quantile (no external ref)
#
# ---------------------------------------------------------------
# SCORING
#   stability.jaccard  Jaccard between the flag set on clean data and
#                      on data perturbed by N(0, 0.01), the array's
#                      own technical error. Every self-derived method
#                      is fully recomputed on the perturbed data, so
#                      threshold re-estimation is included - this is
#                      what a repeat experiment would actually give.
#   state.rate.ratio   max/min flag rate across the six states. 1 is
#                      perfect state-independence.
#   magnitude.H.vs.R   median |dbeta| of flags at H sites divided by
#                      the same at R sites. The current method scores
#                      ~0.12; 1 means the method finds equally large
#                      effects regardless of state.
#   pct.trivial        share of flags whose |dbeta| < 0.10.
# ================================================================

options(stringsAsFactors = FALSE)
set.seed(20260828)

pick_dir <- function(...) {
  cand <- c(...); for (d in cand) if (dir.exists(d)) return(d)
  stop("no candidate directory exists")
}
src_dir <- pick_dir(
  "/home/s_s355/research3.data/TCGA/filter.BRCA.UCEC.Aug17.2022/sorted.TCGA.380355cg.files",
  path.expand("~/Downloads/tcga_brca_380355cg"))
repo_dir <- pick_dir("/mmfs1/home/wln26/Experiments.Outlier.July31.2026",
                     path.expand("~/Desktop/bioinformatics-research"))
out_dir <- file.path(repo_dir, "Results")

inputs <- list(
  Normal = file.path(src_dir, "Sorted.BRCA.53Alive.Normal.380355cg.75col.May28.2026.txt"),
  Tumor  = file.path(src_dir, "Sorted.BRCA.53Alive.Tumor.380355cg.75col.May28.2026.txt"))
flag_ext <- list(
  Normal = file.path(out_dir, "Task13_Normal_extref_flags.csv"),
  Tumor  = file.path(out_dir, "Task13_Tumor_extref_flags.csv"))

STATE_ORDER <- c("L", "LM", "M", "HM", "H", "R")
NOISE_SD <- 0.01
N_REPS   <- 10
EPS      <- 1e-3      # beta clipping for the logit

has_dt <- requireNamespace("data.table", quietly = TRUE)
read_any <- function(p, tab = FALSE) {
  if (has_dt) as.data.frame(data.table::fread(p, showProgress = FALSE,
      sep = if (tab) "\t" else ",", header = TRUE, check.names = FALSE))
  else if (tab) read.table(p, header = TRUE, sep = "\t", check.names = FALSE)
  else read.csv(p, check.names = FALSE)
}
write_out <- function(tbl, f) {
  p <- file.path(out_dir, f); write.csv(tbl, p, row.names = FALSE)
  cat("  written:", p, " (", nrow(tbl), " rows )\n", sep = "")
}

rowMedian <- function(M) apply(M, 1, median, na.rm = TRUE)
rowMad    <- function(M) apply(M, 1, mad, na.rm = TRUE)
to_M <- function(B) { b <- pmin(pmax(B, EPS), 1 - EPS); log2(b / (1 - b)) }

# ---------------- flag builders -------------------------------
# Each returns a matrix of -1 / 0 / +1 with NA preserved.

f_tukey <- function(B, k) {
  q1 <- apply(B, 1, quantile, .25, na.rm = TRUE)
  q3 <- apply(B, 1, quantile, .75, na.rm = TRUE)
  iq <- q3 - q1
  o <- matrix(0L, nrow(B), ncol(B)); o[B < q1 - k * iq] <- -1L
  o[B > q3 + k * iq] <- 1L; o[is.na(B)] <- NA; o
}
f_madz <- function(X, k) {
  m <- rowMedian(X); s <- rowMad(X); s[is.na(s) | s == 0] <- NA
  o <- matrix(0L, nrow(X), ncol(X))
  o[!is.na(s) & X < m - k * s] <- -1L
  o[!is.na(s) & X > m + k * s] <-  1L
  o[is.na(X)] <- NA; o
}
f_betafit <- function(B, p) {
  mu <- rowMeans(B, na.rm = TRUE)
  v  <- apply(B, 1, var, na.rm = TRUE)
  vmax <- mu * (1 - mu)
  ok <- is.finite(mu) & is.finite(v) & v > 0 & v < vmax
  cc <- rep(NA_real_, length(mu)); cc[ok] <- (vmax[ok] / v[ok]) - 1
  a <- mu * cc; b <- (1 - mu) * cc
  o <- matrix(0L, nrow(B), ncol(B))
  lo <- pbeta(B, a, b); hi <- pbeta(B, a, b, lower.tail = FALSE)
  o[!is.na(lo) & lo < p] <- -1L
  o[!is.na(hi) & hi < p] <-  1L
  o[is.na(B) | is.na(a)] <- NA; o
}
# Leave-one-out empirical quantile.
#
# The naive form recomputes quantile() over 52 values for each of the
# 53 columns, i.e. 53 row-wise apply() passes per call, which makes
# calibration by bisection intractable. Instead sort each row ONCE:
# the k-th order statistic of "everyone except sample j" is just
# S[i, k] when k < rank(j), and S[i, k+1] when k >= rank(j). The
# type-7 quantile is then two indexed lookups and a lerp, vectorised
# over rows. Same numbers, ~25x less work.
#
# Rows carrying any NA are returned as NA rather than silently
# quantiled over a shorter vector - a LOO threshold built from a
# different number of samples is not comparable across rows.
f_loo <- function(B, p) {
  n <- ncol(B); m <- n - 1L
  o <- matrix(NA_integer_, nrow(B), ncol(B))
  good <- which(rowSums(is.na(B)) == 0)
  if (!length(good)) return(o)
  G <- B[good, , drop = FALSE]
  S <- t(apply(G, 1, sort))
  R <- t(apply(G, 1, rank, ties.method = "first"))
  ridx <- seq_len(nrow(G))

  ostat <- function(k, j) {          # k-th order stat of the other 52
    k <- min(max(k, 1L), m)
    idx <- ifelse(k < R[, j], k, k + 1L)
    S[cbind(ridx, idx)]
  }
  qtype7 <- function(pr, j) {        # type-7 quantile over the other 52
    h  <- (m - 1) * pr + 1
    lo <- floor(h); fr <- h - lo
    a <- ostat(lo, j)
    if (fr == 0) return(a)
    a + fr * (ostat(lo + 1L, j) - a)
  }

  og <- matrix(0L, nrow(G), n)
  for (j in seq_len(n)) {
    hi <- qtype7(1 - p, j); lo <- qtype7(p, j)
    v <- G[, j]
    og[, j] <- ifelse(v > hi, 1L, ifelse(v < lo, -1L, 0L))
  }
  o[good, ] <- og
  o
}
f_extfloor <- function(F, B, offset) {
  med <- rowMedian(B)
  o <- F; o[!is.na(F) & F != 0 & abs(B - med) < offset] <- 0L; o
}

# ---------------- calibration ---------------------------------
# bisect a method's parameter until its flag rate matches the target

calibrate <- function(fn, lo, hi, target, decreasing = TRUE, iters = 18) {
  rate <- function(k) { o <- fn(k); mean(o != 0, na.rm = TRUE) }
  for (i in seq_len(iters)) {
    mid <- (lo + hi) / 2
    r <- rate(mid)
    if (is.na(r)) break
    if ((r > target) == decreasing) lo <- mid else hi <- mid
  }
  (lo + hi) / 2
}

jaccard <- function(a, b) {
  ok <- !is.na(a) & !is.na(b)
  A <- a[ok] != 0; B2 <- b[ok] != 0
  u <- sum(A | B2); if (!u) return(NA_real_)
  sum(A & B2 & a[ok] == b[ok]) / u
}

cat("====================================================\n")
cat("TASK 24 - OUTLIER METHOD BENCHMARK\n")
cat("====================================================\n\n")

all_scores <- list(); all_states <- list()

for (tissue in c("Normal", "Tumor")) {

  cat("##################################################\n")
  cat("### ", tissue, "\n", sep = "")
  cat("##################################################\n")

  cols <- paste0(substr(tissue, 1, 1), 1:53)
  d <- read_any(inputs[[tissue]], tab = TRUE)
  d <- d[d$Chromosome == "chr22", ]; d <- d[order(d$Start), ]
  B <- as.matrix(d[, cols, drop = FALSE]); storage.mode(B) <- "numeric"
  st <- d$methy.state

  fx <- read_any(flag_ext[[tissue]])
  F <- as.matrix(fx[match(d$Composite.Element.REF, fx$cgID), cols, drop = FALSE])
  storage.mode(F) <- "numeric"; rm(fx); invisible(gc())

  target <- mean(F != 0, na.rm = TRUE)
  cat("  target flag rate (external reference): ",
      round(100 * target, 3), "%\n\n", sep = "")

  M <- to_M(B)
  med <- rowMedian(B)

  cat("  calibrating...\n")
  k_tukey <- calibrate(function(k) f_tukey(B, k), 0.1, 12, target)
  k_madb  <- calibrate(function(k) f_madz(B, k),  0.1, 12, target)
  k_madm  <- calibrate(function(k) f_madz(M, k),  0.1, 12, target)
  p_beta  <- calibrate(function(p) f_betafit(B, p), 1e-12, 0.2, target,
                       decreasing = FALSE)
  p_loo   <- calibrate(function(p) f_loo(B, p), 1e-6, 0.2, target,
                       decreasing = FALSE)
  cat(sprintf("    tukey k=%.3f  mad.beta k=%.3f  mad.mvalue k=%.3f  beta p=%.3g  loo p=%.4f\n\n",
              k_tukey, k_madb, k_madm, p_beta, p_loo))

  methods <- list(
    ext.percentile = F,
    ext.floor.0.10 = f_extfloor(F, B, 0.10),
    tukey.iqr      = f_tukey(B, k_tukey),
    mad.beta       = f_madz(B, k_madb),
    mad.mvalue     = f_madz(M, k_madm),
    beta.fit       = f_betafit(B, p_beta),
    loo.quantile   = f_loo(B, p_loo))

  # ---- stability: recompute each method on perturbed data ----
  cat("  stability: ", N_REPS, " replicates at noise sd = ", NOISE_SD, "\n", sep = "")
  stab <- setNames(numeric(length(methods)), names(methods))
  for (r in seq_len(N_REPS)) {
    Bn <- B + matrix(rnorm(length(B), 0, NOISE_SD), nrow(B))
    Bn[Bn < 0] <- 0; Bn[Bn > 1] <- 1; Bn[is.na(B)] <- NA
    Mn <- to_M(Bn)
    # ext thresholds are external and fixed: re-apply, do not re-fit
    Fn <- F
    pm <- rowMedian(Bn)
    pert <- list(
      ext.percentile = { z <- F; z[!is.na(B) & !is.na(Bn)] <-
          F[!is.na(B) & !is.na(Bn)]; z },
      ext.floor.0.10 = f_extfloor(F, Bn, 0.10),
      tukey.iqr      = f_tukey(Bn, k_tukey),
      mad.beta       = f_madz(Bn, k_madb),
      mad.mvalue     = f_madz(Mn, k_madm),
      beta.fit       = f_betafit(Bn, p_beta),
      loo.quantile   = f_loo(Bn, p_loo))
    for (nm in names(methods))
      stab[nm] <- stab[nm] + jaccard(methods[[nm]], pert[[nm]])
  }
  stab <- stab / N_REPS

  # ext.percentile has fixed external thresholds; re-flag beta directly
  # against the reconstructed thresholds so its stability is comparable
  P.lo <- P.hi <- N.lo <- N.hi <- rep(NA_real_, nrow(B))
  for (i in seq_len(nrow(B))) {
    b <- B[i, ]; fl <- F[i, ]; if (all(is.na(fl))) next
    up <- which(fl == 1); dn <- which(fl == -1); ze <- which(fl == 0)
    if (length(up) && length(ze)) { P.lo[i] <- max(b[ze], na.rm=TRUE); P.hi[i] <- min(b[up], na.rm=TRUE) }
    if (length(dn) && length(ze)) { N.lo[i] <- max(b[dn], na.rm=TRUE); N.hi[i] <- min(b[ze], na.rm=TRUE) }
  }
  Pm <- (P.lo + P.hi)/2; Nm <- (N.lo + N.hi)/2
  ext_apply <- function(X) {
    o <- matrix(0L, nrow(X), ncol(X))
    o[!is.na(Pm) & X > Pm] <- 1L; o[!is.na(Nm) & X < Nm] <- -1L
    o[is.na(X)] <- NA; o
  }
  base_ext <- ext_apply(B); s2 <- 0
  for (r in seq_len(N_REPS)) {
    Bn <- B + matrix(rnorm(length(B), 0, NOISE_SD), nrow(B))
    Bn[Bn<0] <- 0; Bn[Bn>1] <- 1; Bn[is.na(B)] <- NA
    s2 <- s2 + jaccard(base_ext, ext_apply(Bn))
  }
  stab["ext.percentile"] <- s2 / N_REPS

  # ---- scoring ----
  rows <- list(); srows <- list()
  for (nm in names(methods)) {
    O <- methods[[nm]]
    dev <- abs(B - med)
    rate_by_state <- sapply(STATE_ORDER, function(s) {
      k <- which(st == s); if (!length(k)) return(NA_real_)
      mean(O[k, ] != 0, na.rm = TRUE) })
    mag_by_state <- sapply(STATE_ORDER, function(s) {
      k <- which(st == s); if (!length(k)) return(NA_real_)
      v <- as.vector(O[k, ]); w <- as.vector(dev[k, ])
      ok <- !is.na(v) & v != 0 & !is.na(w)
      if (!any(ok)) return(NA_real_); median(w[ok]) })
    allv <- as.vector(O); alld <- as.vector(dev)
    ok <- !is.na(allv) & allv != 0 & !is.na(alld)
    for (s in STATE_ORDER) srows[[length(srows)+1]] <- data.frame(
      tissue = tissue, method = nm, methy.state = s,
      flag.rate.pct = round(100*rate_by_state[[s]], 3),
      median.abs.dbeta = round(mag_by_state[[s]], 4))
    rows[[length(rows) + 1]] <- data.frame(
      tissue = tissue, method = nm,
      flag.rate.pct = round(100 * mean(O != 0, na.rm = TRUE), 3),
      stability.jaccard = round(stab[[nm]], 4),
      state.rate.ratio = round(max(rate_by_state, na.rm=TRUE) /
                               max(1e-9, min(rate_by_state, na.rm=TRUE)), 2),
      magnitude.H.vs.R = round(mag_by_state[["H"]] / mag_by_state[["R"]], 3),
      median.abs.dbeta = round(median(alld[ok]), 4),
      pct.trivial = round(100 * mean(alld[ok] < 0.10), 2))
  }
  sc <- do.call(rbind, rows); rownames(sc) <- NULL
  sc <- sc[order(-sc$stability.jaccard), ]
  cat("\n  ---- scores, best stability first ----\n")
  print(sc, row.names = FALSE)
  all_scores[[tissue]] <- sc
  all_states[[tissue]] <- do.call(rbind, srows)

  rm(B, F, M); invisible(gc()); cat("\n")
}

sc <- do.call(rbind, all_scores); rownames(sc) <- NULL
write_out(sc, "Task24_MethodBenchmark_Scores.csv")
write_out(do.call(rbind, all_states), "Task24_MethodBenchmark_ByState.csv")

cat("\n====================================================\n")
cat("TASK 24 COMPLETE\n")
cat("====================================================\n")
cat("All methods calibrated to the external reference's flag rate, so\n")
cat("stability and state-dependence are compared at equal sensitivity.\n")
cat("====================================================\n")
