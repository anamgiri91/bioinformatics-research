# ================================================================
# Task 32 - Five ways to make the flag mean something
# Date: Sep 3, 2026
#
# Task 31 closed the door on tuning OutlierMeth: the p argument has
# two usable settings on the 747-sample tcga panel, and even the
# strictest reachable threshold leaves 100% of H-site hyper flags
# below |dbeta| = 0.10. Any fix has to add a term the package does
# not have. This script builds and scores the candidates.
#
# ---------------------------------------------------------------
# PART A - a MEASURED noise floor instead of an assumed one
#
# plan.md step 3b: every stability result so far assumes N(0, 0.01)
# because that is what the literature says about 450k technical
# error, and that assumption is the one thing the reproducibility
# argument rests on. Adjacent CpGs give a way to measure it here.
#
# For two probes i, j a few hundred bp apart, let D_s = beta_i,s -
# beta_j,s across the 53 samples. Whatever systematic methylation
# difference exists between the two positions is constant in s and
# drops out of sd(D). What is left is technical noise at both probes
# plus any genuine sample-to-sample divergence between them, so
#
#     sd(D) >= sqrt(2) * sigma_technical
#     sigma_technical <= sd(D) / sqrt(2)
#
# is an UPPER bound, and the tightest bound comes from the most
# tightly co-methylated pairs - hence the low quantile over pairs
# rather than the mean. Distant pairs on the same chromosome are
# scored too, as a null: if close pairs are not much tighter than
# distant ones, the bound is measuring biology and is worthless.
#
# ---------------------------------------------------------------
# PART B - the methods, all at one flag rate
#
#   ext.floor.med.0.10   ext flag AND |beta - cohort median| >= 0.10.
#                        Task 24's winner; the rate every other method
#                        is tuned to, following Task 26 Part C.
#   ext.floor.med.0.05   the same at the 0.05 the 2026-09-03 email
#                        proposed. Reported at its own rate as well,
#                        since 0.05 vs 0.10 is the actual question.
#   ext.floor.state      ext flag AND |beta - median| >= a floor set
#                        PER STATE from that state's own deviation
#                        distribution. One global parameter, six
#                        thresholds. This is what Task 30 showed a
#                        state-pooled rule is really estimating.
#   ext.floor.noise      ext flag AND |beta - median| >= k * the
#                        measured sigma from Part A, per state.
#                        The same idea with the constant read off the
#                        array instead of off the beta distribution.
#   ext.delt             ext flag AND (beta - threshold) >= d, i.e.
#                        the package's own deltMeth used as a filter.
#                        Tests whether flooring on distance-to-
#                        threshold works as well as distance-to-median.
#   bio.stat.dev         Task 30's state-pooled deviation rule, with
#                        no external reference at all.
#   asin.mad             median +/- k*MAD on the angular transform
#                        phi = 2*asin(sqrt(beta)).
#   mad.beta, mad.mvalue carried over from Task 24 so the three
#                        transforms can be compared at one rate.
#
# The three MAD rows exist to test a specific claim. Task 24 found
# M-values the least stable method tested, which is the opposite of
# the textbook advice, and attributed it to the logit amplifying
# noise at the extremes. The angular transform stabilises binomial
# variance too but is bounded on [0, pi], so it amplifies far less.
# If the explanation is right, stability should fall monotonically
# with how much the transform stretches the extremes - and this
# script measures that stretch directly (amplification.at.flags)
# rather than asserting it.
#
# ---------------------------------------------------------------
# SCORING  (Task 24's, plus one)
#   pct.contrary  share of flags that run against the site's state -
#                 -1 at L/LM or +1 at H/HM. This is the 2026-09-03
#                 email's own diagnostic, so every method is scored
#                 on it directly.
# ================================================================

options(stringsAsFactors = FALSE)
set.seed(20260903)

pick_dir <- function(...) { for (d in c(...)) if (dir.exists(d)) return(d)
  stop("no candidate directory exists") }
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
CONTRARY <- list(L = -1L, LM = -1L, H = 1L, HM = 1L)   # flag AGAINST the state
NOISE_SD <- 0.01
N_REPS   <- 10
EPS      <- 1e-3
PAIR_DISTS <- c(50, 100, 200, 500)
FAR_MIN <- 1e6                     # "distant pair" null

has_dt <- requireNamespace("data.table", quietly = TRUE)
read_any <- function(p, tab = FALSE) {
  if (has_dt) as.data.frame(data.table::fread(p, showProgress = FALSE,
      sep = if (tab) "\t" else ",", header = TRUE, check.names = FALSE))
  else if (tab) read.table(p, header = TRUE, sep = "\t", check.names = FALSE)
  else read.csv(p, check.names = FALSE)
}
write_out <- function(tbl, f) {
  p <- file.path(out_dir, f); write.csv(tbl, p, row.names = FALSE)
  cat("  written: ", p, "  (", nrow(tbl), " rows)\n", sep = "")
}
rowMedian <- function(M) apply(M, 1, median, na.rm = TRUE)
rowMad    <- function(M) apply(M, 1, mad, na.rm = TRUE)
to_M    <- function(B) { b <- pmin(pmax(B, EPS), 1 - EPS); log2(b / (1 - b)) }
to_asin <- function(B) 2 * asin(sqrt(pmin(pmax(B, 0), 1)))
jaccard <- function(a, b) {
  ok <- !is.na(a) & !is.na(b)
  A <- a[ok] != 0; B2 <- b[ok] != 0
  u <- sum(A | B2); if (!u) return(NA_real_)
  sum(A & B2 & a[ok] == b[ok]) / u
}
calibrate <- function(fn, lo, hi, target, decreasing = TRUE, iters = 20) {
  for (i in seq_len(iters)) {
    mid <- (lo + hi) / 2
    r <- mean(fn(mid) != 0, na.rm = TRUE)
    if (is.na(r)) break
    if ((r > target) == decreasing) lo <- mid else hi <- mid
  }
  (lo + hi) / 2
}
f_madz <- function(X, k) {
  m <- rowMedian(X); s <- rowMad(X); s[is.na(s) | s == 0] <- NA
  o <- matrix(0L, nrow(X), ncol(X))
  o[!is.na(s) & X < m - k * s] <- -1L
  o[!is.na(s) & X > m + k * s] <-  1L
  o[is.na(X)] <- NA; o
}
TAIL <- c(L = "upper", LM = "upper", M = "both", HM = "lower", H = "lower", R = "both")
state_pooled_flags <- function(X, st, q) {
  o <- matrix(0L, nrow(X), ncol(X))
  for (s in STATE_ORDER) {
    k <- which(st == s); if (!length(k)) next
    pool <- as.vector(X[k, , drop = FALSE]); pool <- pool[!is.na(pool)]
    if (length(pool) < 20) next
    tl <- TAIL[[s]]; hi <- lo <- NA_real_
    if (tl == "upper") hi <- quantile(pool, q, names = FALSE)
    if (tl == "lower") lo <- quantile(pool, 1 - q, names = FALSE)
    if (tl == "both") { hi <- quantile(pool, 1 - (1 - q)/2, names = FALSE)
                        lo <- quantile(pool, (1 - q)/2, names = FALSE) }
    Xs <- X[k, , drop = FALSE]; os <- matrix(0L, nrow(Xs), ncol(Xs))
    if (is.finite(hi)) os[!is.na(Xs) & Xs > hi] <-  1L
    if (is.finite(lo)) os[!is.na(Xs) & Xs < lo] <- -1L
    o[k, ] <- os
  }
  o[is.na(X)] <- NA; o
}

cat("=========================================================\n")
cat("TASK 32 - IMPROVED FLAGGING RULES\n")
cat("=========================================================\n\n")

noise_rows <- list(); score_rows <- list(); floor_rows <- list()

for (tissue in c("Normal", "Tumor")) {

  cat("#########################################################\n")
  cat("### ", tissue, "\n", sep = "")
  cat("#########################################################\n\n")

  cols <- paste0(substr(tissue, 1, 1), 1:53)
  d <- read_any(inputs[[tissue]], tab = TRUE)
  d <- d[d$Chromosome == "chr22", ]; d <- d[order(d$Start), ]
  ids <- d$Composite.Element.REF; pos <- d$Start
  B <- as.matrix(d[, cols, drop = FALSE]); storage.mode(B) <- "numeric"
  st <- d$methy.state; rm(d); invisible(gc())

  fx <- read_any(flag_ext[[tissue]])
  F <- as.matrix(fx[match(ids, fx$cgID), cols, drop = FALSE])
  storage.mode(F) <- "numeric"; rm(fx); invisible(gc())

  med <- rowMedian(B); dev <- abs(B - med)

  # ================= PART A : measured noise =================
  cat("PART A - technical-noise bound from neighbouring CpG pairs\n")
  cat("  sigma_hat = sd(beta_i - beta_j across samples) / sqrt(2), an UPPER bound.\n")
  cat("  The 10th percentile over pairs is the usable bound; distant pairs are the null.\n\n")

  pair_sigma <- function(i, j) {
    D <- B[i, ] - B[j, ]
    if (sum(!is.na(D)) < 30) return(NA_real_)
    sd(D, na.rm = TRUE) / sqrt(2)
  }

  for (dist in c(PAIR_DISTS, NA)) {
    if (is.na(dist)) {
      # distant-pair null: sample random same-chromosome pairs > 1 Mb apart
      set.seed(20260903)
      n <- length(pos); ii <- sample.int(n, 4000, replace = TRUE)
      jj <- sample.int(n, 4000, replace = TRUE)
      keep <- which(abs(pos[ii] - pos[jj]) > FAR_MIN & st[ii] == st[jj])
      ii <- ii[keep]; jj <- jj[keep]; lab <- "distant (>1 Mb, null)"
    } else {
      ii <- integer(0); jj <- integer(0)
      for (g in seq_len(length(pos) - 1L)) {
        h <- g + 1L
        while (h <= length(pos) && pos[h] - pos[g] <= dist) {
          if (st[h] == st[g]) { ii <- c(ii, g); jj <- c(jj, h) }
          h <- h + 1L
        }
      }
      lab <- paste0("<= ", dist, " bp")
    }
    if (!length(ii)) next
    for (s in c("ALL", STATE_ORDER)) {
      sel <- if (s == "ALL") seq_along(ii) else which(st[ii] == s)
      if (length(sel) < 20) next
      sg <- mapply(pair_sigma, ii[sel], jj[sel])
      sg <- sg[is.finite(sg)]
      if (length(sg) < 20) next
      noise_rows[[length(noise_rows) + 1]] <- data.frame(
        tissue = tissue, pair.set = lab, methy.state = s, n.pairs = length(sg),
        sigma.p10 = round(quantile(sg, .10, names = FALSE), 5),
        sigma.median = round(median(sg), 5),
        sigma.p90 = round(quantile(sg, .90, names = FALSE), 5))
      if (s == "ALL" || (!is.na(dist) && dist == 100))
        cat(sprintf("   %-22s %-4s pairs=%5d   sigma p10=%.4f  median=%.4f  p90=%.4f\n",
            lab, s, length(sg), quantile(sg,.10,names=FALSE), median(sg),
            quantile(sg,.90,names=FALSE)))
    }
  }

  nz <- do.call(rbind, noise_rows)
  nz <- nz[nz$tissue == tissue & nz$pair.set == "<= 100 bp", ]
  sigma_state <- setNames(rep(NA_real_, length(STATE_ORDER)), STATE_ORDER)
  for (s in STATE_ORDER) {
    r <- nz[nz$methy.state == s, ]
    sigma_state[s] <- if (nrow(r)) r$sigma.p10[1] else
      nz$sigma.p10[nz$methy.state == "ALL"][1]
  }
  sigma_state[!is.finite(sigma_state)] <-
    nz$sigma.p10[nz$methy.state == "ALL"][1]
  cat("\n  per-state sigma used for ext.floor.noise: ",
      paste(sprintf("%s=%.4f", STATE_ORDER, sigma_state[STATE_ORDER]),
            collapse = "  "), "\n\n", sep = "")

  # ================= reconstruct ext thresholds ==============
  P <- N <- rep(NA_real_, nrow(B))
  for (i in seq_len(nrow(B))) {
    b <- B[i, ]; fl <- F[i, ]; if (all(is.na(fl))) next
    up <- which(fl == 1); dn <- which(fl == -1); ze <- which(fl == 0)
    loP <- if (length(c(ze, dn))) max(b[c(ze, dn)], na.rm = TRUE) else NA_real_
    hiP <- if (length(up)) min(b[up], na.rm = TRUE) else NA_real_
    if (is.finite(loP) && is.finite(hiP)) P[i] <- (loP + hiP)/2
    loN <- if (length(dn)) max(b[dn], na.rm = TRUE) else NA_real_
    hiN <- if (length(c(ze, up))) min(b[c(ze, up)], na.rm = TRUE) else NA_real_
    if (is.finite(loN) && is.finite(hiN)) N[i] <- (loN + hiN)/2
  }
  # re-applying the recovered thresholds is how every ext-based method
  # is re-flagged under noise: the panel is external and does not move.
  ext_apply <- function(X) {
    o <- matrix(0L, nrow(X), ncol(X))
    o[!is.na(P) & X > P] <-  1L
    o[!is.na(N) & X < N] <- -1L
    o[is.na(X)] <- NA; o
  }

  # ================= method builders =========================
  f_floor_const <- function(Fm, X, off) {
    m <- rowMedian(X); o <- Fm
    o[!is.na(Fm) & Fm != 0 & abs(X - m) < off] <- 0L; o
  }
  f_floor_state <- function(Fm, X, qf) {
    m <- rowMedian(X); dv <- abs(X - m); o <- Fm
    for (s in STATE_ORDER) {
      k <- which(st == s); if (!length(k)) next
      pool <- as.vector(dv[k, , drop = FALSE]); pool <- pool[is.finite(pool)]
      if (!length(pool)) next
      th <- quantile(pool, qf, names = FALSE)
      Os <- o[k, , drop = FALSE]; Dv <- dv[k, , drop = FALSE]
      Os[!is.na(Os) & Os != 0 & Dv < th] <- 0L
      o[k, ] <- Os
    }
    o
  }
  f_floor_noise <- function(Fm, X, kk) {
    m <- rowMedian(X); dv <- abs(X - m); o <- Fm
    for (s in STATE_ORDER) {
      k <- which(st == s); if (!length(k)) next
      th <- kk * sigma_state[[s]]
      Os <- o[k, , drop = FALSE]; Dv <- dv[k, , drop = FALSE]
      Os[!is.na(Os) & Os != 0 & Dv < th] <- 0L
      o[k, ] <- Os
    }
    o
  }
  f_delt <- function(Fm, X, dd) {
    o <- Fm
    hi <- !is.na(Fm) & Fm ==  1 & (is.na(P) | (X - P) < dd)
    lo <- !is.na(Fm) & Fm == -1 & (is.na(N) | (N - X) < dd)
    o[hi] <- 0L; o[lo] <- 0L; o
  }

  base_ext <- ext_apply(B)
  target <- mean(f_floor_const(base_ext, B, 0.10) != 0, na.rm = TRUE)
  cat("PART B - all methods tuned to the ext + 0.10 floor rate: ",
      round(100 * target, 3), "%\n\n", sep = "")

  cat("  calibrating...\n")
  qf_state <- calibrate(function(q) f_floor_state(base_ext, B, q), 0, 0.9999,
                        target, decreasing = TRUE)
  k_noise  <- calibrate(function(k) f_floor_noise(base_ext, B, k), 0, 60, target)
  d_delt   <- calibrate(function(dd) f_delt(base_ext, B, dd), 0, 0.9, target)
  # rate DECREASES as q rises (a higher quantile is a stricter cut), so
  # this bisects with decreasing = TRUE. Getting this backwards pins q at
  # the lower bound and silently reports a ~50% flag rate.
  q_pool   <- calibrate(function(q) state_pooled_flags(B - med, st, q),
                        0.5, 0.999999, target, decreasing = TRUE)
  k_madb   <- calibrate(function(k) f_madz(B, k), 0.1, 40, target)
  k_mada   <- calibrate(function(k) f_madz(to_asin(B), k), 0.1, 40, target)
  k_madm   <- calibrate(function(k) f_madz(to_M(B), k), 0.1, 40, target)
  cat(sprintf("    floor.state q=%.4f  floor.noise k=%.2f  delt d=%.4f  pool q=%.6f\n",
              qf_state, k_noise, d_delt, q_pool))
  cat(sprintf("    mad k: beta=%.2f  asin=%.2f  mvalue=%.2f\n\n",
              k_madb, k_mada, k_madm))

  # record what floor each state actually got, in beta units
  dvB <- abs(B - med)
  for (s in STATE_ORDER) {
    k <- which(st == s); if (!length(k)) next
    pool <- as.vector(dvB[k, , drop = FALSE]); pool <- pool[is.finite(pool)]
    floor_rows[[length(floor_rows) + 1]] <- data.frame(
      tissue = tissue, methy.state = s,
      floor.const = 0.10,
      floor.state = round(quantile(pool, qf_state, names = FALSE), 4),
      sigma.measured = round(sigma_state[[s]], 5),
      floor.noise = round(k_noise * sigma_state[[s]], 4),
      k.noise = round(k_noise, 2), q.state = round(qf_state, 4))
  }

  methods <- list(
    ext.percentile     = base_ext,
    ext.floor.med.0.05 = f_floor_const(base_ext, B, 0.05),
    ext.floor.med.0.10 = f_floor_const(base_ext, B, 0.10),
    ext.floor.state    = f_floor_state(base_ext, B, qf_state),
    ext.floor.noise    = f_floor_noise(base_ext, B, k_noise),
    ext.delt           = f_delt(base_ext, B, d_delt),
    bio.stat.dev       = state_pooled_flags(B - med, st, q_pool),
    mad.beta           = f_madz(B, k_madb),
    asin.mad           = f_madz(to_asin(B), k_mada),
    mad.mvalue         = f_madz(to_M(B), k_madm))

  cat("  stability: ", N_REPS, " replicates at sd = ", NOISE_SD, "\n", sep = "")
  stab <- setNames(numeric(length(methods)), names(methods))
  for (r in seq_len(N_REPS)) {
    Bn <- B + matrix(rnorm(length(B), 0, NOISE_SD), nrow(B))
    Bn[Bn < 0] <- 0; Bn[Bn > 1] <- 1; Bn[is.na(B)] <- NA
    en <- ext_apply(Bn)
    pert <- list(
      ext.percentile     = en,
      ext.floor.med.0.05 = f_floor_const(en, Bn, 0.05),
      ext.floor.med.0.10 = f_floor_const(en, Bn, 0.10),
      ext.floor.state    = f_floor_state(en, Bn, qf_state),
      ext.floor.noise    = f_floor_noise(en, Bn, k_noise),
      ext.delt           = f_delt(en, Bn, d_delt),
      bio.stat.dev       = state_pooled_flags(Bn - rowMedian(Bn), st, q_pool),
      mad.beta           = f_madz(Bn, k_madb),
      asin.mad           = f_madz(to_asin(Bn), k_mada),
      mad.mvalue         = f_madz(to_M(Bn), k_madm))
    for (nm in names(methods)) stab[nm] <- stab[nm] + jaccard(methods[[nm]], pert[[nm]])
  }
  stab <- stab / N_REPS

  # How much each transform T stretches beta where the flags actually are:
  # |dT/dbeta| at the flagged cells, divided by |dT/dbeta| at beta = 0.5, so
  # every transform reads 1 at mid-methylation and the columns are comparable.
  #
  #   beta     dT/db = 1                    -> 1 at b = 0.5
  #   arcsine  dT/db = (b(1-b))^-1/2        -> 2 at b = 0.5
  #   M-value  dT/db = 1/(ln2 * b(1-b))     -> 4/ln2 = 5.7708 at b = 0.5
  #
  # An earlier version divided the M-value derivative by 4 rather than 4/ln2,
  # which inflated its amplification by 1/ln2 = 1.443x. Ordering and the
  # monotone conclusion were unaffected; the number was not.
  AMP_AT_HALF <- c(mad.beta = 1, asin.mad = 2, mad.mvalue = 4 / log(2))
  amp <- function(nm) {
    if (!nm %in% names(AMP_AT_HALF)) return(NA_real_)
    O <- methods[[nm]]; sel <- !is.na(O) & O != 0
    b <- B[sel]; b <- b[is.finite(b)]
    if (!length(b)) return(NA_real_)
    bb <- pmax(b * (1 - b), 1e-9)
    g <- switch(nm,
      mad.beta   = rep(1, length(b)),
      asin.mad   = 1 / sqrt(bb),
      mad.mvalue = 1 / (log(2) * bb))
    round(median(g) / AMP_AT_HALF[[nm]], 2)
  }

  for (nm in names(methods)) {
    O <- methods[[nm]]
    rate_by_state <- sapply(STATE_ORDER, function(s) {
      k <- which(st == s); if (!length(k)) return(NA_real_)
      mean(O[k, ] != 0, na.rm = TRUE) })
    mag_by_state <- sapply(STATE_ORDER, function(s) {
      k <- which(st == s); if (!length(k)) return(NA_real_)
      v <- as.vector(O[k, ]); w <- as.vector(dev[k, ])
      ok <- !is.na(v) & v != 0 & !is.na(w)
      if (!any(ok)) return(NA_real_); median(w[ok]) })
    sel <- !is.na(O) & O != 0
    mg <- dev[sel]; mg <- mg[is.finite(mg)]
    ncontrary <- 0L
    for (s in names(CONTRARY)) {
      k <- which(st == s); if (!length(k)) next
      ncontrary <- ncontrary + sum(O[k, , drop = FALSE] == CONTRARY[[s]], na.rm = TRUE)
    }
    score_rows[[length(score_rows) + 1]] <- data.frame(
      tissue = tissue, method = nm,
      flag.rate.pct = round(100 * mean(O != 0, na.rm = TRUE), 3),
      flags = sum(sel),
      stability.jaccard = round(stab[[nm]], 4),
      state.rate.ratio = round(max(rate_by_state, na.rm = TRUE) /
                               max(1e-9, min(rate_by_state, na.rm = TRUE)), 2),
      magnitude.H.vs.R = round(mag_by_state[["H"]] / mag_by_state[["R"]], 3),
      median.abs.dbeta = if (length(mg)) round(median(mg), 4) else NA_real_,
      pct.trivial = if (length(mg)) round(100 * mean(mg < 0.10), 2) else NA_real_,
      pct.contrary = round(100 * ncontrary / max(1, sum(sel)), 2),
      amplification.at.flags = amp(nm))
  }

  sc <- do.call(rbind, score_rows); sc <- sc[sc$tissue == tissue, ]
  sc <- sc[order(-sc$stability.jaccard), ]
  cat("\n  ---- scores, best stability first ----\n")
  print(sc, row.names = FALSE)
  cat("\n")

  rm(B, F, dev); invisible(gc())
}

write_out(do.call(rbind, noise_rows), "Task32_MeasuredNoise.csv")
write_out(do.call(rbind, floor_rows), "Task32_FloorsApplied.csv")
write_out(do.call(rbind, score_rows), "Task32_ImprovedMethodScores.csv")

cat("\n---- transform comparison (all at the same flag rate) ----\n")
sc <- do.call(rbind, score_rows)
print(sc[sc$method %in% c("mad.beta", "asin.mad", "mad.mvalue"),
         c("tissue","method","amplification.at.flags","stability.jaccard",
           "median.abs.dbeta","pct.trivial","pct.contrary")], row.names = FALSE)

cat("\n=========================================================\n")
cat("TASK 32 COMPLETE\n")
cat("=========================================================\n")
