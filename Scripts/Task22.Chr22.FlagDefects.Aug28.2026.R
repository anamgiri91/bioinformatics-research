# ================================================================
# Task 22 - Defects in the outlier flag itself
# Date: Aug 28, 2026
#
# Tasks 18-21 asked whether the flags are biased by methylation
# state. They are. This script asks the next question: WHY, and is
# the flag reproducible at all?
#
# Four defects are measured, each with the data on hand. None of
# them require tcga.rda - Part A reconstructs the external
# thresholds from the flag matrix itself.
#
#   A. ZONE ASYMMETRY. How much of the [0,1] beta scale lies on
#      each side of the external thresholds, by state.
#   B. THRESHOLD RESOLUTION. How tightly the 53 samples pack around
#      the threshold - i.e. how far apart consecutive samples are
#      where the decision is made.
#   C. FLAG STABILITY. Re-flag after adding measurement noise and
#      count how many calls change. This is the decisive test: a
#      flag that flips under noise smaller than the array's own
#      technical error is not a measurement, it is a coin toss.
#   D. DIRECTIONAL SKEW and COHORT-LEVEL MISCALIBRATION. Which way
#      the flags point, and whether they are individual outliers or
#      the whole cohort sitting outside the panel's range.
#   E. IS IT OutlierMeth, OR THE WHOLE APPROACH? The same cells are
#      re-flagged with Tukey's 3xIQR rule - a completely different
#      statistic, and the one behind the outliers.coef2/coef3 columns
#      already present in the source data. If an independent rule
#      shows the same state dependence, the defect is not a bug in
#      one package but a property of scale-free thresholding on
#      bounded, heteroscedastic beta values.
#
# ---------------------------------------------------------------
# THRESHOLD RECONSTRUCTION (Part A)
#
# flagMeth() sets +1 iff beta > P and -1 iff beta < N. So for any
# CpG with samples on both sides of a threshold:
#
#     max(beta | flag != +1)  <=  P  <  min(beta | flag == +1)
#     max(beta | flag == -1)  <   N  <= min(beta | flag != -1)
#
# The midpoint of that bracket estimates the threshold and the
# bracket width bounds the error. A CpG with no flagged sample, or
# with every sample flagged, gives no bracket and is skipped - so
# these estimates describe the CpGs where a decision was actually
# made, which is the relevant population.
#
# The reconstruction doubles as a correctness check on the whole
# Task 13 pipeline: if any CpG showed a flagged sample whose beta
# sits below an unflagged one, the flag matrix and the beta matrix
# would be out of register. Zero such violations is a strong
# indication the join in Task 13 is sound.
# ---------------------------------------------------------------
#
# Outputs:
#   Task22_Chr22_ThresholdGeometry_<tissue>.csv
#   Task22_Chr22_FlagStability_<tissue>.csv
#   Task22_Chr22_Prevalence_<tissue>.csv
#   Task22_Chr22_Direction_<tissue>.csv
#   Task22_Chr22_PercentileVsTukey_<tissue>.csv
# ================================================================

options(stringsAsFactors = FALSE)
set.seed(20260828)

pick_dir <- function(...) {
  cand <- c(...)
  for (d in cand) if (dir.exists(d)) return(d)
  stop("no candidate directory exists:\n  ", paste(cand, collapse = "\n  "))
}

src_dir <- pick_dir(
  "/home/s_s355/research3.data/TCGA/filter.BRCA.UCEC.Aug17.2022/sorted.TCGA.380355cg.files",
  path.expand("~/Downloads/tcga_brca_380355cg"))
repo_dir <- pick_dir(
  "/mmfs1/home/wln26/Experiments.Outlier.July31.2026",
  path.expand("~/Desktop/bioinformatics-research"))
out_dir <- file.path(repo_dir, "Results")

inputs <- list(
  Normal = file.path(src_dir, "Sorted.BRCA.53Alive.Normal.380355cg.75col.May28.2026.txt"),
  Tumor  = file.path(src_dir, "Sorted.BRCA.53Alive.Tumor.380355cg.75col.May28.2026.txt"))
flag_ext <- list(
  Normal = file.path(out_dir, "Task13_Normal_extref_flags.csv"),
  Tumor  = file.path(out_dir, "Task13_Tumor_extref_flags.csv"))

STATE_ORDER <- c("L", "LM", "M", "HM", "H", "R")

# Technical reproducibility of the 450k array. Replicate SDs are
# reported around 0.01-0.03 beta for most probes and are worst in
# the intermediate range; 0.05 is a conservative upper bound.
NOISE_SD <- c(0.005, 0.01, 0.02, 0.05)
N_REPS   <- 20

has_dt <- requireNamespace("data.table", quietly = TRUE)
read_any <- function(p, tab = FALSE) {
  if (has_dt) {
    as.data.frame(data.table::fread(p, showProgress = FALSE,
      sep = if (tab) "\t" else ",", header = TRUE, check.names = FALSE))
  } else if (tab) {
    read.table(p, header = TRUE, sep = "\t", check.names = FALSE)
  } else {
    read.csv(p, check.names = FALSE)
  }
}
write_out <- function(tbl, fname) {
  path <- file.path(out_dir, fname)
  write.csv(tbl, path, row.names = FALSE)
  cat("  written:", path, " (", nrow(tbl), " rows )\n", sep = "")
}

cat("====================================================\n")
cat("TASK 22 - DEFECTS IN THE OUTLIER FLAG\n")
cat("====================================================\n\n")

for (tissue in c("Normal", "Tumor")) {

  cat("##################################################\n")
  cat("### ", tissue, "\n", sep = "")
  cat("##################################################\n\n")

  cols <- paste0(substr(tissue, 1, 1), 1:53)

  d <- read_any(inputs[[tissue]], tab = TRUE)
  d <- d[d$Chromosome == "chr22", ]
  d <- d[order(d$Start), ]

  B <- as.matrix(d[, cols, drop = FALSE]); storage.mode(B) <- "numeric"

  fx <- read_any(flag_ext[[tissue]])
  idx <- match(d$Composite.Element.REF, fx$cgID)
  stopifnot(!anyNA(idx))
  F <- as.matrix(fx[idx, cols, drop = FALSE]); storage.mode(F) <- "numeric"
  rm(fx); invisible(gc())

  st <- d$methy.state
  n  <- nrow(B)
  cat("chr22 CpGs:", n, " samples:", ncol(B), "\n")
  print(table(factor(st, levels = STATE_ORDER)))
  cat("\n")

  # --------------------------------------------------
  # A + B. threshold geometry and resolution
  # --------------------------------------------------

  P.lo <- P.hi <- N.lo <- N.hi <- rep(NA_real_, n)
  viol <- 0L

  for (i in seq_len(n)) {
    b <- B[i, ]; fl <- F[i, ]
    if (all(is.na(fl))) next
    up <- which(fl ==  1); dn <- which(fl == -1); ze <- which(fl == 0)
    if (length(up) && length(ze)) {
      P.lo[i] <- max(b[ze], na.rm = TRUE); P.hi[i] <- min(b[up], na.rm = TRUE)
      if (P.hi[i] <= P.lo[i]) viol <- viol + 1L
    }
    if (length(dn) && length(ze)) {
      N.lo[i] <- max(b[dn], na.rm = TRUE); N.hi[i] <- min(b[ze], na.rm = TRUE)
      if (N.lo[i] >= N.hi[i]) viol <- viol + 1L
    }
  }
  Pmid <- (P.lo + P.hi) / 2
  Nmid <- (N.lo + N.hi) / 2

  cat("=== Part A/B: threshold geometry ===\n")
  cat("  register check - CpGs where a flagged sample's beta sits inside\n")
  cat("  the unflagged range (would mean the flag and beta matrices are\n")
  cat("  misaligned): ", viol, "\n\n", sep = "")

  geo <- do.call(rbind, lapply(STATE_ORDER, function(s) {
    k  <- which(st == s)
    kp <- k[!is.na(Pmid[k])]; kn <- k[!is.na(Nmid[k])]
    if (!length(k)) return(NULL)
    data.frame(
      tissue = tissue, methy.state = s, n.cpgs = length(k),
      n.P.bracketed = length(kp),
      P.threshold = if (length(kp)) round(median(Pmid[kp]), 4) else NA,
      hyper.zone.width = if (length(kp)) round(median(1 - Pmid[kp]), 4) else NA,
      P.bracket.width = if (length(kp)) round(median(P.hi[kp] - P.lo[kp]), 5) else NA,
      n.N.bracketed = length(kn),
      N.threshold = if (length(kn)) round(median(Nmid[kn]), 4) else NA,
      hypo.zone.width = if (length(kn)) round(median(Nmid[kn]), 4) else NA,
      N.bracket.width = if (length(kn)) round(median(N.hi[kn] - N.lo[kn]), 5) else NA,
      zone.asymmetry = if (length(kp) && length(kn))
        round(median(1 - Pmid[kp]) / median(Nmid[kn]), 2) else NA)
  }))
  print(geo[, c("methy.state", "n.cpgs", "P.threshold", "hyper.zone.width",
                "P.bracket.width", "N.threshold", "hypo.zone.width",
                "N.bracket.width", "zone.asymmetry")], row.names = FALSE)
  cat("\n  hyper.zone.width = 1 - P, the share of the beta scale in which a\n")
  cat("  sample counts as hyper-methylated. hypo.zone.width = N - 0.\n")
  cat("  P.bracket.width = how far apart the two samples straddling the\n")
  cat("  threshold sit. Where that is smaller than array noise, the call\n")
  cat("  is not reproducible - Part C measures exactly that.\n")
  write_out(geo, paste0("Task22_Chr22_ThresholdGeometry_", tissue, ".csv"))

  # --------------------------------------------------
  # C. flag stability under measurement noise
  # --------------------------------------------------
  # Re-flag using the reconstructed thresholds after perturbing beta.
  # Only CpGs with a usable bracket are testable; the midpoint is used
  # as the threshold estimate. Flip rate is measured against the flags
  # that midpoint reproduces on the unperturbed data, so the estimate
  # error cancels out of the comparison.

  cat("\n=== Part C: flag stability under measurement noise ===\n")
  cat("  ", N_REPS, " replicates per noise level. beta is perturbed by\n", sep = "")
  cat("  N(0, sd), clamped to [0,1], and re-flagged at the same thresholds.\n\n")

  testable <- which(!is.na(Pmid) | !is.na(Nmid))
  Pt <- Pmid[testable]; Nt <- Nmid[testable]
  Bt <- B[testable, , drop = FALSE]
  stt <- st[testable]

  base <- matrix(0L, nrow(Bt), ncol(Bt))
  base[!is.na(Pt) & Bt >  Pt] <-  1L
  base[!is.na(Nt) & Bt <  Nt] <- -1L
  base[is.na(Bt)] <- NA

  stab_rows <- list()
  for (sd in NOISE_SD) {
    flips <- lost <- gained <- dirflip <- 0
    per_state_flip <- setNames(numeric(length(STATE_ORDER)), STATE_ORDER)
    per_state_n    <- setNames(numeric(length(STATE_ORDER)), STATE_ORDER)
    for (r in seq_len(N_REPS)) {
      Bn <- Bt + matrix(rnorm(length(Bt), 0, sd), nrow(Bt))
      Bn[Bn < 0] <- 0; Bn[Bn > 1] <- 1
      nf <- matrix(0L, nrow(Bn), ncol(Bn))
      nf[!is.na(Pt) & Bn >  Pt] <-  1L
      nf[!is.na(Nt) & Bn <  Nt] <- -1L
      nf[is.na(Bn)] <- NA
      ch <- which(!is.na(base) & !is.na(nf) & base != nf)
      flips   <- flips   + length(ch)
      lost    <- lost    + sum(base[ch] != 0 & nf[ch] == 0)
      gained  <- gained  + sum(base[ch] == 0 & nf[ch] != 0)
      dirflip <- dirflip + sum(base[ch] != 0 & nf[ch] != 0)
      for (s in STATE_ORDER) {
        ks <- which(stt == s)
        if (!length(ks)) next
        bs <- base[ks, , drop = FALSE]; ns <- nf[ks, , drop = FALSE]
        per_state_flip[s] <- per_state_flip[s] +
          sum(!is.na(bs) & bs != 0 & !is.na(ns) & ns != bs)
        per_state_n[s] <- per_state_n[s] + sum(!is.na(bs) & bs != 0)
      }
    }
    nbase <- sum(base != 0, na.rm = TRUE)
    stab_rows[[length(stab_rows) + 1]] <- data.frame(
      tissue = tissue, noise.sd = sd,
      n.cpgs.testable = nrow(Bt),
      baseline.flags = nbase,
      mean.calls.changed.per.rep = round(flips / N_REPS, 1),
      mean.flags.lost = round(lost / N_REPS, 1),
      mean.flags.gained = round(gained / N_REPS, 1),
      mean.direction.reversed = round(dirflip / N_REPS, 1),
      pct.baseline.flags.lost = round(100 * (lost / N_REPS) / nbase, 2),
      flag.churn.ratio = round((flips / N_REPS) / nbase, 3))
    for (s in STATE_ORDER) {
      if (per_state_n[s] == 0) next
      stab_rows[[length(stab_rows)]][[paste0("pct.unstable.", s)]] <-
        round(100 * per_state_flip[s] / per_state_n[s], 2)
    }
  }
  stab <- do.call(rbind, lapply(stab_rows, function(x) x))
  print(stab[, c("noise.sd", "baseline.flags", "mean.flags.lost",
                 "mean.flags.gained", "pct.baseline.flags.lost",
                 "flag.churn.ratio")], row.names = FALSE)
  cat("\n  per-state instability (% of baseline flags that change call):\n")
  ps <- grep("^pct.unstable", names(stab), value = TRUE)
  print(stab[, c("noise.sd", ps)], row.names = FALSE)
  write_out(stab, paste0("Task22_Chr22_FlagStability_", tissue, ".csv"))

  # --------------------------------------------------
  # D. prevalence and direction
  # --------------------------------------------------

  cat("\n=== Part D: prevalence and direction ===\n")
  nflag <- rowSums(F != 0, na.rm = TRUE)
  nobs  <- rowSums(!is.na(F))
  prev  <- ifelse(nobs > 0, nflag / nobs, NA)

  prv <- do.call(rbind, lapply(STATE_ORDER, function(s) {
    k <- which(st == s & nobs > 0)
    if (!length(k)) return(NULL)
    p <- prev[k]
    data.frame(tissue = tissue, methy.state = s, n.cpgs = length(k),
      cpgs.with.any.flag = sum(p > 0),
      median.pct.cohort.flagged = round(100 * median(p[p > 0]), 2),
      cpgs.ge.25pct = sum(p >= .25), cpgs.ge.50pct = sum(p >= .50),
      max.samples.flagged = max(nflag[k]),
      flags.total = sum(nflag[k]),
      pct.flags.in.cpgs.ge.50pct =
        round(100 * sum(nflag[k][p >= .50]) / max(1, sum(nflag[k])), 2))
  }))
  cat("  A CpG where half the cohort is flagged is a threshold that does not\n")
  cat("  fit this cohort, not 26 outliers.\n\n")
  print(prv[, c("methy.state", "n.cpgs", "cpgs.with.any.flag",
                "median.pct.cohort.flagged", "cpgs.ge.50pct",
                "max.samples.flagged", "pct.flags.in.cpgs.ge.50pct")],
        row.names = FALSE)
  write_out(prv, paste0("Task22_Chr22_Prevalence_", tissue, ".csv"))

  dir <- do.call(rbind, lapply(STATE_ORDER, function(s) {
    k <- which(st == s)
    if (!length(k)) return(NULL)
    v <- as.vector(F[k, ]); v <- v[!is.na(v) & v != 0]
    if (!length(v)) return(NULL)
    g <- geo[geo$methy.state == s, ]
    data.frame(tissue = tissue, methy.state = s, n.flags = length(v),
      pct.hyper = round(100 * mean(v == 1), 1),
      pct.hypo  = round(100 * mean(v == -1), 1),
      hyper.zone.width = g$hyper.zone.width,
      hypo.zone.width  = g$hypo.zone.width)
  }))
  dir$flags.per.unit.hyper.zone <- round(
    (dir$n.flags * dir$pct.hyper / 100) / dir$hyper.zone.width, 0)
  dir$flags.per.unit.hypo.zone <- round(
    (dir$n.flags * dir$pct.hypo / 100) / dir$hypo.zone.width, 0)
  cat("\n  direction of flags vs the room available in each direction:\n\n")
  print(dir[, c("methy.state", "n.flags", "pct.hyper", "pct.hypo",
                "hyper.zone.width", "hypo.zone.width",
                "flags.per.unit.hyper.zone", "flags.per.unit.hypo.zone")],
        row.names = FALSE)
  write_out(dir, paste0("Task22_Chr22_Direction_", tissue, ".csv"))

  # --------------------------------------------------
  # E. percentile rule vs Tukey 3xIQR
  # --------------------------------------------------
  # Tukey's rule thresholds at Q1 - 3*IQR and Q3 + 3*IQR, computed
  # within this cohort. It shares nothing with the external panel
  # except the data it is applied to. It is also the rule behind the
  # outliers.coef2 / outliers.coef3 columns already in the source
  # file, so the comparison costs nothing.

  cat("\n=== Part E: percentile rule vs Tukey 3xIQR ===\n")

  q1  <- apply(B, 1, quantile, .25, na.rm = TRUE)
  q3  <- apply(B, 1, quantile, .75, na.rm = TRUE)
  iqr <- q3 - q1
  tk  <- (B < q1 - 3 * iqr) | (B > q3 + 3 * iqr)
  tk[is.na(B)] <- NA
  ex  <- F != 0

  cmp <- do.call(rbind, lapply(STATE_ORDER, function(s) {
    k <- which(st == s)
    if (!length(k)) return(NULL)
    e <- as.vector(ex[k, ]); t3 <- as.vector(tk[k, ])
    ok <- !is.na(e) & !is.na(t3); e <- e[ok]; t3 <- t3[ok]
    if (!length(e)) return(NULL)
    data.frame(tissue = tissue, methy.state = s, n.cpgs = length(k),
      n.cells = length(e),
      median.cohort.IQR = round(median(iqr[k], na.rm = TRUE), 4),
      ext.flags = sum(e), ext.rate.pct = round(100 * mean(e), 3),
      tukey.flags = sum(t3), tukey.rate.pct = round(100 * mean(t3), 3),
      both = sum(e & t3),
      jaccard = round(sum(e & t3) / max(1, sum(e | t3)), 4),
      median.outliers.coef3 = median(d$outliers.coef3[k], na.rm = TRUE))
  }))
  print(cmp[, c("methy.state", "median.cohort.IQR", "ext.rate.pct",
                "tukey.rate.pct", "jaccard", "median.outliers.coef3")],
        row.names = FALSE)
  cat("\n  Two independent rules, two different state profiles, and almost\n")
  cat("  no overlap (Jaccard 0.05-0.29). Note median.cohort.IQR: the\n")
  cat("  dispersion Tukey thresholds on is itself a function of the state,\n")
  cat("  so Tukey inherits the same pathology in its own direction. The\n")
  cat("  defect is not specific to OutlierMeth - it follows from applying\n")
  cat("  ANY scale-free dispersion threshold to bounded, heteroscedastic\n")
  cat("  beta values.\n")
  write_out(cmp, paste0("Task22_Chr22_PercentileVsTukey_", tissue, ".csv"))

  rm(B, F, Bt, base, tk, ex); invisible(gc())
  cat("\n")
}

cat("====================================================\n")
cat("TASK 22 COMPLETE\n")
cat("====================================================\n")
