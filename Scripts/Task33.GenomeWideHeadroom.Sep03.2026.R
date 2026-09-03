# ================================================================
# Task 33 - Tasks 29-32 taken off chromosome 22
# Date: Sep 3, 2026
#
# Three of this project's load-bearing results are still chr22-only:
# the threshold geometry and stability of results.md section 7, and
# now the effect-size ceiling and p-envelope of Task 31. plan.md steps
# 3 and 7b both ask for them genome-wide, and the Task 13 flag
# matrices already cover all 380,355 CpGs for both references and both
# tissues, so this is a wider denominator on existing logic rather
# than new flagging.
#
# It also closes the half of plan.md step 4 that was never run: does
# the per-sample burden ranking survive a magnitude floor? That is the
# question that decides whether N15 / N17 / N14 / N48 carry real
# burden or heteroscedasticity.
#
# ---------------------------------------------------------------
# PART A  effect-size ceiling per state (Task 31 Part A, genome-wide)
# PART B  what any setting of p could give (Task 31 Part B)
# PART C  threshold geometry and directional asymmetry (section 7)
# PART D  flag stability under measured-range noise (section 7)
# PART E  per-sample burden with and without the 0.10 floor
#
# ---------------------------------------------------------------
# COST
# 380,355 x 53 doubles is ~161 MB per matrix and two are held at once,
# so tissues are processed one at a time with an explicit gc() between
# them. Part D uses 5 replicates rather than the 10-20 used on chr22;
# at 380k CpGs the Monte Carlo error on a Jaccard is already small and
# the run would otherwise take over an hour.
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
R_GRID   <- c(NA, 3L, 1L)
NOISE_SD <- c(0.005, 0.010, 0.020)
N_REPS   <- 5
FLOOR    <- 0.10

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
jaccard <- function(a, b) {
  ok <- !is.na(a) & !is.na(b)
  A <- a[ok] != 0; B2 <- b[ok] != 0
  u <- sum(A | B2); if (!u) return(NA_real_)
  sum(A & B2 & a[ok] == b[ok]) / u
}

cat("=========================================================\n")
cat("TASK 33 - GENOME-WIDE HEADROOM, GEOMETRY, STABILITY, BURDEN\n")
cat("=========================================================\n\n")

ceil_rows <- list(); env_rows <- list(); geom_rows <- list()
stab_rows <- list(); burden_rows <- list()

for (tissue in c("Normal", "Tumor")) {

  cat("#########################################################\n")
  cat("### ", tissue, "   (all 380,355 CpGs x 53 samples)\n", sep = "")
  cat("#########################################################\n\n")

  cols <- paste0(substr(tissue, 1, 1), 1:53)
  d <- read_any(inputs[[tissue]], tab = TRUE)
  ids <- d$Composite.Element.REF
  B <- as.matrix(d[, cols, drop = FALSE]); storage.mode(B) <- "numeric"
  colnames(B) <- cols
  st <- d$methy.state
  stopifnot(all(st %in% STATE_ORDER))
  rm(d); invisible(gc())

  fx <- read_any(flag_ext[[tissue]])
  F <- as.matrix(fx[match(ids, fx$cgID), cols, drop = FALSE])
  storage.mode(F) <- "numeric"; rm(fx); invisible(gc())
  cat("  loaded ", nrow(B), " CpGs; ", sum(!is.na(F)),
      " evaluable cells (", round(100 * mean(!is.na(F)), 2), "%)\n\n", sep = "")

  med <- apply(B, 1, median, na.rm = TRUE)
  dev <- abs(B - med)

  # ---------------- PART A ----------------------------------
  cat("PART A - effect-size ceiling per state, genome-wide\n\n")
  for (s in STATE_ORDER) {
    k <- which(st == s); if (!length(k)) next
    ru <- 1 - med[k]; rd <- med[k]; dv <- dev[k, , drop = FALSE]
    row <- data.frame(
      tissue = tissue, methy.state = s, sites = length(k),
      median.cohort.beta = round(median(med[k], na.rm = TRUE), 4),
      room.up = round(median(ru, na.rm = TRUE), 4),
      room.down = round(median(rd, na.rm = TRUE), 4),
      pct.sites.roomup.under.0.10 = round(100 * mean(ru < 0.10, na.rm = TRUE), 1),
      pct.sites.roomdown.under.0.10 = round(100 * mean(rd < 0.10, na.rm = TRUE), 1),
      cells.over.0.05 = round(100 * mean(dv >= 0.05, na.rm = TRUE), 3),
      cells.over.0.10 = round(100 * mean(dv >= 0.10, na.rm = TRUE), 3),
      cells.over.0.15 = round(100 * mean(dv >= 0.15, na.rm = TRUE), 3))
    ceil_rows[[length(ceil_rows) + 1]] <- row
    cat(sprintf("   %-3s sites=%6d  median beta=%.3f  room up=%.3f down=%.3f  | cells >=0.05: %6.3f%%  >=0.10: %6.3f%%  >=0.15: %6.3f%%\n",
        s, length(k), row$median.cohort.beta, row$room.up, row$room.down,
        row$cells.over.0.05, row$cells.over.0.10, row$cells.over.0.15))
    rm(dv); invisible(gc())
  }

  # ---------------- PART B ----------------------------------
  # Same nesting argument as Task 31: any stricter p keeps a prefix of
  # each CpG's flags ordered by beta. Vectorised here rather than
  # looped, because 380k row-wise order() calls is not viable.
  cat("\nPART B - the p envelope, genome-wide\n")
  cat("  r = flags kept per CpG per direction; r=1 is the strictest reachable p.\n\n")

  keep_top <- function(r) {
    if (is.na(r)) return(F)
    O <- matrix(0L, nrow(F), ncol(F))
    O[is.na(F)] <- NA_integer_
    for (dir in c(1L, -1L)) {
      # rank each flagged cell within its CpG: hyper by descending beta,
      # hypo by ascending. A stricter threshold keeps the top-ranked r.
      M <- ifelse(!is.na(F) & F == dir, B, NA_real_)
      if (dir == -1L) M <- -M
      # rank via row-wise ordering, done in one pass over columns
      R <- matrix(NA_integer_, nrow(M), ncol(M))
      ord <- t(apply(M, 1, function(v) order(v, decreasing = TRUE, na.last = TRUE)))
      for (j in seq_len(ncol(M))) R[cbind(seq_len(nrow(M)), ord[, j])] <- j
      O[!is.na(M) & R <= r] <- dir
    }
    O
  }

  for (r in R_GRID) {
    O <- keep_top(r)
    lab <- if (is.na(r)) "all (p=0.01 as run)" else sprintf("top %d per CpG", r)
    cat("  r = ", lab, ":  ", sum(O != 0, na.rm = TRUE), " flags total\n", sep = "")
    for (s in STATE_ORDER) {
      k <- which(st == s); if (!length(k)) next
      Os <- O[k, , drop = FALSE]; dv <- dev[k, , drop = FALSE]
      # "hyper"/"hypo" rather than "+1"/"-1": a column holding only
      # "+1" and "-1" is read back by read.csv as NUMERIC 1 and -1, so
      # a downstream filter on direction == "+1" silently matches
      # nothing. Task 31's equivalent column escapes this only because
      # it also carries the value "both".
      for (dl in list(list(l = "hyper", m = 1L), list(l = "hypo", m = -1L))) {
        sel <- !is.na(Os) & Os == dl$m & !is.na(dv)
        n <- sum(sel); mg <- dv[sel]
        env_rows[[length(env_rows) + 1]] <- data.frame(
          tissue = tissue, keep.top.per.cpg = if (is.na(r)) NA_integer_ else r,
          setting = lab, methy.state = s, direction = dl$l, flags = n,
          median.abs.dbeta = if (n) round(median(mg), 4) else NA_real_,
          max.abs.dbeta = if (n) round(max(mg), 4) else NA_real_,
          pct.under.0.10 = if (n) round(100 * mean(mg < 0.10), 2) else NA_real_,
          surviving.floor.0.10 = if (n) sum(mg >= 0.10) else 0L)
        if (s %in% c("L", "H") )
          cat(sprintf("       %-3s %-5s n=%7d  median|db|=%.4f  under 0.10: %6.2f%%  surviving 0.10: %6d\n",
              s, dl$l, n, if (n) median(mg) else NA,
              if (n) 100*mean(mg < 0.10) else NA, if (n) sum(mg >= 0.10) else 0L))
      }
      rm(Os, dv)
    }
    rm(O); invisible(gc())
  }

  # ---------------- PART C ----------------------------------
  cat("\nPART C - threshold geometry and directional asymmetry, genome-wide\n")
  cat("  P, N recovered from the flag matrix (Task 22 Part A). flag.density is\n")
  cat("  flags per unit of available beta space in that direction.\n\n")

  nr <- nrow(B)
  P <- N <- rep(NA_real_, nr)
  for (i in seq_len(nr)) {
    fl <- F[i, ]; if (all(is.na(fl))) next
    b <- B[i, ]
    up <- which(fl == 1); dn <- which(fl == -1); ze <- which(fl == 0)
    if (length(up)) {
      lo <- if (length(c(ze, dn))) max(b[c(ze, dn)], na.rm = TRUE) else NA_real_
      hi <- min(b[up], na.rm = TRUE)
      if (is.finite(lo) && is.finite(hi)) P[i] <- (lo + hi) / 2
    }
    if (length(dn)) {
      lo <- max(b[dn], na.rm = TRUE)
      hi <- if (length(c(ze, up))) min(b[c(ze, up)], na.rm = TRUE) else NA_real_
      if (is.finite(lo) && is.finite(hi)) N[i] <- (lo + hi) / 2
    }
  }
  cat("  thresholds recovered: ", sum(is.finite(P)), " hyper, ",
      sum(is.finite(N)), " hypo, of ", nr, " CpGs\n\n", sep = "")

  for (s in STATE_ORDER) {
    k <- which(st == s); if (!length(k)) next
    npos <- sum(F[k, , drop = FALSE] ==  1, na.rm = TRUE)
    nneg <- sum(F[k, , drop = FALSE] == -1, na.rm = TRUE)
    zoneU <- median(1 - P[k], na.rm = TRUE)
    zoneD <- median(N[k], na.rm = TRUE)
    dU <- if (is.finite(zoneU) && zoneU > 0) npos / zoneU else NA_real_
    dD <- if (is.finite(zoneD) && zoneD > 0) nneg / zoneD else NA_real_
    # L and LM are pinned against 0, so the direction the state constrains is
    # hypo; H and HM are pinned against 1, so it is hyper. M and R are pinned
    # against neither and get no constrained direction.
    cdir <- if (s %in% c("L", "LM")) "hypo" else
            if (s %in% c("H", "HM")) "hyper" else NA_character_
    ratio_c <- if (is.na(cdir)) NA_real_ else
               if (cdir == "hypo") dD / dU else dU / dD
    geom_rows[[length(geom_rows) + 1]] <- data.frame(
      tissue = tissue, methy.state = s, sites = length(k),
      median.thr.hyper.P = round(median(P[k], na.rm = TRUE), 4),
      median.thr.hypo.N = round(median(N[k], na.rm = TRUE), 4),
      hyper.zone.1minusP = round(zoneU, 4),
      hypo.zone.N = round(zoneD, 4),
      zone.asymmetry.hyper.over.hypo = round(zoneU / zoneD, 3),
      flags.pos1 = npos, flags.neg1 = nneg,
      pct.flags.pos1 = round(100 * npos / max(1, npos + nneg), 1),
      density.hyper = round(dU, 0), density.hypo = round(dD, 0),
      constrained.direction = cdir,
      density.ratio.constrained.over.free = round(ratio_c, 1))
    cat(sprintf("   %-3s P=%.3f N=%.3f | hyper zone %.3f hypo zone %.3f | +1 %6d  -1 %6d  (%.1f%% +1) | constrained: %-5s  density %8.1fx higher there\n",
        s, median(P[k], na.rm = TRUE), median(N[k], na.rm = TRUE),
        zoneU, zoneD, npos, nneg,
        100 * npos / max(1, npos + nneg),
        if (is.na(cdir)) "-" else cdir, ratio_c))
  }

  # ---------------- PART D ----------------------------------
  cat("\nPART D - flag stability under noise, genome-wide\n")
  cat("  Noise levels bracket the range Task 32 measured from adjacent probes\n")
  cat("  (p10 ", if (tissue == "Normal") "0.0029" else "0.0036",
      ", median ", if (tissue == "Normal") "0.0142" else "0.0171", ").\n\n", sep = "")

  ext_apply <- function(X) {
    o <- matrix(0L, nrow(X), ncol(X), dimnames = dimnames(X))
    o[!is.na(P) & X > P] <-  1L
    o[!is.na(N) & X < N] <- -1L
    o[is.na(X)] <- NA; o
  }
  floor_apply <- function(Fm, X) {
    m <- apply(X, 1, median, na.rm = TRUE)
    o <- Fm; o[!is.na(Fm) & Fm != 0 & abs(X - m) < FLOOR] <- 0L; o
  }
  base_ext   <- ext_apply(B)
  base_floor <- floor_apply(base_ext, B)
  cat("  baseline flags: ext ", sum(base_ext != 0, na.rm = TRUE),
      "   ext+", FLOOR, " floor ", sum(base_floor != 0, na.rm = TRUE), "\n", sep = "")

  for (sd in NOISE_SD) {
    je <- jf <- 0; lost <- gain <- 0
    for (r in seq_len(N_REPS)) {
      Bn <- B + matrix(rnorm(length(B), 0, sd), nrow(B))
      Bn[Bn < 0] <- 0; Bn[Bn > 1] <- 1; Bn[is.na(B)] <- NA
      en <- ext_apply(Bn); fn <- floor_apply(en, Bn)
      je <- je + jaccard(base_ext, en)
      jf <- jf + jaccard(base_floor, fn)
      lost <- lost + sum(base_ext != 0 & en == 0, na.rm = TRUE)
      gain <- gain + sum(base_ext == 0 & en != 0, na.rm = TRUE)
      rm(Bn, en, fn); invisible(gc())
    }
    nb <- sum(base_ext != 0, na.rm = TRUE)
    stab_rows[[length(stab_rows) + 1]] <- data.frame(
      tissue = tissue, noise.sd = sd, reps = N_REPS,
      baseline.flags = nb,
      flags.lost = round(lost / N_REPS),
      flags.gained = round(gain / N_REPS),
      pct.baseline.lost = round(100 * (lost / N_REPS) / nb, 1),
      churn.over.baseline = round(((lost + gain) / N_REPS) / nb, 2),
      jaccard.ext = round(je / N_REPS, 4),
      jaccard.ext.floor = round(jf / N_REPS, 4))
    cat(sprintf("   sd=%.3f  lost %7.0f  gained %8.0f  %5.1f%% of baseline lost  churn %5.2f  |  Jaccard ext %.4f  ext+floor %.4f\n",
        sd, lost / N_REPS, gain / N_REPS, 100 * (lost / N_REPS) / nb,
        ((lost + gain) / N_REPS) / nb, je / N_REPS, jf / N_REPS))
  }

  # ---------------- PART E ----------------------------------
  # plan.md step 4, the half that was never run. If the per-sample
  # ranking survives a magnitude floor, the burden is real; if the
  # leaders collapse, it was heteroscedasticity.
  cat("\nPART E - per-sample burden with and without the ", FLOOR, " floor\n\n", sep = "")

  raw <- colSums(base_ext   != 0, na.rm = TRUE)
  flr <- colSums(base_floor != 0, na.rm = TRUE)
  rk_raw <- rank(-raw, ties.method = "min")
  rk_flr <- rank(-flr, ties.method = "min")
  rho <- cor(raw, flr, method = "spearman")
  burden_rows[[length(burden_rows) + 1]] <- data.frame(
    tissue = tissue, sample = names(raw),
    flags.raw = as.integer(raw), rank.raw = as.integer(rk_raw),
    flags.floored = as.integer(flr), rank.floored = as.integer(rk_flr),
    retained.pct = round(100 * flr / pmax(1, raw), 1),
    rank.change = as.integer(rk_raw - rk_flr))
  cat("  Spearman between raw and floored burden: ", round(rho, 4), "\n", sep = "")
  cat("  overall flags retained: ", round(100 * sum(flr) / sum(raw), 1), "%\n\n", sep = "")
  o <- order(rk_raw)
  cat("  top 8 by raw burden:\n")
  cat(sprintf("    %-5s %9s %6s %10s %6s %9s\n",
              "samp", "raw", "rank", "floored", "rank", "retained"))
  for (i in head(o, 8))
    cat(sprintf("    %-5s %9d %6d %10d %6d %8.1f%%\n",
        names(raw)[i], raw[i], rk_raw[i], flr[i], rk_flr[i],
        100 * flr[i] / max(1, raw[i])))
  cat("\n  top 8 by floored burden:\n")
  for (i in head(order(rk_flr), 8))
    cat(sprintf("    %-5s %9d %6d %10d %6d %8.1f%%\n",
        names(raw)[i], raw[i], rk_raw[i], flr[i], rk_flr[i],
        100 * flr[i] / max(1, raw[i])))

  rm(B, F, dev, base_ext, base_floor, P, N); invisible(gc())
  cat("\n")
}

write_out(do.call(rbind, ceil_rows),   "Task33_GenomeWide_EffectSizeCeiling.csv")
write_out(do.call(rbind, env_rows),    "Task33_GenomeWide_PEnvelope.csv")
write_out(do.call(rbind, geom_rows),   "Task33_GenomeWide_ThresholdGeometry.csv")
write_out(do.call(rbind, stab_rows),   "Task33_GenomeWide_FlagStability.csv")
write_out(do.call(rbind, burden_rows), "Task33_GenomeWide_SampleBurdenFloor.csv")

cat("\n=========================================================\n")
cat("TASK 33 COMPLETE\n")
cat("=========================================================\n")
