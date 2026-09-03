# ================================================================
# Task 31 - Can the OutlierMeth parameters be tuned out of trouble?
# Date: Sep 3, 2026
#
# THE QUESTION (2026-09-03 email, "Literature review" item 1)
#
#   "Check to see if you can adjust the 'parameter' setting to get
#    results. Currently, the 'ext' reference-based methods seem to
#    be very sensitive in identifying outliers."
#
# OutlierMeth exposes exactly two knobs (Packages/OutlierMeth/R):
#
#   reference  one of four packaged panels - tcga, geo, all, blood
#   p          one of four significance levels - 0.01, 0.001,
#              0.0001, 0.00001, which select the pre-computed
#              columns P<p> / N<p> in the panel
#
#   flagMeth(beta, reference, p):  +1 iff beta > P<p>,  -1 iff beta < N<p>
#
# There is nothing else. No magnitude argument, no per-sample call,
# no multiple-testing control (see Docs/LITERATURE.md).
#
# ---------------------------------------------------------------
# WHY THIS CAN BE ANSWERED WITHOUT tcga.rda
#
# tcga.rda is not available on this machine, so P0.001 cannot be
# read off the panel. It does not need to be. The four p levels are
# quantiles of one distribution, so
#
#     P0.01 <= P0.001 <= P0.0001 <= P0.00001
#
# and the flag set is therefore MONOTONE AND NESTED in p: lowering p
# can only raise the threshold, which can only remove flags, and the
# flags it removes are always the ones closest to the threshold.
# At each CpG the survivors under any stricter p are exactly a
# prefix of that CpG's currently-flagged cells sorted by beta.
#
# So the entire family of outcomes reachable by tuning p can be
# enumerated from the p = 0.01 flag matrix alone, by keeping the top
# r flags per CpG for r = 1, 2, 3, ... The r = 1 case is the
# strictest setting that still flags anything at all - a lower bound
# on sensitivity and an UPPER bound on achievable effect size.
#
# If the r = 1 magnitudes at H and L sites are still trivial, then no
# choice of p fixes the problem, and that is a statement about the
# method rather than about this particular panel.
#
# ---------------------------------------------------------------
# PART A - the ceiling
#   For each state: how much beta space is there between the cohort
#   median and the boundary the state pushes against? A hyper flag
#   at a site whose median is 0.96 can never exceed |dbeta| = 0.04,
#   whatever p is set to. This is the hard ceiling on effect size,
#   and it is a property of the state, not of the threshold.
#
# PART B - the p envelope
#   Sweep r = all, 5, 3, 2, 1 and report flags retained and the
#   magnitude distribution of the survivors, per state.
#
# PART C - deltMeth and relMeth
#   The package already computes beta - threshold (deltMeth) and
#   (beta - P)/(1 - P) (relMeth). deltMeth is one line away from the
#   magnitude floor this project has been arguing for. relMeth is
#   the opposite: it divides by the remaining head-room, which is
#   0.037 at H sites, so it INFLATES exactly the flags that should be
#   discounted. Both are quantified here because both ship in the
#   package and a user could reasonably reach for either.
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
R_GRID <- c(NA, 5L, 3L, 2L, 1L)      # NA = keep all (p = 0.01 as run)
FLOORS <- c(0.05, 0.10, 0.15, 0.20)

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

cat("=========================================================\n")
cat("TASK 31 - PARAMETER HEADROOM IN OutlierMeth\n")
cat("=========================================================\n\n")

ceil_rows <- list(); env_rows <- list(); rel_rows <- list()

for (tissue in c("Normal", "Tumor")) {

  cat("#########################################################\n")
  cat("### ", tissue, "   (chr22, 6,809 CpGs x 53 samples)\n", sep = "")
  cat("#########################################################\n\n")

  cols <- paste0(substr(tissue, 1, 1), 1:53)
  d <- read_any(inputs[[tissue]], tab = TRUE)
  d <- d[d$Chromosome == "chr22", ]; d <- d[order(d$Start), ]
  ids <- d$Composite.Element.REF
  B <- as.matrix(d[, cols, drop = FALSE]); storage.mode(B) <- "numeric"
  st <- d$methy.state; rm(d); invisible(gc())

  fx <- read_any(flag_ext[[tissue]])
  F <- as.matrix(fx[match(ids, fx$cgID), cols, drop = FALSE])
  storage.mode(F) <- "numeric"; rm(fx); invisible(gc())

  med <- apply(B, 1, median, na.rm = TRUE)
  dev <- abs(B - med)

  # ---------------- PART A : the ceiling --------------------
  # head-room above / below the cohort median at each site, and how
  # many cells could EVER clear a given floor there.
  cat("PART A - effect-size ceiling per state\n")
  cat("  room.up   = 1 - cohort median  (largest possible hyper |dbeta| at that site)\n")
  cat("  room.down =     cohort median  (largest possible hypo  |dbeta|)\n")
  cat("  cells.over.0.10 = share of ALL cells at that state already >= 0.10 from their median\n\n")
  for (s in STATE_ORDER) {
    k <- which(st == s); if (!length(k)) next
    ru <- 1 - med[k]; rd <- med[k]
    dv <- dev[k, , drop = FALSE]
    row <- data.frame(
      tissue = tissue, methy.state = s, sites = length(k),
      median.cohort.beta = round(median(med[k], na.rm = TRUE), 4),
      room.up = round(median(ru, na.rm = TRUE), 4),
      room.down = round(median(rd, na.rm = TRUE), 4),
      pct.sites.roomup.under.0.10 = round(100 * mean(ru < 0.10, na.rm = TRUE), 1),
      pct.sites.roomdown.under.0.10 = round(100 * mean(rd < 0.10, na.rm = TRUE), 1),
      cells.over.0.05 = round(100 * mean(dv >= 0.05, na.rm = TRUE), 2),
      cells.over.0.10 = round(100 * mean(dv >= 0.10, na.rm = TRUE), 2),
      cells.over.0.15 = round(100 * mean(dv >= 0.15, na.rm = TRUE), 2))
    ceil_rows[[length(ceil_rows) + 1]] <- row
    cat(sprintf("   %-3s sites=%5d  median beta=%.3f  room up=%.3f down=%.3f  | cells >=0.05: %5.2f%%  >=0.10: %5.2f%%  >=0.15: %5.2f%%\n",
        s, length(k), row$median.cohort.beta, row$room.up, row$room.down,
        row$cells.over.0.05, row$cells.over.0.10, row$cells.over.0.15))
  }

  # ---------------- PART B : the p envelope -----------------
  # keep the top r flags per CpG per direction; r = 1 is the
  # strictest reachable non-empty setting.
  cat("\nPART B - what any setting of p could give\n")
  cat("  r = number of flags kept per CpG per direction; r=1 is the strictest\n")
  cat("  non-empty threshold reachable by lowering p. Magnitudes are |beta - cohort median|.\n\n")

  keep_top <- function(F, B, r) {
    if (is.na(r)) return(F)
    O <- matrix(0L, nrow(F), ncol(F))
    for (i in seq_len(nrow(F))) {
      fl <- F[i, ]; if (all(is.na(fl))) { O[i, ] <- NA; next }
      b <- B[i, ]
      up <- which(!is.na(fl) & fl == 1)
      if (length(up)) { o <- up[order(b[up], decreasing = TRUE)]
                        O[i, head(o, r)] <- 1L }
      dn <- which(!is.na(fl) & fl == -1)
      if (length(dn)) { o <- dn[order(b[dn], decreasing = FALSE)]
                        O[i, head(o, r)] <- -1L }
      O[i, is.na(fl)] <- NA
    }
    O
  }

  for (r in R_GRID) {
    O <- keep_top(F, B, r)
    lab <- if (is.na(r)) "all (p=0.01 as run)" else sprintf("top %d per CpG", r)
    tot <- sum(O != 0, na.rm = TRUE)
    cat("  r = ", lab, ":  ", tot, " flags total\n", sep = "")
    for (s in STATE_ORDER) {
      k <- which(st == s); if (!length(k)) next
      Os <- O[k, , drop = FALSE]; dv <- dev[k, , drop = FALSE]
      # Split by direction. Pooling them is misleading at L and H: the
      # two directions have wildly different head-room there (at L,
      # 0.027 down and 0.973 up), so a change in the hyper/hypo MIX
      # moves the pooled median without any flag changing magnitude.
      for (dsel in list(list(lab = "both", m = Os != 0),
                        list(lab = "+1",   m = Os ==  1),
                        list(lab = "-1",   m = Os == -1))) {
        sel <- !is.na(Os) & dsel$m & !is.na(dv)
        n <- sum(sel); mg <- dv[sel]
        env_rows[[length(env_rows) + 1]] <- data.frame(
          tissue = tissue, keep.top.per.cpg = if (is.na(r)) NA_integer_ else r,
          setting = lab, methy.state = s, direction = dsel$lab, flags = n,
          flag.rate.pct = round(100 * n / max(1, sum(!is.na(Os))), 4),
          median.abs.dbeta = if (n) round(median(mg), 4) else NA_real_,
          p90.abs.dbeta = if (n) round(quantile(mg, .9, names = FALSE), 4) else NA_real_,
          max.abs.dbeta = if (n) round(max(mg), 4) else NA_real_,
          pct.under.0.05 = if (n) round(100 * mean(mg < 0.05), 1) else NA_real_,
          pct.under.0.10 = if (n) round(100 * mean(mg < 0.10), 1) else NA_real_,
          surviving.floor.0.10 = if (n) sum(mg >= 0.10) else 0L,
          surviving.floor.0.15 = if (n) sum(mg >= 0.15) else 0L)
        if (s %in% c("L", "H", "R") && dsel$lab != "both")
          cat(sprintf("       %-3s %-2s n=%6d  median|db|=%.4f  p90=%.4f  max=%.4f  under 0.10: %5.1f%%\n",
              s, dsel$lab, n, if (n) median(mg) else NA,
              if (n) quantile(mg,.9,names=FALSE) else NA,
              if (n) max(mg) else NA, if (n) 100*mean(mg < 0.10) else NA))
      }
    }
  }

  # ---------------- PART C : deltMeth vs relMeth ------------
  # reconstruct P and N from the flag matrix (Task 22 Part A) so the
  # package's own two magnitude measures can be computed here.
  cat("\nPART C - the package's own magnitude functions\n")
  P <- N <- rep(NA_real_, nrow(B))
  for (i in seq_len(nrow(B))) {
    b <- B[i, ]; fl <- F[i, ]; if (all(is.na(fl))) next
    up <- which(fl == 1); dn <- which(fl == -1); ze <- which(fl == 0)
    loP <- if (length(c(ze, dn))) max(b[c(ze, dn)], na.rm = TRUE) else NA_real_
    hiP <- if (length(up)) min(b[up], na.rm = TRUE) else NA_real_
    if (is.finite(loP) && is.finite(hiP)) P[i] <- (loP + hiP) / 2
    loN <- if (length(dn)) max(b[dn], na.rm = TRUE) else NA_real_
    hiN <- if (length(c(ze, up))) min(b[c(ze, up)], na.rm = TRUE) else NA_real_
    if (is.finite(loN) && is.finite(hiN)) N[i] <- (loN + hiN) / 2
  }
  cat("  thresholds recovered for ", sum(is.finite(P)), " CpGs (hyper) and ",
      sum(is.finite(N)), " (hypo) of ", nrow(B), "\n", sep = "")

  for (s in STATE_ORDER) {
    k <- which(st == s); if (!length(k)) next
    up <- which(!is.na(F[k, , drop = FALSE]) & F[k, , drop = FALSE] == 1, arr.ind = TRUE)
    if (!nrow(up)) next
    Bs <- B[k, , drop = FALSE]; Ps <- P[k]
    b  <- Bs[up]; pp <- Ps[up[, 1]]
    ok <- is.finite(b) & is.finite(pp)
    if (!any(ok)) next
    delt <- b[ok] - pp[ok]                    # deltMeth
    rel  <- delt / (1 - pp[ok])               # relMeth
    rel_rows[[length(rel_rows) + 1]] <- data.frame(
      tissue = tissue, methy.state = s, hyper.flags = sum(ok),
      median.threshold.P = round(median(pp[ok]), 4),
      median.headroom.1minusP = round(median(1 - pp[ok]), 4),
      median.deltMeth = round(median(delt), 4),
      median.relMeth = round(median(rel), 4),
      inflation.rel.over.delt = round(median(rel) / max(1e-9, median(delt)), 1),
      pct.deltMeth.under.0.05 = round(100 * mean(delt < 0.05), 1),
      pct.relMeth.under.0.05 = round(100 * mean(rel < 0.05), 1))
    cat(sprintf("   %-3s hyper n=%6d  P=%.3f  1-P=%.3f  deltMeth=%.4f  relMeth=%.4f  (relMeth inflates %.0fx)\n",
        s, sum(ok), median(pp[ok]), median(1 - pp[ok]),
        median(delt), median(rel), median(rel) / max(1e-9, median(delt))))
  }

  cat("\n")
  rm(B, F, dev); invisible(gc())
}

write_out(do.call(rbind, ceil_rows), "Task31_EffectSizeCeiling.csv")
write_out(do.call(rbind, env_rows),  "Task31_PEnvelope.csv")
write_out(do.call(rbind, rel_rows),  "Task31_DeltVsRelMeth.csv")

cat("\n---- bottom line: the strictest reachable p (r = 1) ----\n")
ev <- do.call(rbind, env_rows)
print(ev[!is.na(ev$keep.top.per.cpg) & ev$keep.top.per.cpg == 1 &
           ev$direction != "both",
         c("tissue","methy.state","direction","flags","median.abs.dbeta",
           "max.abs.dbeta","pct.under.0.10","surviving.floor.0.10")],
      row.names = FALSE)

cat("\n---- the same cells at p = 0.01 as actually run, for contrast ----\n")
print(ev[is.na(ev$keep.top.per.cpg) & ev$direction != "both",
         c("tissue","methy.state","direction","flags","median.abs.dbeta",
           "pct.under.0.10","surviving.floor.0.10")], row.names = FALSE)

cat("\n=========================================================\n")
cat("TASK 31 COMPLETE\n")
cat("=========================================================\n")

# ================================================================
# PART D - how many of the four p levels are actually distinct?
#
# This is arithmetic on the packaged panels, not on our data, and it
# turned out to be the single most important thing in Task 31.
#
# referenceMeth() builds every threshold with the plain R call
#   quantile(x, 1 - p, na.rm = TRUE)
# which is type 7: the value sits at position h = (n-1)(1-p) + 1 in
# the sorted reference, interpolating between order statistics
# floor(h) and floor(h)+1. So level p is an actual estimate of a
# tail only when n is large enough for h to fall clear of the top of
# the sample - roughly n >= 1/p + 2. Below that, the "threshold" is
# the reference MAXIMUM with a small interpolation, and successive p
# levels differ by a fraction of the gap between the two most
# extreme reference samples.
#
# This is the identical order-statistic degeneracy that README
# documents for the n = 53 self-reference arm. What was not
# previously noticed is that it also bites the EXTERNAL arm at the
# strict end, because none of the packaged panels is large enough:
#
#   tcga   747 samples, 21 tissue types   (this project's reference)
#   geo  1,268 samples, 15 tissue types
#   all  2,015 samples, 25 tissue types   (tcga + geo)
#   blood  656 samples, whole blood
#
# Sample counts are from Downs, Thursby & Cope 2023, Epigenetics
# 18(1):2213874. Note that 2,015 is the COMBINED "all" panel; this
# project loads tcga.rda and calls flagMeth(..., reference = tcga),
# so its reference is the 747-sample one.
# ================================================================

cat("\n\n=========================================================\n")
cat("TASK 31 PART D - ARE THE FOUR p LEVELS DISTINCT?\n")
cat("=========================================================\n\n")

PANELS <- c(tcga = 747, geo = 1268, all = 2015, blood = 656, `self(n=53)` = 53)
PLEVELS <- c(0.01, 0.001, 0.0001, 0.00001)

plev_rows <- list()
cat(sprintf("%-12s %6s %10s %12s %14s %16s %s\n", "panel", "n", "p",
            "position h", "interpolates", "rank from top", "expected n above"))
for (nm in names(PANELS)) {
  n <- PANELS[[nm]]
  for (p in PLEVELS) {
    h <- (n - 1) * (1 - p) + 1
    lo <- floor(h)
    rank_top <- n - lo
    usable <- n >= (1 / p) + 2
    plev_rows[[length(plev_rows) + 1]] <- data.frame(
      panel = nm, n.reference = n, p = p,
      position.h = round(h, 3),
      order.stats.spanned = sprintf("%d-%d", lo, min(lo + 1L, n)),
      rank.from.top = rank_top,
      expected.n.above = round(n * p, 3),
      distinct.from.max = usable)
    cat(sprintf("%-12s %6d %10g %12.2f %14s %16d %16.3f%s\n",
        nm, n, p, h, sprintf("%d-%d", lo, min(lo + 1L, n)),
        rank_top, n * p, if (usable) "" else "   <- pinned to the top"))
  }
}
pl <- do.call(rbind, plev_rows); rownames(pl) <- NULL
write_out(pl, "Task31_PLevelResolution.csv")

cat("\nsmallest reference n at which each p level clears the top two order statistics:\n")
for (p in PLEVELS) {
  n <- 3L; while (floor((n - 1) * (1 - p) + 1) > n - 2L) n <- n + 1L
  cat(sprintf("   p = %-8g  needs n >= %6d   (largest packaged panel is 2,015)\n", p, n))
}

cat("\nCONSEQUENCE\n")
cat("  With reference = tcga (n = 747), which is what this project ran:\n")
cat("    p = 0.01     -> the 8th largest of 747. A real tail estimate.\n")
cat("    p = 0.001    -> between the 746th and 747th, i.e. between the two\n")
cat("                    most extreme reference samples.\n")
cat("    p = 0.0001   -> 93 percent of the way from the 746th to the 747th.\n")
cat("    p = 0.00001  -> 99 percent of the way. Indistinguishable from p = 0.0001.\n")
cat("  So the p argument offers TWO usable settings on this panel, not four,\n")
cat("  and the strict end is the reference maximum however it is spelled.\n")
cat("  Even the largest packaged panel (all, n = 2,015) reaches only the 3rd\n")
cat("  largest at p = 0.001 and is pinned to the top at 0.0001 and below.\n")
cat("  Combined with Part B - where even the strictest reachable threshold\n")
cat("  leaves H-site hyper flags 100 percent below |dbeta| = 0.10 - the answer to\n")
cat("  \"can the parameters be adjusted to fix this\" is no.\n")

cat("\n=========================================================\n")
cat("TASK 31 PART D COMPLETE\n")
cat("=========================================================\n")
