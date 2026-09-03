# ================================================================
# Task 29 - Where do the contrary flags actually come from?
# Date: Sep 3, 2026
#
# THE QUESTION (2026-09-03 email, "Analysis to-do" item 1)
#
#   For the normal 100-CpG window on chr22:
#     L  sites - 8  "-1" flags from the ext reference
#     LM sites - 21 "-1"
#     H  sites - 10 "+1"   (the email says 11; see NOTE below)
#     HM sites - 46 "+1"
#   and the same for tumour.
#
#   "to see if they are from all CG sites, or from 1 or 2 CG sites.
#    to see if they are from 1 or 2-3 samples, or if every sample
#    can have a chance."
#
# That is a concentration question, and it decides how the contrary
# flags should be read:
#
#   * concentrated on FEW SITES -> the sites are the problem. A
#     handful of badly behaved probes (cross-reactive, SNP-under-
#     probe, poorly mapped) manufacture the whole signal, and probe
#     masking removes it.
#   * concentrated in FEW SAMPLES -> the samples are the problem.
#     A few arrays carry global burden (batch, purity, conversion),
#     which is the McCartney et al. confound.
#   * SPREAD over both -> neither. The flags are a property of the
#     threshold rule itself, which is what Task 22's geometry
#     argument predicts, and no amount of QC will remove them.
#
# The three readings imply three different fixes, so this has to be
# measured rather than assumed.
#
# ---------------------------------------------------------------
# NOTE ON THE H COUNT
#
# The email lists 11 "+1" at H sites. The window has 11 H SITES and
# 10 "+1" flags across them (Task27_Chr22_Full100CpG_SiteTable.csv,
# sum of ext.pos1.N over state.normal == "H"). This script recounts
# from the flag matrix and prints its own totals; the 10 is what the
# data gives. Every other number in the email reproduces exactly.
#
# ---------------------------------------------------------------
# WHAT IS MEASURED
#
# For each tissue x state x direction cell, over the evaluable
# (non-NA) cells only:
#
#   sites   - how many sites of that state exist, how many carry at
#             least one flag, the largest single-site count, the
#             share held by the top 1 and top 2 sites, and Gini.
#   samples - the same over the 53 samples.
#   p.conc  - Monte Carlo p for the observed top-1 share under a
#             uniform multinomial allocation of the same number of
#             flags over the same number of bins. Small p = genuinely
#             concentrated; large p = as spread as chance allows.
#   magnitude - median |beta - cohort median| of the flagged cells,
#             and the share below 0.05 / 0.10, so the concentration
#             answer is read next to the effect size.
#
# A per-flag detail table is written as well: every contrary flag
# with its site, sample, beta, cohort median and deviation.
#
# PRIVACY: the detail table carries per-sample beta values for
# flagged cells only, under de-identified column labels (N1..N53),
# exactly as Task18/Task20's ContraryDetail files already do. It
# stays out of git for the same reason - see CLASSIFICATION.md.
# The summary table is counts only and is safe to commit.
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

N_SITES <- 100
STATE_ORDER <- c("L", "LM", "M", "HM", "H", "R")
# the four cells the email asks about, plus their mirror images so the
# "contrary" numbers can be read against the concordant ones
CELLS <- list(
  list(state = "L",  dir = -1L, label = "contrary"),
  list(state = "LM", dir = -1L, label = "contrary"),
  list(state = "H",  dir =  1L, label = "contrary"),
  list(state = "HM", dir =  1L, label = "contrary"),
  list(state = "L",  dir =  1L, label = "concordant"),
  list(state = "LM", dir =  1L, label = "concordant"),
  list(state = "H",  dir = -1L, label = "concordant"),
  list(state = "HM", dir = -1L, label = "concordant"),
  list(state = "M",  dir = -1L, label = "neutral"),
  list(state = "M",  dir =  1L, label = "neutral"),
  list(state = "R",  dir = -1L, label = "neutral"),
  list(state = "R",  dir =  1L, label = "neutral"))
N_MC <- 20000

# Emitted direction values are "hyper"/"hypo", not "+1"/"-1". A CSV column
# holding only "+1" and "-1" is read back by read.csv as NUMERIC 1 and -1,
# so a downstream filter on direction == "+1" silently matches nothing.
# The tables in results.md keep the +1 / -1 notation; only the machine-read
# column changes.

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

# Gini on a count vector, 0 = perfectly even over the bins supplied.
gini <- function(x) {
  if (!length(x) || sum(x) == 0) return(NA_real_)
  x <- sort(as.numeric(x)); n <- length(x)
  sum((2 * seq_len(n) - n - 1) * x) / (n * sum(x))
}

# Monte Carlo p for "the top bin holds at least this many", when k
# flags are thrown uniformly and independently into b bins. This is
# the right null for "could every site / every sample have had an
# equal chance?" - it is exactly the email's phrasing.
p_top1 <- function(k, b, observed_max, reps = N_MC) {
  if (k < 1 || b < 2) return(NA_real_)
  hits <- 0L
  for (i in seq_len(reps)) {
    tb <- tabulate(sample.int(b, k, replace = TRUE), nbins = b)
    if (max(tb) >= observed_max) hits <- hits + 1L
  }
  (hits + 1) / (reps + 1)
}

cat("=========================================================\n")
cat("TASK 29 - CONTRARY-FLAG DECOMPOSITION, 100-CpG WINDOW\n")
cat("=========================================================\n\n")

summary_rows <- list(); site_rows <- list(); sample_rows <- list()
detail_rows  <- list()

for (tissue in c("Normal", "Tumor")) {

  cat("#########################################################\n")
  cat("### ", tissue, "\n", sep = "")
  cat("#########################################################\n\n")

  cols <- paste0(substr(tissue, 1, 1), 1:53)

  # the window is defined by the NORMAL file's chr22 ordering, which
  # is what the email pasted; the same 100 cgIDs are then looked up in
  # the tumour file so both tissues describe the same sites
  nrm <- read_any(inputs$Normal, tab = TRUE)
  n22 <- nrm[nrm$Chromosome == "chr22", ]; n22 <- n22[order(n22$Start), ]
  ids <- n22$Composite.Element.REF[seq_len(N_SITES)]
  pos <- n22$Start[seq_len(N_SITES)]
  rm(nrm, n22); invisible(gc())

  d <- read_any(inputs[[tissue]], tab = TRUE)
  d <- d[match(ids, d$Composite.Element.REF), ]
  B <- as.matrix(d[, cols, drop = FALSE]); storage.mode(B) <- "numeric"
  rownames(B) <- ids; colnames(B) <- cols
  st <- d$methy.state
  stopifnot(all(st %in% STATE_ORDER))
  rm(d); invisible(gc())

  fx <- read_any(flag_ext[[tissue]])
  F <- as.matrix(fx[match(ids, fx$cgID), cols, drop = FALSE])
  storage.mode(F) <- "numeric"; rownames(F) <- ids; colnames(F) <- cols
  rm(fx); invisible(gc())

  med <- apply(B, 1, median, na.rm = TRUE)
  dev <- abs(B - med)

  cat("sites by state:\n"); print(table(factor(st, levels = STATE_ORDER)))
  cat("\next flags by state and direction (this is the email's table):\n")
  tab <- t(sapply(STATE_ORDER, function(s) {
    k <- which(st == s)
    if (!length(k)) return(c(sites = 0, evaluable = 0, neg1 = 0, pos1 = 0))
    c(sites = length(k), evaluable = sum(!is.na(F[k, , drop = FALSE])),
      neg1 = sum(F[k, , drop = FALSE] == -1, na.rm = TRUE),
      pos1 = sum(F[k, , drop = FALSE] ==  1, na.rm = TRUE)) }))
  print(tab)
  cat("\n")

  for (cell in CELLS) {
    s <- cell$state; dr <- cell$dir
    k <- which(st == s)
    if (!length(k)) next
    Fs <- F[k, , drop = FALSE]; Bs <- B[k, , drop = FALSE]
    Ds <- dev[k, , drop = FALSE]
    hit <- which(!is.na(Fs) & Fs == dr, arr.ind = TRUE)
    nflag <- nrow(hit)

    per_site   <- tabulate(hit[, 1], nbins = length(k))
    per_sample <- tabulate(hit[, 2], nbins = ncol(Fs))
    ssort <- sort(per_site, decreasing = TRUE)
    psort <- sort(per_sample, decreasing = TRUE)
    shr <- function(v, n) if (nflag == 0) NA_real_ else
      round(100 * sum(head(v, n)) / nflag, 1)

    mags <- if (nflag) Ds[hit] else numeric(0)
    mags <- mags[!is.na(mags)]

    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      tissue = tissue, methy.state = s,
      direction = if (dr > 0) "hyper" else "hypo", cell = cell$label,
      sites.in.state = length(k),
      evaluable.cells = sum(!is.na(Fs)),
      flags = nflag,
      flag.rate.pct = round(100 * nflag / max(1, sum(!is.na(Fs))), 3),
      # --- concentration over SITES ---
      sites.with.flag = sum(per_site > 0),
      max.per.site = if (nflag) max(per_site) else 0L,
      top1.site.share.pct = shr(ssort, 1),
      top2.site.share.pct = shr(ssort, 2),
      gini.sites = round(gini(per_site), 3),
      p.site.concentration = if (nflag > 1)
        round(p_top1(nflag, length(k), max(per_site)), 4) else NA_real_,
      # --- concentration over SAMPLES ---
      samples.with.flag = sum(per_sample > 0),
      max.per.sample = if (nflag) max(per_sample) else 0L,
      top1.sample.share.pct = shr(psort, 1),
      top3.sample.share.pct = shr(psort, 3),
      gini.samples = round(gini(per_sample), 3),
      p.sample.concentration = if (nflag > 1)
        round(p_top1(nflag, ncol(Fs), max(per_sample)), 4) else NA_real_,
      # --- magnitude, so concentration is read next to effect size ---
      median.abs.dbeta = if (length(mags)) round(median(mags), 4) else NA_real_,
      pct.under.0.05 = if (length(mags)) round(100 * mean(mags < 0.05), 1) else NA_real_,
      pct.under.0.10 = if (length(mags)) round(100 * mean(mags < 0.10), 1) else NA_real_)

    if (cell$label == "contrary" && nflag > 0) {
      cat("--- ", tissue, " ", s, " ", if (dr > 0) "+1" else "-1",
          ": ", nflag, " flags over ", length(k), " sites x 53 samples\n", sep = "")
      cat("      sites carrying >=1: ", sum(per_site > 0), "/", length(k),
          "   max at one site: ", max(per_site),
          "   top-2 sites hold ", shr(ssort, 2), "%\n", sep = "")
      cat("      samples carrying >=1: ", sum(per_sample > 0), "/53",
          "   max in one sample: ", max(per_sample),
          "   top-3 samples hold ", shr(psort, 3), "%\n", sep = "")
      cat("      median |dbeta| ", round(median(mags), 4),
          "   under 0.10: ", round(100 * mean(mags < 0.10), 1), "%\n", sep = "")

      site_rows[[length(site_rows) + 1]] <- data.frame(
        tissue = tissue, methy.state = s,
        direction = if (dr > 0) "hyper" else "hypo",
        cgID = rownames(Fs), pos = pos[k],
        cohort.median.beta = round(med[k], 4),
        contrary.flags = per_site,
        pct.of.cell = round(100 * per_site / nflag, 1))
      sample_rows[[length(sample_rows) + 1]] <- data.frame(
        tissue = tissue, methy.state = s,
        direction = if (dr > 0) "hyper" else "hypo",
        sample = colnames(Fs), contrary.flags = per_sample,
        pct.of.cell = round(100 * per_sample / nflag, 1))
      detail_rows[[length(detail_rows) + 1]] <- data.frame(
        tissue = tissue, methy.state = s,
        direction = if (dr > 0) "hyper" else "hypo",
        cgID = rownames(Fs)[hit[, 1]], pos = pos[k][hit[, 1]],
        sample = colnames(Fs)[hit[, 2]],
        beta = round(Bs[hit], 4),
        cohort.median.beta = round(med[k][hit[, 1]], 4),
        abs.dbeta = round(Ds[hit], 4))
    }
  }
  cat("\n")

  # ---- how many samples carry the union of all four contrary cells ----
  cc <- do.call(rbind, Filter(function(z) z$tissue == tissue,
                              sample_rows))
  if (!is.null(cc)) {
    agg <- aggregate(contrary.flags ~ sample, cc, sum)
    agg <- agg[order(-agg$contrary.flags), ]
    tot <- sum(agg$contrary.flags)
    cat("all four contrary cells pooled: ", tot, " flags over ",
        sum(agg$contrary.flags > 0), " of 53 samples\n", sep = "")
    cat("  top 5 samples: ",
        paste(sprintf("%s=%d", head(agg$sample, 5), head(agg$contrary.flags, 5)),
              collapse = "  "), "\n", sep = "")
    cat("  top 3 hold ", round(100 * sum(head(agg$contrary.flags, 3)) / tot, 1),
        "%  (uniform expectation ", round(100 * 3 / 53, 1), "%)\n", sep = "")
    cat("  Gini over samples: ", round(gini(agg$contrary.flags), 3), "\n\n", sep = "")
  }

  rm(B, F, dev); invisible(gc())
}

sm <- do.call(rbind, summary_rows); rownames(sm) <- NULL
write_out(sm, "Task29_Window100_ContraryDecomposition.csv")
write_out(do.call(rbind, site_rows),   "Task29_Window100_ContraryBySite.csv")
write_out(do.call(rbind, sample_rows), "Task29_Window100_ContraryBySample.csv")
write_out(do.call(rbind, detail_rows), "Task29_Window100_ContraryDetail.csv")

cat("\n---- the four cells the email asked about ----\n")
show <- sm[sm$cell == "contrary",
  c("tissue","methy.state","direction","sites.in.state","flags",
    "sites.with.flag","max.per.site","top2.site.share.pct","p.site.concentration",
    "samples.with.flag","max.per.sample","top3.sample.share.pct","p.sample.concentration",
    "median.abs.dbeta","pct.under.0.10")]
print(show, row.names = FALSE)

cat("\n=========================================================\n")
cat("TASK 29 COMPLETE\n")
cat("=========================================================\n")

# ================================================================
# PART B - why THOSE sites?
#
# Part A shows the contrary flags concentrate on a minority of sites
# (4 of 8 L sites, 8 of 14 LM, 4 of 11 H) while spreading across
# samples. The obvious next question is whether the concentrating
# sites are biologically special or merely geometrically unlucky.
#
# flagMeth sets +1 iff beta > P and -1 iff beta < N, where P and N
# come from the external panel and know nothing about this cohort.
# So the number of samples a site flags is decided by ONE quantity:
# how far the external threshold sits from this cohort's own
# distribution at that site. Task 22 showed the thresholds can be
# recovered from the flag matrix itself - the unflagged cells bracket
# them - and that is reused here.
#
# If contrary-flag count is explained by the median-to-threshold
# distance measured in cohort dispersion units, the concentrating
# sites carry no biology: they are the sites where an external
# threshold happened to land inside this cohort's spread. If it is
# NOT explained, the sites are worth looking at individually.
#
# Reported per site: reconstructed P and N, the distance from the
# cohort median to each in raw beta and in cohort-MAD units, and the
# contrary-flag count. Spearman correlation is taken within state, so
# it cannot be driven by the between-state differences Task 22
# already documented.
# ================================================================

cat("\n\n=========================================================\n")
cat("TASK 29 PART B - IS SITE CONCENTRATION JUST THRESHOLD DISTANCE?\n")
cat("=========================================================\n\n")

geom_rows <- list()

for (tissue in c("Normal", "Tumor")) {
  cols <- paste0(substr(tissue, 1, 1), 1:53)

  nrm <- read_any(inputs$Normal, tab = TRUE)
  n22 <- nrm[nrm$Chromosome == "chr22", ]; n22 <- n22[order(n22$Start), ]
  ids <- n22$Composite.Element.REF[seq_len(N_SITES)]
  pos <- n22$Start[seq_len(N_SITES)]
  rm(nrm, n22); invisible(gc())

  d <- read_any(inputs[[tissue]], tab = TRUE)
  d <- d[match(ids, d$Composite.Element.REF), ]
  B <- as.matrix(d[, cols, drop = FALSE]); storage.mode(B) <- "numeric"
  rownames(B) <- ids
  st <- d$methy.state; rm(d); invisible(gc())

  fx <- read_any(flag_ext[[tissue]])
  F <- as.matrix(fx[match(ids, fx$cgID), cols, drop = FALSE])
  storage.mode(F) <- "numeric"; rm(fx); invisible(gc())

  for (i in seq_len(N_SITES)) {
    b <- B[i, ]; fl <- F[i, ]
    if (all(is.na(fl))) next
    up <- which(fl == 1); dn <- which(fl == -1); ze <- which(fl == 0)
    P <- NA_real_; N <- NA_real_
    # bracket P: highest unflagged/hypo beta <= P < lowest hyper beta
    lowP <- if (length(c(ze, dn))) max(b[c(ze, dn)], na.rm = TRUE) else NA_real_
    hiP  <- if (length(up)) min(b[up], na.rm = TRUE) else NA_real_
    if (is.finite(lowP) && is.finite(hiP)) P <- (lowP + hiP) / 2
    lowN <- if (length(dn)) max(b[dn], na.rm = TRUE) else NA_real_
    hiN  <- if (length(c(ze, up))) min(b[c(ze, up)], na.rm = TRUE) else NA_real_
    if (is.finite(lowN) && is.finite(hiN)) N <- (lowN + hiN) / 2

    m <- median(b, na.rm = TRUE); s <- mad(b, na.rm = TRUE)
    contrary_dir <- if (st[i] %in% c("L", "LM")) -1L else
                    if (st[i] %in% c("H", "HM"))  1L else NA_integer_
    geom_rows[[length(geom_rows) + 1]] <- data.frame(
      tissue = tissue, cgID = ids[i], pos = pos[i], methy.state = st[i],
      cohort.median = round(m, 4), cohort.mad = round(s, 4),
      thr.hyper.P = round(P, 4), thr.hypo.N = round(N, 4),
      dist.median.to.P = round(P - m, 4),
      dist.N.to.median = round(m - N, 4),
      mad.units.to.P = if (is.finite(s) && s > 0) round((P - m) / s, 2) else NA_real_,
      mad.units.to.N = if (is.finite(s) && s > 0) round((m - N) / s, 2) else NA_real_,
      n.pos1 = sum(fl == 1, na.rm = TRUE),
      n.neg1 = sum(fl == -1, na.rm = TRUE),
      contrary.dir = contrary_dir,
      contrary.flags = if (is.na(contrary_dir)) NA_integer_ else
        sum(fl == contrary_dir, na.rm = TRUE))
  }
  rm(B, F); invisible(gc())
}

gm <- do.call(rbind, geom_rows); rownames(gm) <- NULL
write_out(gm, "Task29_Window100_SiteThresholdGeometry.csv")

cat("Spearman rho between a site's contrary-flag count and how close the\n")
cat("external threshold sits to this cohort's median (in cohort-MAD units).\n")
cat("A strong NEGATIVE rho means: the more sites flag, the closer the\n")
cat("threshold - i.e. geometry, not biology.\n\n")

for (tissue in c("Normal", "Tumor")) {
  cat("--- ", tissue, "\n", sep = "")
  for (s in c("L", "LM", "HM", "H")) {
    g <- gm[gm$tissue == tissue & gm$methy.state == s, ]
    if (nrow(g) < 4) { cat(sprintf("    %-3s  n=%d  (too few sites)\n", s, nrow(g))); next }
    dist <- if (s %in% c("L", "LM")) g$mad.units.to.N else g$mad.units.to.P
    ok <- is.finite(dist) & is.finite(g$contrary.flags)
    if (sum(ok) < 4 || length(unique(g$contrary.flags[ok])) < 2) {
      cat(sprintf("    %-3s  n=%d  (no variation)\n", s, sum(ok))); next }
    rho <- suppressWarnings(cor(dist[ok], g$contrary.flags[ok], method = "spearman"))
    pv  <- suppressWarnings(cor.test(dist[ok], g$contrary.flags[ok],
                                     method = "spearman", exact = FALSE)$p.value)
    cat(sprintf("    %-3s  sites=%2d  contrary=%2d  rho=%+.3f  p=%.4f   median dist=%.2f MAD\n",
                s, sum(ok), sum(g$contrary.flags[ok]), rho, pv,
                median(dist[ok], na.rm = TRUE)))
  }
  cat("\n")
}

cat("the ten sites carrying the most contrary flags:\n")
top <- gm[!is.na(gm$contrary.flags) & gm$contrary.flags > 0, ]
top <- top[order(-top$contrary.flags), ]
print(head(top[, c("tissue","cgID","pos","methy.state","cohort.median",
                   "cohort.mad","thr.hypo.N","thr.hyper.P",
                   "mad.units.to.N","mad.units.to.P","contrary.flags")], 10),
      row.names = FALSE)

cat("\n=========================================================\n")
cat("TASK 29 PART B COMPLETE\n")
cat("=========================================================\n")
