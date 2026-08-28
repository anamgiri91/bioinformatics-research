# ================================================================
# Task 21 - Is a top-ranked sample a real outlier, or an artifact
#           of the window it was measured in?
# Date: Aug 28, 2026
#
# Continues Task 19. Task 19 established two things about N37, which
# tops the Task 17 per-sample summary with 11 of 71 biological flags:
#
#   (a) the total is STRUCTURALLY FIXED. At n = 53 the bio rule
#       flags exactly one sample per eligible CpG, so the summary
#       counts how often a sample is the cohort extreme, not how
#       many outliers it carries.
#   (b) N37's flags are SPATIALLY CLUSTERED. Ten of its eleven fall
#       on LM sites inside one ~9.2 kb window.
#
# Both are true but neither settles the question the debug was for:
# is N37 an unusual sample? This script answers that by widening the
# denominator, which is the one thing Task 19 did not do.
#
# THREE NESTED SCALES are compared for every sample and method:
#   window  - the 100 CpGs of Task 17/18/19       (index100)
#   tight   - the 100 CpGs of minimum span        (tight100, Task 20)
#   chr22   - all 6,809 chr22 CpGs
# and, for context, the genome-wide self-referential burden already
# computed in Task 15 is joined on.
#
# A sample that leads at one scale and sits mid-pack at the next is
# not an outlier sample; it won a 100-site lottery. A sample that
# leads at every scale is worth investigating.
#
# ---------------------------------------------------------------
# EVENT COLLAPSING
#
# Flags at neighbouring CpGs are not independent. Following the
# epimutation definition used by epimutacions (Barbosa et al.) - a
# run of at least 3 outlier CpGs no more than 1 kb apart - this
# script reports, per sample:
#   raw.flags     - cells flagged
#   events.1kb    - runs of flags merged at 1 kb
#   epimutations  - those runs containing >= 3 flagged CpGs
# The gap between raw.flags and events is the amount of double
# counting in the per-sample summary.
# ---------------------------------------------------------------
#
# Outputs:
#   Task21_Chr22_SampleBurden_<tissue>.csv     per sample x method x scale
#   Task21_Chr22_ScaleComparison_<tissue>.csv  rank of each sample per scale
#   Task21_Chr22_ContraryDrivers_<tissue>.csv  which samples drive the
#                                              state-contrary ext flags
#   Task21_Chr22_N37_Cluster_Normal.csv        the LM cluster, all 53 samples
# ================================================================

options(stringsAsFactors = FALSE)

pick_dir <- function(...) {
  cand <- c(...)
  for (d in cand) if (dir.exists(d)) return(d)
  stop("no candidate directory exists:\n  ", paste(cand, collapse = "\n  "))
}

src_dir <- pick_dir(
  "/home/s_s355/research3.data/TCGA/filter.BRCA.UCEC.Aug17.2022/sorted.TCGA.380355cg.files",
  path.expand("~/Downloads/tcga_brca_380355cg")
)
repo_dir <- pick_dir(
  "/mmfs1/home/wln26/Experiments.Outlier.July31.2026",
  path.expand("~/Desktop/bioinformatics-research")
)

input_normal <- file.path(src_dir, "Sorted.BRCA.53Alive.Normal.380355cg.75col.May28.2026.txt")
input_tumor  <- file.path(src_dir, "Sorted.BRCA.53Alive.Tumor.380355cg.75col.May28.2026.txt")
out_dir <- file.path(repo_dir, "Results")

flag_files <- list(
  Normal = list(self = file.path(out_dir, "Task13_Normal_selfref_flags.csv"),
                ext  = file.path(out_dir, "Task13_Normal_extref_flags.csv")),
  Tumor  = list(self = file.path(out_dir, "Task13_Tumor_selfref_flags.csv"),
                ext  = file.path(out_dir, "Task13_Tumor_extref_flags.csv"))
)
task15_files <- list(
  Normal = file.path(out_dir, "Task15_Normal_flagbias_by_sample.csv"),
  Tumor  = file.path(out_dir, "Task15_Tumor_flagbias_by_sample.csv")
)

N_SITES      <- 100
STATE_ORDER  <- c("L", "LM", "M", "HM", "H", "R")
P_LEVEL      <- 0.01
EVENT_GAP    <- 1000    # bp, epimutacions' window
MIN_RUN      <- 3       # CpGs, epimutacions' minimum run length
FOCUS_SAMPLE <- "N37"

has_dt <- requireNamespace("data.table", quietly = TRUE)

read_src <- function(p) {
  if (has_dt) {
    as.data.frame(data.table::fread(p, sep = "\t", header = TRUE,
                                    check.names = FALSE, showProgress = FALSE))
  } else {
    read.table(p, header = TRUE, sep = "\t", check.names = FALSE)
  }
}

read_big <- function(p) {
  if (has_dt) {
    as.data.frame(data.table::fread(p, showProgress = FALSE))
  } else {
    read.csv(p, check.names = FALSE)
  }
}

write_out <- function(tbl, fname) {
  path <- file.path(out_dir, fname)
  write.csv(tbl, path, row.names = FALSE)
  cat("  written:", path, " (", nrow(tbl), " rows )\n", sep = "")
  invisible(path)
}

build_bio_flag <- function(meth, states, leave_one_out = FALSE) {
  n_cpg <- nrow(meth); n_smp <- ncol(meth)
  out <- matrix(0, nrow = n_cpg, ncol = n_smp,
                dimnames = list(rownames(meth), colnames(meth)))
  for (i in seq_len(n_cpg)) {
    state <- states[i]; values <- as.numeric(meth[i, ])
    if (all(is.na(values))) { out[i, ] <- NA; next }
    if (!state %in% c("L", "LM", "H", "HM")) next
    hi <- state %in% c("L", "LM")
    p  <- if (hi) 1 - P_LEVEL else P_LEVEL
    if (!leave_one_out) {
      thr <- quantile(values, probs = p, na.rm = TRUE, names = FALSE)
      out[i, ] <- if (hi) as.integer(values > thr) else -as.integer(values < thr)
      out[i, is.na(values)] <- NA
    } else {
      for (j in seq_len(n_smp)) {
        if (is.na(values[j])) { out[i, j] <- NA; next }
        others <- values[-j]
        if (all(is.na(others))) { out[i, j] <- NA; next }
        thr <- quantile(others, probs = p, na.rm = TRUE, names = FALSE)
        out[i, j] <- if (hi) as.integer(values[j] > thr) else -as.integer(values[j] < thr)
      }
    }
  }
  out
}

# runs of flagged positions separated by <= gap bp
runs_of <- function(positions, gap = EVENT_GAP) {
  if (length(positions) == 0) return(integer(0))
  p <- sort(positions)
  grp <- cumsum(c(1L, as.integer(diff(p) > gap)))
  as.integer(table(grp))
}

cat("====================================================\n")
cat("TASK 21 - PER-SAMPLE BURDEN ACROSS THREE SCALES\n")
cat("====================================================\n\n")

cat("Loading source datasets...\n")
normal <- read_src(input_normal)
tumor  <- read_src(input_tumor)

normal22 <- normal[normal$Chromosome == "chr22", ]
normal22 <- normal22[order(normal22$Start), ]
n22 <- nrow(normal22)

spans  <- normal22$Start[N_SITES:n22] - normal22$Start[1:(n22 - N_SITES + 1)]
startB <- which.min(spans)

scales <- list(
  index100 = seq_len(N_SITES),
  tight100 = startB:(startB + N_SITES - 1),
  chr22    = seq_len(n22)
)
cat("scales: index100 (100 CpGs), tight100 (100 CpGs, span ",
    spans[startB], " bp), chr22 (", n22, " CpGs)\n\n", sep = "")

# ------------------------------------------------
# per-sample burden at one scale
# ------------------------------------------------

burden_at_scale <- function(tissue, scale_name, idx, flag_self, flag_ext) {

  src         <- if (tissue == "Normal") normal else tumor
  sample_cols <- paste0(substr(tissue, 1, 1), 1:53)
  sel_ids     <- normal22$Composite.Element.REF[idx]

  m   <- match(sel_ids, src$Composite.Element.REF)
  meta <- src[m, ]

  meth <- as.matrix(meta[, sample_cols, drop = FALSE])
  rownames(meth) <- sel_ids
  storage.mode(meth) <- "numeric"
  states <- meta$methy.state
  pos    <- meta$Start

  cat("  ", tissue, "/", scale_name, ": bio flags over ", nrow(meth),
      " CpGs...\n", sep = "")
  bio     <- build_bio_flag(meth, states, FALSE)
  bio_loo <- build_bio_flag(meth, states, TRUE)

  pick <- function(cache) {
    k <- match(sel_ids, cache$cgID)
    mm <- as.matrix(cache[k, sample_cols, drop = FALSE])
    rownames(mm) <- sel_ids; storage.mode(mm) <- "numeric"; mm
  }
  self <- pick(flag_self)
  ext  <- pick(flag_ext)

  flags <- list(bio = bio, bio.loo = bio_loo, self = self, ext = ext)

  rows <- list()
  for (method in names(flags)) {
    fm    <- flags[[method]]
    total <- sum(fm != 0, na.rm = TRUE)
    for (j in seq_len(ncol(fm))) {
      col <- fm[, j]
      neg <- sum(col == -1, na.rm = TRUE)
      pos1 <- sum(col == 1, na.rm = TRUE)
      fp  <- pos[which(col != 0 & !is.na(col))]
      rl  <- runs_of(fp)
      pval <- if (total > 0)
        stats::binom.test(neg + pos1, total, p = 1/ncol(fm),
                          alternative = "greater")$p.value else NA_real_
      rows[[length(rows) + 1]] <- data.frame(
        tissue = tissue, scale = scale_name, method = method,
        sample = colnames(fm)[j],
        n.cpgs = nrow(fm),
        count.neg1 = neg, count.pos1 = pos1, raw.flags = neg + pos1,
        pct.of.all = if (total > 0) round(100 * (neg + pos1) / total, 3) else NA,
        expected = if (total > 0) round(total / ncol(fm), 3) else NA,
        events.1kb = length(rl),
        epimutations = sum(rl >= MIN_RUN),
        max.run.cpgs = if (length(rl)) max(rl) else 0L,
        flags.per.event = if (length(rl)) round((neg + pos1) / length(rl), 3) else NA,
        binom.p.bonf = signif(min(1, pval * ncol(fm)), 4))
    }
  }
  out <- do.call(rbind, rows)

  # rank within method (1 = most flags)
  out$rank <- ave(-out$raw.flags, out$method, FUN = function(x) rank(x, ties.method = "min"))
  out
}

# ------------------------------------------------
# which samples drive the state-CONTRARY ext flags
# ------------------------------------------------

CONTRARY <- list(H = 1, HM = 1, L = -1, LM = -1)

contrary_drivers <- function(tissue, idx, flag_ext) {

  src         <- if (tissue == "Normal") normal else tumor
  sample_cols <- paste0(substr(tissue, 1, 1), 1:53)
  sel_ids     <- normal22$Composite.Element.REF[idx]
  meta        <- src[match(sel_ids, src$Composite.Element.REF), ]
  states      <- meta$methy.state

  k  <- match(sel_ids, flag_ext$cgID)
  fm <- as.matrix(flag_ext[k, sample_cols, drop = FALSE])
  rownames(fm) <- sel_ids; storage.mode(fm) <- "numeric"

  rows <- list()
  for (j in seq_len(ncol(fm))) {
    tot_contra <- 0L
    per_state  <- setNames(integer(length(CONTRARY)), names(CONTRARY))
    for (st in names(CONTRARY)) {
      keep <- which(states == st)
      if (!length(keep)) next
      hits <- sum(fm[keep, j] == CONTRARY[[st]], na.rm = TRUE)
      per_state[st] <- hits
      tot_contra <- tot_contra + hits
    }
    concordant <- sum(fm[, j] != 0, na.rm = TRUE) - tot_contra
    rows[[length(rows) + 1]] <- data.frame(
      tissue = tissue, sample = colnames(fm)[j],
      contrary.L = per_state[["L"]], contrary.LM = per_state[["LM"]],
      contrary.HM = per_state[["HM"]], contrary.H = per_state[["H"]],
      contrary.total = tot_contra,
      # NOT "concordant": this is every other non-zero flag, which
      # includes all flags at M and R sites, where no state direction
      # is defined at all. On chr22 that is 308 sites in Normal and
      # 2,095 in Tumor, so the label matters.
      other.flags = concordant,
      all.ext.flags = concordant + tot_contra,
      pct.contrary = if (concordant + tot_contra > 0)
        round(100 * tot_contra / (concordant + tot_contra), 2) else NA)
  }
  out <- do.call(rbind, rows)
  out <- out[order(-out$contrary.total), ]
  rownames(out) <- NULL
  out
}

# ------------------------------------------------
# Run
# ------------------------------------------------

for (tissue in c("Normal", "Tumor")) {

  cat("--------------------------------------------------\n")
  cat(tissue, "\n")
  cat("--------------------------------------------------\n")

  cat("  reading flag matrices...\n")
  fself <- read_big(flag_files[[tissue]]$self)
  fext  <- read_big(flag_files[[tissue]]$ext)

  bur <- do.call(rbind, lapply(names(scales), function(s)
    burden_at_scale(tissue, s, scales[[s]], fself, fext)))

  write_out(bur, paste0("Task21_Chr22_SampleBurden_", tissue, ".csv"))

  # ---- scale comparison: rank of every sample at each scale
  cat("\n===== Rank stability across scales -", tissue, "=====\n")
  cat("  A sample that leads at one scale and not the next is a\n")
  cat("  window artifact, not an outlier sample.\n\n")

  cmp_rows <- list()
  for (method in c("bio", "bio.loo", "self", "ext")) {
    b <- bur[bur$method == method, ]
    w <- reshape(b[, c("sample", "scale", "raw.flags", "rank")],
                 idvar = "sample", timevar = "scale", direction = "wide")
    w$method <- method
    w$tissue <- tissue
    cmp_rows[[length(cmp_rows) + 1]] <- w
  }
  cmp <- do.call(rbind, cmp_rows)
  rownames(cmp) <- NULL
  cmp <- cmp[order(cmp$method, cmp$rank.chr22), ]
  write_out(cmp, paste0("Task21_Chr22_ScaleComparison_", tissue, ".csv"))

  for (method in c("bio", "ext")) {
    cat("\n  -- method", method, ": top 8 by chr22-wide burden --\n")
    s <- cmp[cmp$method == method, ]
    s <- s[order(s$rank.chr22), ]
    print(utils::head(s[, c("sample", "raw.flags.index100", "rank.index100",
                            "raw.flags.tight100", "rank.tight100",
                            "raw.flags.chr22", "rank.chr22")], 8),
          row.names = FALSE)

    cat("\n  -- method", method, ": top 8 by the index100 window (the Task 17 view) --\n")
    s2 <- s[order(s$rank.index100), ]
    print(utils::head(s2[, c("sample", "raw.flags.index100", "rank.index100",
                             "raw.flags.tight100", "rank.tight100",
                             "raw.flags.chr22", "rank.chr22")], 8),
          row.names = FALSE)

    # rank correlation between scales
    s3 <- cmp[cmp$method == method, ]
    cat("\n  Spearman rank correlation, index100 vs chr22: ",
        round(cor(s3$rank.index100, s3$rank.chr22, method = "spearman"), 3), "\n", sep = "")
    cat("  Spearman rank correlation, tight100 vs chr22: ",
        round(cor(s3$rank.tight100, s3$rank.chr22, method = "spearman"), 3), "\n", sep = "")
  }

  # ---- genome-wide context from Task 15
  t15f <- task15_files[[tissue]]
  if (file.exists(t15f)) {
    t15 <- read.csv(t15f)
    t15$genomewide.selfref.flags <- t15$count.neg1 + t15$count.1
    t15$genomewide.rank <- rank(-t15$genomewide.selfref.flags, ties.method = "min")
    cat("\n===== Genome-wide context (Task 15, self-reference) -", tissue, "=====\n")
    sub <- bur[bur$method == "ext" & bur$scale == "chr22", ]
    j <- merge(sub[, c("sample", "raw.flags", "rank")],
               t15[, c("sample", "genomewide.selfref.flags", "genomewide.rank")],
               by = "sample")
    names(j)[2:3] <- c("chr22.ext.flags", "chr22.ext.rank")
    j <- j[order(j$chr22.ext.rank), ]
    print(utils::head(j, 10), row.names = FALSE)
    cat("\n  Spearman, chr22 ext burden vs genome-wide self burden: ",
        round(cor(j$chr22.ext.rank, j$genomewide.rank, method = "spearman"), 3), "\n", sep = "")
  }

  # ---- who drives the contrary flags
  cat("\n===== Samples driving state-CONTRARY external flags -", tissue, "=====\n")
  cat("  (chr22-wide; contrary = +1 at an H/HM site or -1 at an L/LM site)\n\n")
  drv <- contrary_drivers(tissue, scales$chr22, fext)
  print(utils::head(drv, 10), row.names = FALSE)
  write_out(drv, paste0("Task21_Chr22_ContraryDrivers_", tissue, ".csv"))

  # ---- focus sample position
  fs <- paste0(substr(tissue, 1, 1), sub("^[NT]", "", FOCUS_SAMPLE))
  cat("\n===== Focus sample", fs, "-", tissue, "=====\n")
  print(bur[bur$sample == fs,
            c("scale", "method", "raw.flags", "expected", "rank",
              "events.1kb", "epimutations", "max.run.cpgs",
              "flags.per.event", "binom.p.bonf")], row.names = FALSE)
  cat("\n")
  rm(fself, fext); invisible(gc())
}

# ------------------------------------------------
# The N37 cluster, all 53 samples, for inspection
# ------------------------------------------------

cat("===== N37's LM cluster in the index100 window =====\n")

sel_ids <- normal22$Composite.Element.REF[scales$index100]
meta    <- normal22[scales$index100, ]
smp     <- paste0("N", 1:53)
meth    <- as.matrix(meta[, smp, drop = FALSE]); storage.mode(meth) <- "numeric"
bio     <- build_bio_flag(meth, meta$methy.state, FALSE)

hit <- which(bio[, "N37"] != 0 & !is.na(bio[, "N37"]))
if (length(hit)) {
  clus <- data.frame(
    cgID = sel_ids[hit], pos = meta$Start[hit],
    methy.state = meta$methy.state[hit],
    N37.beta = round(meth[hit, "N37"], 5),
    cohort.min = round(apply(meth[hit, , drop = FALSE], 1, min, na.rm = TRUE), 5),
    cohort.median = round(apply(meth[hit, , drop = FALSE], 1, median, na.rm = TRUE), 5),
    cohort.max = round(apply(meth[hit, , drop = FALSE], 1, max, na.rm = TRUE), 5),
    second.highest = round(apply(meth[hit, , drop = FALSE], 1,
                                 function(r) sort(r, decreasing = TRUE)[2]), 5))
  clus$margin.over.second <- round(clus$N37.beta - clus$second.highest, 5)
  clus$gap.to.prev.bp <- c(NA, diff(clus$pos))
  clus <- clus[order(clus$pos), ]
  print(clus, row.names = FALSE)
  cat("\n  N37's flags span ", max(clus$pos) - min(clus$pos), " bp in ",
      length(runs_of(clus$pos)), " run(s) at a 1 kb gap.\n", sep = "")
  cat("  Median margin over the second-highest sample: ",
      round(median(clus$margin.over.second), 5), "\n", sep = "")
  cat("\n  Read the margin column before concluding anything. At n = 53 the\n")
  cat("  bio rule flags the cohort maximum whether it wins by 0.001 or 0.4,\n")
  cat("  so the flag alone says nothing. The margin says whether the win\n")
  cat("  was real. Where the margin is large AND the flags form a run of\n")
  cat("  neighbouring CpGs, that is a focal epimutation and the per-sample\n")
  cat("  summary is right for the wrong reason: it counts one event N times.\n")
  cat("  Where the margin is small, the flag is noise.\n")
  write_out(clus, "Task21_Chr22_N37_Cluster_Normal.csv")
}

# ------------------------------------------------
# Part E - magnitude of the flagged differences
# ------------------------------------------------
# The external reference is state-blind, so it can flag a sample at an
# L-state site for a beta difference of 0.01 - statistically past the
# 1st percentile of 2,015 normals, biologically nothing. Beta values
# are heteroscedastic: their SD is compressed below 0.2 and above 0.8
# (Du et al. 2010, BMC Bioinformatics 11:587), which is exactly where
# the L/LM and H/HM sites live. epimutacions guards against this with
# an absolute-difference floor (offset_abs = 0.15).
#
# This block measures how much of the external-reference flag burden
# would survive such a floor, per state.

cat("\n===== Magnitude of external-reference flags by state =====\n")
cat("  |beta - cohort median| for every cell the ext method flagged.\n")
cat("  A flag whose magnitude is under 0.10 is statistically real and\n")
cat("  biologically negligible.\n\n")

mag_rows <- list()
for (tissue in c("Normal", "Tumor")) {

  src         <- if (tissue == "Normal") normal else tumor
  sample_cols <- paste0(substr(tissue, 1, 1), 1:53)
  ids         <- normal22$Composite.Element.REF
  meta        <- src[match(ids, src$Composite.Element.REF), ]

  mth <- as.matrix(meta[, sample_cols, drop = FALSE])
  storage.mode(mth) <- "numeric"
  med <- apply(mth, 1, median, na.rm = TRUE)

  fx <- read_big(flag_files[[tissue]]$ext)
  fm <- as.matrix(fx[match(ids, fx$cgID), sample_cols, drop = FALSE])
  storage.mode(fm) <- "numeric"
  rm(fx); invisible(gc())

  st <- meta$methy.state
  dev <- abs(mth - med)

  for (state in STATE_ORDER) {
    keep <- which(st == state)
    if (!length(keep)) next
    f <- as.vector(fm[keep, , drop = FALSE])
    d <- as.vector(dev[keep, , drop = FALSE])
    ok <- !is.na(f) & f != 0 & !is.na(d)
    if (!any(ok)) next
    d <- d[ok]
    mag_rows[[length(mag_rows) + 1]] <- data.frame(
      tissue = tissue, methy.state = state, n.flags = length(d),
      median.abs.dev = round(median(d), 4),
      q90.abs.dev = round(quantile(d, .9, names = FALSE), 4),
      pct.under.0.05 = round(100 * mean(d < 0.05), 2),
      pct.under.0.10 = round(100 * mean(d < 0.10), 2),
      pct.under.0.15 = round(100 * mean(d < 0.15), 2))
  }
  rm(fm, mth, dev); invisible(gc())
}
mag <- do.call(rbind, mag_rows)
mag$methy.state <- factor(mag$methy.state, levels = STATE_ORDER)
mag <- mag[order(mag$tissue, mag$methy.state), ]
print(mag, row.names = FALSE)
write_out(mag, "Task21_Chr22_ExtFlagMagnitude.csv")

cat("\n====================================================\n")
cat("TASK 21 COMPLETE\n")
cat("====================================================\n")
