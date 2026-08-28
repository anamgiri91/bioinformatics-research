# ================================================================
# Task 19 - Debug the chr22 per-sample flag summary (N37 on top)
# Date: Aug 27, 2026
#
# Addresses item 2 of the pre-meeting email:
#
#   "Debug the sample summary for the 100 CG sites with N37 on top."
#
# Task 17 reported N37 carrying 11 biological flags out of 71 across
# 53 samples, roughly 8x the uniform expectation of 1.34. This script
# establishes whether that is a real per-sample outlier burden, a
# counting artifact, or a single regional event counted many times.
#
# ---------------------------------------------------------------
# WHAT THIS SCRIPT ESTABLISHES
#
# Two separate effects are conflated in the Task 17 summary:
#
# (1) STRUCTURAL. The bio rule thresholds at the 99th/1st percentile
#     of the same 53 values, so at n = 53 exactly one sample clears
#     the threshold at every eligible CpG. The 71 eligible CpGs
#     (L 8 + LM 14 + H 11 + HM 38) therefore produce exactly 71 flags,
#     always. The summary is not counting outlier burden; it is
#     counting how often a sample happens to be the cohort extreme.
#
# (2) SPATIAL. Flags at neighbouring CpGs are not independent events.
#     All 10 of N37's positive flags fall on LM sites inside a single
#     ~9.2 kb window (16,601,097 - 16,610,333), which contains 11 of
#     the 14 LM sites in the selection. N37 is the most-methylated
#     sample at 10 of those 11. That is one regional hypermethylation
#     event, not ten independent outliers - precisely the "OutlierMeth
#     does not use local spatial correlation" concern from the Aug 6
#     notes.
#
# The script therefore reports, per sample, BOTH the raw flag count
# and a collapsed EVENT count that merges flags at nearby sites, at
# three clustering distances so the threshold is not arbitrary.
# ---------------------------------------------------------------
#
# Methods compared (same four as Task 18):
#   bio     - state-aware percentile rule, literal (as specified)
#   bio.loo - same rule, threshold from the other 52 samples
#   self    - OutlierMeth self-referential flags
#   ext     - OutlierMeth external-reference flags
#
# Run for Normal (N1-N53) and Tumor (T1-T53).
#
# NOTE: helper functions are duplicated from Task 18 rather than
# sourced, so each Slurm job remains standalone.
# ================================================================

options(stringsAsFactors = FALSE)

# ------------------------------------------------
# 1. Paths and parameters
# ------------------------------------------------

# Directories resolve to leap2 when present, otherwise to a local
# working copy, so the same script runs unmodified in both places.
pick_dir <- function(...) {
  cand <- c(...)
  for (d in cand) if (dir.exists(d)) return(d)
  stop("no candidate directory exists:\n  ", paste(cand, collapse = "\n  "))
}

src_dir <- pick_dir(
  "/home/s_s355/research3.data/TCGA/filter.BRCA.UCEC.Aug17.2022/sorted.TCGA.380355cg.files",
  path.expand("~/Downloads/tcga_brca_380355cg")
)

input_normal <- file.path(src_dir, "Sorted.BRCA.53Alive.Normal.380355cg.75col.May28.2026.txt")
input_tumor  <- file.path(src_dir, "Sorted.BRCA.53Alive.Tumor.380355cg.75col.May28.2026.txt")

out_dir <- pick_dir(
  "/mmfs1/home/wln26/Experiments.Outlier.July31.2026/Results",
  path.expand("~/Desktop/bioinformatics-research/Results")
)

cat("src_dir:", src_dir, "\nout_dir:", out_dir, "\n\n")

flag_files <- list(
  Normal = list(
    self = file.path(out_dir, "Task13_Normal_selfref_flags.csv"),
    ext  = file.path(out_dir, "Task13_Normal_extref_flags.csv")
  ),
  Tumor = list(
    self = file.path(out_dir, "Task13_Tumor_selfref_flags.csv"),
    ext  = file.path(out_dir, "Task13_Tumor_extref_flags.csv")
  )
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

N_SITES      <- 100
STATE_ORDER  <- c("L", "LM", "M", "HM", "H", "R")
P_LEVEL      <- 0.01
CLUSTER_GAPS <- c(1000, 5000, 10000)   # bp; flags within this distance merge
N_DETAIL     <- 3                       # how many top samples to detail

# ------------------------------------------------
# 2. Helpers
# ------------------------------------------------

read_src <- function(path) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    as.data.frame(data.table::fread(path, sep = "\t", header = TRUE,
                                    check.names = FALSE, showProgress = FALSE))
  } else {
    read.table(path, header = TRUE, sep = "\t", check.names = FALSE)
  }
}

read_big <- function(path) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    as.data.frame(data.table::fread(path, showProgress = FALSE))
  } else {
    read.csv(path, check.names = FALSE)
  }
}

write_out <- function(tbl, fname) {
  path <- file.path(out_dir, fname)
  write.csv(tbl, path, row.names = FALSE)
  cat("  written:", path, "\n")
  invisible(path)
}

build_bio_flag <- function(meth, states, leave_one_out = FALSE) {

  n_cpg <- nrow(meth); n_smp <- ncol(meth)
  out <- matrix(0, nrow = n_cpg, ncol = n_smp,
                dimnames = list(rownames(meth), colnames(meth)))

  for (i in seq_len(n_cpg)) {
    state  <- states[i]
    values <- as.numeric(meth[i, ])
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

subset_flag_file <- function(path, ids, sample_cols) {
  fl  <- read_big(path)
  idx <- match(ids, fl$cgID)
  if (anyNA(idx)) stop("cgIDs missing from ", basename(path), ": ", sum(is.na(idx)))
  m <- as.matrix(fl[idx, sample_cols, drop = FALSE])
  rownames(m) <- ids
  storage.mode(m) <- "numeric"
  m
}

# collapse a sample's flagged positions into events: consecutive
# flagged sites separated by <= max_gap bp count as one event
count_events <- function(positions, max_gap) {
  if (length(positions) == 0) return(0L)
  p <- sort(positions)
  1L + sum(diff(p) > max_gap)
}

cat("====================================================\n")
cat("TASK 19 - CHR22 PER-SAMPLE FLAG SUMMARY DEBUG\n")
cat("====================================================\n\n")

# ------------------------------------------------
# 3. Load and select the same 100 CpGs as Task 17 / 18
# ------------------------------------------------

cat("Loading source datasets...\n")
normal <- read_src(input_normal)
tumor  <- read_src(input_tumor)

normal22 <- normal[normal$Chromosome == "chr22", ]
normal22 <- normal22[order(normal22$Start), ]

n_sel   <- min(N_SITES, nrow(normal22))
sel     <- normal22[seq_len(n_sel), ]
sel_ids <- sel$Composite.Element.REF

cat("  ", n_sel, " chr22 CpGs, ", sel$Start[1], " - ", sel$Start[n_sel],
    " (", sel$Start[n_sel] - sel$Start[1], " bp)\n\n", sep = "")

assemble <- function(tissue) {

  src         <- if (tissue == "Normal") normal else tumor
  sample_cols <- paste0(substr(tissue, 1, 1), 1:53)

  idx  <- match(sel_ids, src$Composite.Element.REF)
  meta <- src[idx, ]

  meth <- as.matrix(meta[, sample_cols, drop = FALSE])
  rownames(meth) <- sel_ids
  storage.mode(meth) <- "numeric"

  states <- meta$methy.state
  names(states) <- sel_ids

  cat("  ", tissue, ": building bio flags...\n", sep = "")
  bio     <- build_bio_flag(meth, states, FALSE)
  bio_loo <- build_bio_flag(meth, states, TRUE)

  cat("  ", tissue, ": loading self / ext flag matrices...\n", sep = "")
  self <- subset_flag_file(flag_files[[tissue]]$self, sel_ids, sample_cols)
  ext  <- subset_flag_file(flag_files[[tissue]]$ext,  sel_ids, sample_cols)

  list(tissue = tissue, meta = meta, meth = meth, states = states,
       pos = meta$Start, sample_cols = sample_cols,
       flags = list(bio = bio, bio.loo = bio_loo, self = self, ext = ext))
}

# ------------------------------------------------
# 4. Structural check - is the total forced?
# ------------------------------------------------

structural_check <- function(A) {

  eligible <- sum(A$states %in% c("L", "LM", "H", "HM"))
  rows <- list()

  for (method in names(A$flags)) {
    fm <- A$flags[[method]]
    tot <- sum(fm != 0, na.rm = TRUE)

    per_cpg <- apply(fm, 1, function(r) sum(r != 0, na.rm = TRUE))

    rows[[length(rows) + 1]] <- data.frame(
      tissue            = A$tissue,
      method            = method,
      n.cpgs            = nrow(fm),
      n.eligible.cpgs   = eligible,
      total.flags       = tot,
      per.cpg.min       = min(per_cpg),
      per.cpg.max       = max(per_cpg),
      per.cpg.distinct  = length(unique(per_cpg)),
      forced.constant   = length(unique(per_cpg[per_cpg > 0])) == 1
    )
  }
  do.call(rbind, rows)
}

# ------------------------------------------------
# 5. Per-sample summary, with events and a binomial test
# ------------------------------------------------

sample_summary <- function(A) {

  rows <- list()

  for (method in names(A$flags)) {

    fm    <- A$flags[[method]]
    total <- sum(fm != 0, na.rm = TRUE)
    n_smp <- ncol(fm)

    for (j in seq_len(n_smp)) {

      col  <- fm[, j]
      neg  <- sum(col == -1, na.rm = TRUE)
      pos  <- sum(col ==  1, na.rm = TRUE)
      tot  <- neg + pos

      flagged_pos <- A$pos[which(col != 0 & !is.na(col))]

      ev <- vapply(CLUSTER_GAPS,
                   function(g) count_events(flagged_pos, g),
                   integer(1))

      # one-sided binomial test against a uniform share of the flags
      pval <- if (total > 0) {
        stats::binom.test(tot, total, p = 1 / n_smp,
                          alternative = "greater")$p.value
      } else NA_real_

      rows[[length(rows) + 1]] <- data.frame(
        tissue        = A$tissue,
        method        = method,
        sample        = colnames(fm)[j],
        count.neg1    = neg,
        count.pos1    = pos,
        total.flags   = tot,
        pct.of.all    = if (total > 0) round(100 * tot / total, 3) else NA,
        expected      = if (total > 0) round(total / n_smp, 3) else NA,
        events.1kb    = ev[1],
        events.5kb    = ev[2],
        events.10kb   = ev[3],
        span.bp       = if (length(flagged_pos) > 1)
                          max(flagged_pos) - min(flagged_pos) else 0L,
        binom.p       = signif(pval, 4),
        binom.p.bonf  = signif(min(1, pval * n_smp), 4)
      )
    }
  }

  out <- do.call(rbind, rows)
  out <- out[order(out$method, -out$total.flags, out$sample), ]
  rownames(out) <- NULL
  out
}

# ------------------------------------------------
# 6. Site-level detail for the top samples
# ------------------------------------------------

top_sample_detail <- function(A, summ, n_top = N_DETAIL) {

  rows <- list()

  for (method in names(A$flags)) {

    fm <- A$flags[[method]]
    sm <- summ[summ$method == method, ]
    sm <- sm[order(-sm$total.flags), ]
    tops <- utils::head(sm$sample[sm$total.flags > 0], n_top)

    for (s in tops) {

      col <- fm[, s]
      hit <- which(col != 0 & !is.na(col))
      if (length(hit) == 0) next

      ord <- order(A$pos[hit])
      hit <- hit[ord]

      prev_pos <- c(NA, A$pos[hit][-length(hit)])

      for (k in seq_along(hit)) {

        i    <- hit[k]
        vals <- as.numeric(A$meth[i, ])
        v    <- vals[which(colnames(fm) == s)]

        rows[[length(rows) + 1]] <- data.frame(
          tissue         = A$tissue,
          method         = method,
          sample         = s,
          cgID           = rownames(fm)[i],
          chr            = A$meta$Chromosome[i],
          pos            = A$pos[i],
          gap.to.prev.bp = A$pos[i] - prev_pos[k],
          methy.state    = A$states[i],
          flag           = col[i],
          beta           = round(v, 5),
          cohort.min     = round(min(vals, na.rm = TRUE), 5),
          cohort.median  = round(median(vals, na.rm = TRUE), 5),
          cohort.max     = round(max(vals, na.rm = TRUE), 5),
          rank.in.cohort = rank(vals, na.last = "keep")[which(colnames(fm) == s)],
          is.cohort.max  = isTRUE(v == max(vals, na.rm = TRUE)),
          is.cohort.min  = isTRUE(v == min(vals, na.rm = TRUE))
        )
      }
    }
  }

  if (length(rows) == 0) return(data.frame(tissue = character(0)))
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# ------------------------------------------------
# 7. Run
# ------------------------------------------------

for (tissue in c("Normal", "Tumor")) {

  cat("--------------------------------------------------\n")
  cat(tissue, "\n")
  cat("--------------------------------------------------\n")

  A <- assemble(tissue)

  cat("\n===== Structural check -", tissue, "=====\n")
  cat("  Is the per-CpG flag count forced to a constant?\n")
  sc <- structural_check(A)
  print(sc, row.names = FALSE)
  write_out(sc, paste0("Task19_Chr22_StructuralCheck_", tissue, ".csv"))

  cat("\n===== Per-sample summary -", tissue, "=====\n")
  summ <- sample_summary(A)

  for (method in names(A$flags)) {
    cat("\n  -- method:", method, "(top 8) --\n")
    sub <- summ[summ$method == method, ]
    print(utils::head(sub[, c("sample", "count.neg1", "count.pos1",
                              "total.flags", "expected",
                              "events.1kb", "events.5kb", "events.10kb",
                              "span.bp", "binom.p.bonf")], 8),
          row.names = FALSE)
  }
  write_out(summ, paste0("Task19_Chr22_SampleSummary_", tissue, ".csv"))

  cat("\n===== Site-level detail, top", N_DETAIL, "samples -", tissue, "=====\n")
  detail <- top_sample_detail(A, summ)
  if (nrow(detail) == 0) {
    cat("  none\n")
  } else {
    cat("  ", nrow(detail), " rows. bio-method rows:\n", sep = "")
    d <- detail[detail$method == "bio", ]
    print(d[, c("sample", "pos", "gap.to.prev.bp", "methy.state",
                "flag", "beta", "cohort.median", "cohort.max",
                "rank.in.cohort")], row.names = FALSE)
  }
  write_out(detail, paste0("Task19_Chr22_TopSampleDetail_", tissue, ".csv"))

  # ---- interpretation aid: raw flags vs collapsed events, bio method
  cat("\n===== Flags vs events (bio method) -", tissue, "=====\n")
  b <- summ[summ$method == "bio" & summ$total.flags > 0, ]
  b <- b[order(-b$total.flags), ]
  cat("  A sample whose flags collapse to few events is showing ONE\n")
  cat("  regional effect, not many independent outliers.\n\n")
  print(utils::head(b[, c("sample", "total.flags", "events.1kb",
                          "events.5kb", "events.10kb", "span.bp")], 8),
        row.names = FALSE)
  cat("\n")
}

cat("\n====================================================\n")
cat("TASK 19 COMPLETE\n")
cat("====================================================\n")
cat("Sites:", n_sel, "consecutive chr22 CpGs\n")
cat("Methods: bio, bio.loo, self, ext\n")
cat("Clustering distances:", paste(CLUSTER_GAPS, collapse = ", "), "bp\n")
cat("\nInterpretation: compare total.flags against events.5kb. Where the\n")
cat("two diverge, the per-sample summary is counting one regional event\n")
cat("many times over, and the sample is not an outlier in the sense the\n")
cat("summary implies.\n")
cat("====================================================\n")
