# ================================================================
# Task 18 - Chromosome 22 flag summary by methylation state
# Date: Aug 27, 2026
#
# Addresses item 1 of the pre-meeting email:
#
#   "For chr22, summarize by state ... example of 100 CG sites
#    H, 11 CG sites, % of -1, 0, 1 for bio/self, % of -1, 0, 1 for
#    external-ref, counts of these two methods.
#    For the H and HM, see if there is any flag.count for 1.
#    For the L and LM, see if there is any flag.count for -1.
#    If there are non-zero counts, zoom in to see which sample and
#    which site, why?"
#
# Uses the SAME 100 consecutive chr22 CpGs selected in Task 17 - the
# email's "H, 11 CG sites" matches that selection exactly (H = 11).
#
# THREE flagging methods are compared per CpG x sample cell:
#   bio  - state-aware percentile rule from Task 17
#            L/LM -> +1 if beta > 99th percentile of the 53 values
#            H/HM -> -1 if beta <  1st percentile of the 53 values
#            M/R  -> no rule, always 0
#   self - OutlierMeth flagMeth() against thresholds derived from the
#          same 53 samples          (Task13_*_selfref_flags.csv)
#   ext  - OutlierMeth flagMeth() against the external TCGA panel
#          of ~2,000 independent samples (Task13_*_extref_flags.csv)
#
# A fourth column, bio.loo, repeats the bio rule leave-one-out: the
# threshold for sample i is computed from the other 52 samples. See
# the note on degeneracy below for why this is included.
#
# ---------------------------------------------------------------
# NOTE ON DEGENERACY (read before interpreting the bio/self columns)
#
# quantile(x, 0.99) over n = 53 values interpolates between the 52nd
# and 53rd order statistics, so exactly one sample can exceed it; the
# 0.01 quantile likewise admits exactly one below. The literal bio
# rule therefore assigns exactly ONE flag to every eligible CpG, and
# the self-referential OutlierMeth flags assign exactly one -1 and one
# +1 to every CpG, regardless of the data. Task 17 shows this: 71
# eligible CpGs (L 8 + LM 14 + H 11 + HM 38) produced exactly 71 flags.
#
# The bio and self columns are reproduced here as specified so the
# comparison the email asks for can be made, but their per-CpG counts
# carry no biological signal. bio.loo removes the guaranteed self-flag.
# The ext column is unaffected and is the one that can be interpreted.
# ---------------------------------------------------------------
#
# NOTE ON methy.state ACROSS TISSUES: methy.state is annotated per
# tissue, so the same cgID can sit in different states in Normal and
# Tumor. Each tissue is summarized under its OWN state annotation;
# the site set (100 cgIDs) is held constant so the two are comparable.
# Both state distributions are printed for reference.
#
# NOTE ON PHYSICAL SPAN: these 100 CpGs are consecutive by index, not
# tightly clustered. They span 11,915,060 - 17,085,407 (5.17 Mb) with
# a median inter-site gap of 450 bp but a single 3.41 Mb gap. The
# email asked for sites "close to each other"; this is flagged in the
# output so a tighter window can be chosen for follow-up.
#
# Run for Normal (N1-N53) and Tumor (T1-T53).
# ================================================================

options(stringsAsFactors = FALSE)

# ------------------------------------------------
# 1. Paths
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

tcga_ref <- Filter(file.exists, c(
  "/mmfs1/home/wln26/OutlierMeth_tarball_check/OutlierMeth/data/tcga.rda",
  path.expand("~/Downloads/tcga.rda")
))
tcga_ref <- if (length(tcga_ref)) tcga_ref[1] else NA_character_

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

N_SITES     <- 100
STATE_ORDER <- c("L", "LM", "M", "HM", "H", "R")   # "R", not "Rc"
P_LEVEL     <- 0.01                                 # matches Task 13 flagging

# ------------------------------------------------
# 2. Helpers
# ------------------------------------------------

# fread when available (the flag CSVs are ~61 MB each), base R otherwise
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

cat("====================================================\n")
cat("TASK 18 - CHR22 FLAG SUMMARY BY METHYLATION STATE\n")
cat("====================================================\n\n")

# ------------------------------------------------
# 3. Select the 100 CpGs (from Normal, by position)
# ------------------------------------------------

cat("Loading Normal source dataset...\n")
normal <- read_src(input_normal)

cat("Loading Tumor source dataset...\n")
tumor <- read_src(input_tumor)

cat("  Normal:", nrow(normal), "rows x", ncol(normal), "cols\n")
cat("  Tumor :", nrow(tumor),  "rows x", ncol(tumor),  "cols\n\n")

normal22 <- normal[normal$Chromosome == "chr22", ]
normal22 <- normal22[order(normal22$Start), ]

n_sel <- min(N_SITES, nrow(normal22))
sel   <- normal22[seq_len(n_sel), ]
sel_ids <- sel$Composite.Element.REF

cat("Selected", n_sel, "consecutive chr22 CpGs (by position, from Normal).\n")
cat("  first position:", sel$Start[1], "\n")
cat("  last  position:", sel$Start[n_sel], "\n")
cat("  physical span :", sel$Start[n_sel] - sel$Start[1], "bp\n")

gaps <- diff(sel$Start)
cat("  inter-site gaps: min", min(gaps), " median", median(gaps),
    " max", max(gaps), "bp\n")
if (max(gaps) > 1e6) {
  cat("  WARNING: largest gap exceeds 1 Mb - these sites are consecutive by\n")
  cat("           index but not tightly clustered. Consider a narrower window.\n")
}
cat("\n")

# ------------------------------------------------
# 4. Per-tissue assembly
# ------------------------------------------------

# state-aware biological rule -----------------------------------
#   literal: threshold from all 53 values (as specified in the email)
#   loo    : threshold for sample i from the other 52
build_bio_flag <- function(meth, states, leave_one_out = FALSE) {

  n_cpg <- nrow(meth)
  n_smp <- ncol(meth)

  out <- matrix(0, nrow = n_cpg, ncol = n_smp,
                dimnames = list(rownames(meth), colnames(meth)))

  for (i in seq_len(n_cpg)) {

    state  <- states[i]
    values <- as.numeric(meth[i, ])

    if (all(is.na(values))) { out[i, ] <- NA; next }
    if (!state %in% c("L", "LM", "H", "HM")) next   # M / R: no rule

    hi <- state %in% c("L", "LM")   # L/LM look for unusually HIGH
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

# pull the 100 selected rows out of a Task 13 flag matrix
subset_flag_file <- function(path, ids, sample_cols) {
  fl  <- read_big(path)
  idx <- match(ids, fl$cgID)
  if (anyNA(idx)) {
    stop("cgIDs missing from ", basename(path), ": ", sum(is.na(idx)))
  }
  m <- as.matrix(fl[idx, sample_cols, drop = FALSE])
  rownames(m) <- ids
  storage.mode(m) <- "numeric"
  m
}

assemble <- function(tissue) {

  cat("--------------------------------------------------\n")
  cat("Assembling", tissue, "\n")
  cat("--------------------------------------------------\n")

  src         <- if (tissue == "Normal") normal else tumor
  sample_cols <- paste0(substr(tissue, 1, 1), 1:53)   # N1..N53 / T1..T53

  idx <- match(sel_ids, src$Composite.Element.REF)
  if (anyNA(idx)) stop("selected cgIDs missing from ", tissue, " source")
  meta <- src[idx, ]

  meth <- as.matrix(meta[, sample_cols, drop = FALSE])
  rownames(meth) <- sel_ids
  storage.mode(meth) <- "numeric"

  states <- meta$methy.state
  names(states) <- sel_ids

  cat("  state distribution over the ", n_sel, " sites:\n", sep = "")
  print(table(factor(states, levels = STATE_ORDER)))

  cat("  building bio flags (literal)...\n")
  bio <- build_bio_flag(meth, states, leave_one_out = FALSE)

  cat("  building bio flags (leave-one-out)...\n")
  bio_loo <- build_bio_flag(meth, states, leave_one_out = TRUE)

  cat("  loading self-reference flags...\n")
  self <- subset_flag_file(flag_files[[tissue]]$self, sel_ids, sample_cols)

  cat("  loading external-reference flags...\n")
  ext <- subset_flag_file(flag_files[[tissue]]$ext, sel_ids, sample_cols)

  n_na_ext <- sum(apply(ext, 1, function(r) all(is.na(r))))
  cat("  CpGs absent from the external panel (all-NA rows):", n_na_ext, "\n\n")

  list(tissue = tissue, meta = meta, meth = meth, states = states,
       sample_cols = sample_cols,
       flags = list(bio = bio, bio.loo = bio_loo, self = self, ext = ext))
}

# ------------------------------------------------
# 5. Part A - summarize by state, per method
# ------------------------------------------------

state_summary <- function(A) {

  rows <- list()

  for (method in names(A$flags)) {

    fm <- A$flags[[method]]

    for (st in STATE_ORDER) {

      keep <- which(A$states == st)
      if (length(keep) == 0) next

      cells <- as.vector(fm[keep, , drop = FALSE])
      n_obs <- sum(!is.na(cells))
      n_na  <- sum(is.na(cells))

      c_neg <- sum(cells == -1, na.rm = TRUE)
      c_zer <- sum(cells ==  0, na.rm = TRUE)
      c_pos <- sum(cells ==  1, na.rm = TRUE)

      pct <- function(x) if (n_obs > 0) round(100 * x / n_obs, 3) else NA

      rows[[length(rows) + 1]] <- data.frame(
        tissue      = A$tissue,
        method      = method,
        methy.state = st,
        n.cpgs      = length(keep),
        n.samples   = ncol(fm),
        n.cells     = length(cells),
        n.evaluated = n_obs,
        n.NA        = n_na,
        count.neg1  = c_neg,
        count.0     = c_zer,
        count.pos1  = c_pos,
        pct.neg1    = pct(c_neg),
        pct.0       = pct(c_zer),
        pct.pos1    = pct(c_pos)
      )
    }
  }

  out <- do.call(rbind, rows)
  out$method      <- factor(out$method, levels = names(A$flags))
  out$methy.state <- factor(out$methy.state, levels = STATE_ORDER)
  out <- out[order(out$methy.state, out$method), ]
  rownames(out) <- NULL
  out
}

# ------------------------------------------------
# 6. Part B - flags that contradict the state
# ------------------------------------------------
#   H / HM  ->  a +1 is contrary  (already high; "higher" is expected)
#   L / LM  ->  a -1 is contrary  (already low;  "lower"  is expected)
#
# The bio rule cannot produce a contrary flag by construction - it only
# ever assigns +1 at L/LM and -1 at H/HM - so bio and bio.loo are
# expected to be structurally zero here. self and ext are the informative
# columns, and ext is the only one free of the degeneracy.

CONTRARY <- list(H = 1, HM = 1, L = -1, LM = -1)

contrary_counts <- function(A) {

  rows <- list()

  for (method in names(A$flags)) {
    fm <- A$flags[[method]]
    for (st in names(CONTRARY)) {

      keep <- which(A$states == st)
      if (length(keep) == 0) next

      want  <- CONTRARY[[st]]
      cells <- as.vector(fm[keep, , drop = FALSE])
      n_obs <- sum(!is.na(cells))
      hits  <- sum(cells == want, na.rm = TRUE)

      rows[[length(rows) + 1]] <- data.frame(
        tissue         = A$tissue,
        method         = method,
        methy.state    = st,
        contrary.flag  = want,
        n.cpgs         = length(keep),
        n.evaluated    = n_obs,
        contrary.count = hits,
        contrary.pct   = if (n_obs > 0) round(100 * hits / n_obs, 3) else NA
      )
    }
  }

  out <- do.call(rbind, rows)
  out$methy.state <- factor(out$methy.state, levels = c("L", "LM", "HM", "H"))
  out <- out[order(out$methy.state, out$method), ]
  rownames(out) <- NULL
  out
}

# ------------------------------------------------
# 7. Part C - drill-down: which sample, which site, why
# ------------------------------------------------

contrary_detail <- function(A, tcga = NULL) {

  rows <- list()

  for (method in names(A$flags)) {
    fm <- A$flags[[method]]

    for (st in names(CONTRARY)) {
      want <- CONTRARY[[st]]
      keep <- which(A$states == st)
      if (length(keep) == 0) next

      hit <- which(fm[keep, , drop = FALSE] == want, arr.ind = TRUE)
      if (nrow(hit) == 0) next

      for (k in seq_len(nrow(hit))) {

        i <- keep[hit[k, "row"]]
        j <- hit[k, "col"]

        cg   <- rownames(fm)[i]
        smp  <- colnames(fm)[j]
        vals <- as.numeric(A$meth[i, ])
        v    <- vals[j]

        # where the external thresholds sat for this CpG, when available
        pthr <- nthr <- NA_real_
        if (!is.null(tcga) && cg %in% rownames(tcga)) {
          pthr <- as.numeric(tcga[cg, paste0("P", P_LEVEL)])
          nthr <- as.numeric(tcga[cg, paste0("N", P_LEVEL)])
        }

        rows[[length(rows) + 1]] <- data.frame(
          tissue        = A$tissue,
          method        = method,
          cgID          = cg,
          chr           = A$meta$Chromosome[i],
          pos           = A$meta$Start[i],
          methy.state   = st,
          sample        = smp,
          flag          = want,
          beta          = round(v, 5),
          cohort.min    = round(min(vals, na.rm = TRUE), 5),
          cohort.median = round(median(vals, na.rm = TRUE), 5),
          cohort.max    = round(max(vals, na.rm = TRUE), 5),
          rank.in.cohort = rank(vals, na.last = "keep")[j],
          ext.P.thresh  = round(pthr, 5),
          ext.N.thresh  = round(nthr, 5),
          margin.beyond = round(
            if (want == 1 && !is.na(pthr)) v - pthr
            else if (want == -1 && !is.na(nthr)) nthr - v
            else NA_real_, 5)
        )
      }
    }
  }

  if (length(rows) == 0) {
    return(data.frame(tissue = character(0), method = character(0),
                      cgID = character(0)))
  }

  out <- do.call(rbind, rows)
  out <- out[order(out$method, out$methy.state, out$pos, out$sample), ]
  rownames(out) <- NULL
  out
}

# ------------------------------------------------
# 8. Run
# ------------------------------------------------

tcga <- NULL
if (!is.na(tcga_ref)) {
  cat("Loading external TCGA reference panel for threshold lookup...\n")
  load(tcga_ref)              # provides `tcga`
  cat("  tcga:", nrow(tcga), "CpGs x", ncol(tcga), "threshold columns\n\n")
} else {
  cat("NOTE: tcga.rda not found; ext.P.thresh / ext.N.thresh will be NA\n\n")
}

for (tissue in c("Normal", "Tumor")) {

  A <- assemble(tissue)

  cat("===== Part A: flag distribution by state -", tissue, "=====\n")
  summ <- state_summary(A)
  print(summ, row.names = FALSE)
  write_out(summ, paste0("Task18_Chr22_StateSummary_", tissue, ".csv"))

  cat("\n===== Part B: flags contrary to state -", tissue, "=====\n")
  cat("  (H/HM flagged +1, or L/LM flagged -1)\n")
  contra <- contrary_counts(A)
  print(contra, row.names = FALSE)
  write_out(contra, paste0("Task18_Chr22_ContraryFlagCounts_", tissue, ".csv"))

  cat("\n===== Part C: contrary-flag detail -", tissue, "=====\n")
  detail <- contrary_detail(A, tcga)
  if (nrow(detail) == 0) {
    cat("  none found\n")
  } else {
    cat("  ", nrow(detail), " contrary flags; by method:\n", sep = "")
    print(table(detail$method, detail$methy.state))
    cat("\n  first 20 rows:\n")
    print(utils::head(detail, 20), row.names = FALSE)
  }
  write_out(detail, paste0("Task18_Chr22_ContraryFlagDetail_", tissue, ".csv"))
  cat("\n")
}

# ------------------------------------------------
# 9. Selected region reference copy
# ------------------------------------------------

region <- data.frame(
  cgID                = sel_ids,
  chr                 = sel$Chromosome,
  pos                 = sel$Start,
  end                 = sel$End,
  methy.state.normal  = sel$methy.state,
  methy.state.tumor   = tumor$methy.state[match(sel_ids, tumor$Composite.Element.REF)],
  gap.to.previous.bp  = c(NA, diff(sel$Start))
)
write_out(region, "Task18_Chr22_Selected100CpGs_BothTissues.csv")

cat("\n====================================================\n")
cat("TASK 18 COMPLETE\n")
cat("====================================================\n")
cat("Sites:", n_sel, "consecutive chr22 CpGs\n")
cat("Methods: bio (literal), bio.loo (leave-one-out), self, ext\n")
cat("Tissues: Normal, Tumor\n")
cat("Significance level: p =", P_LEVEL, "\n")
cat("\nReminder: bio and self per-CpG counts are fixed by n and p and\n")
cat("carry no biological signal. Interpret the ext column.\n")
cat("====================================================\n")
