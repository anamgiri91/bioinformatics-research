# ================================================================
# Task 20 - Chr22 state summary on a TIGHT window, with method
#           concordance and biological annotation of contrary flags
# Date: Aug 28, 2026
#
# Continues Task 18. Task 18 answered the email's item 1 but left
# three gaps, all of which this script closes:
#
#   GAP 1 - the "100 consecutive CG sites" used by Tasks 17/18 are
#           consecutive by INDEX, not close together. They span
#           5,170,347 bp with a single 3.41 Mb gap, so the "explore
#           the local neighbourhood" premise of the email does not
#           hold for them. This script adds a second window: the
#           100 consecutive chr22 CpGs with the SMALLEST genomic
#           span (73,399 bp, chr22:50,459,312-50,532,711). Both
#           windows are analysed side by side so the legacy result
#           stays reproducible and the difference is visible.
#
#   GAP 2 - Task 18 reports counts per method but never measures
#           AGREEMENT between methods. Raw percent-agreement is
#           useless here because all four matrices are >95% zeros
#           (REVIEW.md Finding 4: Task 17's agreement metric counts
#           shared zeros and is near-total by construction). This
#           script reports Cohen's kappa and a non-zero conditional
#           agreement instead.
#
#   GAP 3 - the email asks "if there are non-zero counts, zoom in to
#           see which sample and which site, WHY?". Task 18 produced
#           the which; it did not produce the why. This script
#           annotates every contrary flag with CpG-island context,
#           nearby SNVs, cohort rank, distance past the threshold,
#           and whether neighbouring CpGs in the same sample are
#           flagged too (a lone flag and a flag inside a run of
#           flagged neighbours are different kinds of evidence).
#
# ---------------------------------------------------------------
# METHODS COMPARED (unchanged from Task 18)
#   bio     - state-aware percentile rule, literal
#               L/LM -> +1 if beta > 99th pct of the 53 values
#               H/HM -> -1 if beta <  1st pct of the 53 values
#               M/R  -> no rule, always 0
#   bio.loo - same rule, threshold for sample i from the other 52
#   self    - OutlierMeth flagMeth() vs thresholds from these 53
#   ext     - OutlierMeth flagMeth() vs the external TCGA-GEO panel
#             (747 independent TCGA normal samples, 21 tissue types,
#             Downs, Thursby & Cope 2023, Epigenetics 18:2213874)
#
# DEGENERACY, restated because it governs how to read every table:
# quantile(x, p) over n = 53 is an order statistic of the same 53
# values, so bio flags exactly one sample per eligible CpG and self
# flags exactly one -1 and one +1 per CpG, always, for any p from
# 0.01 to 0.00001. Those columns cannot show a state effect. bio.loo
# removes the guaranteed self-flag; ext is reference-free of the
# problem entirely. Interpret ext; read bio/self as a control that
# shows what "no signal" looks like.
# ---------------------------------------------------------------
#
# Outputs (per tissue x window):
#   Task20_Chr22_<window>_StateSummary_<tissue>.csv
#   Task20_Chr22_<window>_ContraryCounts_<tissue>.csv
#   Task20_Chr22_<window>_ContraryDetail_<tissue>.csv
#   Task20_Chr22_<window>_Concordance_<tissue>.csv
# plus one shared:
#   Task20_Chr22_WindowComparison.csv
# ================================================================

options(stringsAsFactors = FALSE)

# ------------------------------------------------
# 1. Paths
# ------------------------------------------------

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

out_dir  <- file.path(repo_dir, "Results")
data_dir <- file.path(repo_dir, "Data")

cgi_file <- file.path(data_dir, "UCSC.CpGI.chr22.txt")
snv_file <- file.path(data_dir, "SNV.all.txt")

flag_files <- list(
  Normal = list(self = file.path(out_dir, "Task13_Normal_selfref_flags.csv"),
                ext  = file.path(out_dir, "Task13_Normal_extref_flags.csv")),
  Tumor  = list(self = file.path(out_dir, "Task13_Tumor_selfref_flags.csv"),
                ext  = file.path(out_dir, "Task13_Tumor_extref_flags.csv"))
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

N_SITES     <- 100
STATE_ORDER <- c("L", "LM", "M", "HM", "H", "R")   # "R", not "Rc" (REVIEW Finding 2)
P_LEVEL     <- 0.01
NEIGHBOUR_BP <- 1000   # epimutacions' epimutation window (Barbosa et al.)

cat("src_dir :", src_dir,  "\n")
cat("out_dir :", out_dir,  "\n")
cat("data_dir:", data_dir, "\n\n")

# ------------------------------------------------
# 2. Helpers
# ------------------------------------------------

has_dt <- requireNamespace("data.table", quietly = TRUE)

read_src <- function(path) {
  if (has_dt) as.data.frame(data.table::fread(path, sep = "\t", header = TRUE,
                                              check.names = FALSE, showProgress = FALSE))
  else read.table(path, header = TRUE, sep = "\t", check.names = FALSE)
}

read_big <- function(path) {
  if (has_dt) as.data.frame(data.table::fread(path, showProgress = FALSE))
  else read.csv(path, check.names = FALSE)
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

subset_flag_file <- function(cache, ids, sample_cols) {
  idx <- match(ids, cache$cgID)
  if (anyNA(idx)) stop("cgIDs missing from flag matrix: ", sum(is.na(idx)))
  m <- as.matrix(cache[idx, sample_cols, drop = FALSE])
  rownames(m) <- ids
  storage.mode(m) <- "numeric"
  m
}

# Cohen's kappa on paired flag vectors taking values in {-1, 0, 1}.
# Cells where either method is NA are dropped pairwise.
cohens_kappa <- function(a, b) {
  keep <- !is.na(a) & !is.na(b)
  a <- a[keep]; b <- b[keep]
  if (length(a) == 0) return(c(kappa = NA_real_, po = NA_real_, pe = NA_real_))
  lev <- c(-1, 0, 1)
  tab <- table(factor(a, levels = lev), factor(b, levels = lev))
  n   <- sum(tab)
  po  <- sum(diag(tab)) / n
  pe  <- sum(rowSums(tab) * colSums(tab)) / n^2
  k   <- if (pe == 1) NA_real_ else (po - pe) / (1 - pe)
  c(kappa = k, po = po, pe = pe)
}

cat("====================================================\n")
cat("TASK 20 - TIGHT-WINDOW STATE SUMMARY + CONCORDANCE\n")
cat("====================================================\n\n")

# ------------------------------------------------
# 3. Load sources and define the two windows
# ------------------------------------------------

cat("Loading source datasets...\n")
normal <- read_src(input_normal)
tumor  <- read_src(input_tumor)
cat("  Normal:", nrow(normal), "x", ncol(normal),
    " Tumor:", nrow(tumor), "x", ncol(tumor), "\n\n")

normal22 <- normal[normal$Chromosome == "chr22", ]
normal22 <- normal22[order(normal22$Start), ]
cat("chr22 CpGs:", nrow(normal22), " span ",
    min(normal22$Start), "-", max(normal22$Start), "\n\n")

# window A: first 100 by position (what Tasks 17/18 used)
idxA <- seq_len(N_SITES)

# window B: the 100 consecutive CpGs with the smallest genomic span
spans <- normal22$Start[N_SITES:nrow(normal22)] -
         normal22$Start[1:(nrow(normal22) - N_SITES + 1)]
startB <- which.min(spans)
idxB   <- startB:(startB + N_SITES - 1)

windows <- list(
  index100 = list(label = "index100", idx = idxA,
                  note = "first 100 chr22 CpGs by position (Task 17/18 legacy)"),
  tight100 = list(label = "tight100", idx = idxB,
                  note = "100 consecutive chr22 CpGs of minimum genomic span")
)

wrows <- list()
for (w in windows) {
  s <- normal22[w$idx, ]
  g <- diff(s$Start)
  cat("window ", w$label, ": ", w$note, "\n", sep = "")
  cat("  chr22:", s$Start[1], "-", s$Start[N_SITES],
      "  span ", s$Start[N_SITES] - s$Start[1], " bp\n", sep = "")
  cat("  gaps: min ", min(g), "  median ", median(g), "  max ", max(g), " bp\n", sep = "")
  cat("  Normal states: ")
  print(table(factor(s$methy.state, levels = STATE_ORDER)))
  wrows[[length(wrows) + 1]] <- data.frame(
    window = w$label, note = w$note,
    first.pos = s$Start[1], last.pos = s$Start[N_SITES],
    span.bp = s$Start[N_SITES] - s$Start[1],
    gap.min = min(g), gap.median = median(g), gap.max = max(g),
    n.within.1kb.of.neighbour = sum(g <= NEIGHBOUR_BP)
  )
  cat("\n")
}
write_out(do.call(rbind, wrows), "Task20_Chr22_WindowComparison.csv")
cat("\n")

# ------------------------------------------------
# 4. Annotation tracks (the "why")
# ------------------------------------------------

cgi <- NULL
if (file.exists(cgi_file)) {
  cgi <- read.table(cgi_file, header = FALSE, sep = "\t")[, 1:3]
  names(cgi) <- c("chr", "start", "end")
  cgi <- cgi[cgi$chr == "chr22", ]
  cat("CpG islands (chr22):", nrow(cgi), "\n")
} else cat("NOTE: CpG island track not found; island columns will be NA\n")

snv <- NULL
if (file.exists(snv_file)) {
  s <- read.table(snv_file, header = FALSE, sep = "\t")[, 1:2]
  names(s) <- c("chr", "pos")
  snv <- sort(s$pos[s$chr == "chr22"])
  cat("SNV positions (chr22):", length(snv), "\n")
} else cat("NOTE: SNV track not found; SNV columns will be NA\n")
cat("\n")

in_island <- function(pos) {
  if (is.null(cgi)) return(rep(NA, length(pos)))
  vapply(pos, function(p) any(p >= cgi$start & p <= cgi$end), logical(1))
}

dist_to_island <- function(pos) {
  if (is.null(cgi)) return(rep(NA_real_, length(pos)))
  vapply(pos, function(p) {
    if (any(p >= cgi$start & p <= cgi$end)) return(0)
    min(c(abs(cgi$start - p), abs(cgi$end - p)))
  }, numeric(1))
}

snv_within <- function(pos, bp = 10) {
  if (is.null(snv)) return(rep(NA_integer_, length(pos)))
  vapply(pos, function(p) sum(snv >= p - bp & snv <= p + bp), integer(1))
}

# ------------------------------------------------
# 5. Assemble one tissue x window
# ------------------------------------------------

flag_cache <- list()
get_flags <- function(tissue, which) {
  key <- paste(tissue, which)
  if (is.null(flag_cache[[key]])) {
    cat("  reading ", basename(flag_files[[tissue]][[which]]), "...\n", sep = "")
    flag_cache[[key]] <<- read_big(flag_files[[tissue]][[which]])
  }
  flag_cache[[key]]
}

assemble <- function(tissue, w) {

  src         <- if (tissue == "Normal") normal else tumor
  sample_cols <- paste0(substr(tissue, 1, 1), 1:53)

  sel_ids <- normal22$Composite.Element.REF[w$idx]

  idx <- match(sel_ids, src$Composite.Element.REF)
  if (anyNA(idx)) stop("selected cgIDs missing from ", tissue)
  meta <- src[idx, ]

  meth <- as.matrix(meta[, sample_cols, drop = FALSE])
  rownames(meth) <- sel_ids
  storage.mode(meth) <- "numeric"

  states <- meta$methy.state
  names(states) <- sel_ids

  bio     <- build_bio_flag(meth, states, FALSE)
  bio_loo <- build_bio_flag(meth, states, TRUE)
  self    <- subset_flag_file(get_flags(tissue, "self"), sel_ids, sample_cols)
  ext     <- subset_flag_file(get_flags(tissue, "ext"),  sel_ids, sample_cols)

  list(tissue = tissue, window = w$label, meta = meta, meth = meth,
       states = states, pos = meta$Start, sample_cols = sample_cols,
       flags = list(bio = bio, bio.loo = bio_loo, self = self, ext = ext))
}

# ------------------------------------------------
# 6. Part A - state summary
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
      pct   <- function(x) if (n_obs > 0) round(100 * x / n_obs, 3) else NA
      c_neg <- sum(cells == -1, na.rm = TRUE)
      c_zer <- sum(cells ==  0, na.rm = TRUE)
      c_pos <- sum(cells ==  1, na.rm = TRUE)
      rows[[length(rows) + 1]] <- data.frame(
        tissue = A$tissue, window = A$window, method = method, methy.state = st,
        n.cpgs = length(keep), n.samples = ncol(fm), n.cells = length(cells),
        n.evaluated = n_obs, n.NA = sum(is.na(cells)),
        count.neg1 = c_neg, count.0 = c_zer, count.pos1 = c_pos,
        pct.neg1 = pct(c_neg), pct.0 = pct(c_zer), pct.pos1 = pct(c_pos))
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
# 7. Part B - contrary flags
# ------------------------------------------------
# H / HM -> a +1 is contrary (already high; higher is expected)
# L / LM -> a -1 is contrary (already low;  lower  is expected)

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
      # how many DISTINCT CpGs and samples carry a contrary flag
      sub   <- fm[keep, , drop = FALSE]
      rows[[length(rows) + 1]] <- data.frame(
        tissue = A$tissue, window = A$window, method = method,
        methy.state = st, contrary.flag = want,
        n.cpgs = length(keep), n.evaluated = n_obs,
        contrary.count = hits,
        contrary.pct = if (n_obs > 0) round(100 * hits / n_obs, 3) else NA,
        n.cpgs.affected = sum(apply(sub, 1, function(r) any(r == want, na.rm = TRUE))),
        n.samples.affected = sum(apply(sub, 2, function(c) any(c == want, na.rm = TRUE))))
    }
  }
  out <- do.call(rbind, rows)
  out$methy.state <- factor(out$methy.state, levels = c("L", "LM", "HM", "H"))
  out <- out[order(out$methy.state, out$method), ]
  rownames(out) <- NULL
  out
}

# ------------------------------------------------
# 8. Part C - annotated drill-down ("why")
# ------------------------------------------------

contrary_detail <- function(A) {

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
        p    <- A$pos[i]

        # regional context: same sample, same method, CpGs within
        # NEIGHBOUR_BP that also carry any non-zero flag
        near <- which(abs(A$pos - p) <= NEIGHBOUR_BP & seq_along(A$pos) != i)
        n_near_flagged <- if (length(near)) sum(fm[near, j] != 0, na.rm = TRUE) else 0L

        rows[[length(rows) + 1]] <- data.frame(
          tissue = A$tissue, window = A$window, method = method,
          cgID = cg, chr = A$meta$Chromosome[i], pos = p,
          methy.state = st, sample = smp, flag = want,
          beta = round(v, 5),
          cohort.min = round(min(vals, na.rm = TRUE), 5),
          cohort.q25 = round(quantile(vals, .25, na.rm = TRUE, names = FALSE), 5),
          cohort.median = round(median(vals, na.rm = TRUE), 5),
          cohort.q75 = round(quantile(vals, .75, na.rm = TRUE, names = FALSE), 5),
          cohort.max = round(max(vals, na.rm = TRUE), 5),
          cohort.sd  = round(sd(vals, na.rm = TRUE), 5),
          rank.in.cohort = rank(vals, na.last = "keep")[j],
          # robust distance from the cohort centre, the standard
          # Tukey / MAD-z view of "how far out is this really"
          tukey.outlier = {
            q <- quantile(vals, c(.25, .75), na.rm = TRUE, names = FALSE)
            iqr <- q[2] - q[1]
            isTRUE(v < q[1] - 3 * iqr || v > q[2] + 3 * iqr)
          },
          mad.z = {
            m <- median(vals, na.rm = TRUE)
            md <- mad(vals, na.rm = TRUE)
            if (is.na(md) || md == 0) NA_real_ else round((v - m) / md, 3)
          },
          in.cpg.island = in_island(p),
          dist.to.island.bp = dist_to_island(p),
          snv.within.10bp = snv_within(p),
          n.neighbours.1kb = length(near),
          n.neighbours.1kb.flagged = n_near_flagged,
          isolated.flag = n_near_flagged == 0)
      }
    }
  }

  if (length(rows) == 0) return(data.frame(tissue = character(0)))
  out <- do.call(rbind, rows)
  out <- out[order(out$method, out$methy.state, out$pos, out$sample), ]
  rownames(out) <- NULL
  out
}

# ------------------------------------------------
# 9. Part D - concordance between methods
# ------------------------------------------------
# Raw percent agreement is meaningless on matrices that are >95%
# zeros. Two metrics are reported instead:
#   kappa        - Cohen's kappa, chance-corrected
#   agree.nonzero - of the cells where AT LEAST ONE method flagged,
#                   the fraction where both agreed exactly
# The second is the number to quote when asking "do the biological
# rule and the external reference find the same outliers".

concordance <- function(A) {

  methods <- names(A$flags)
  pairs   <- utils::combn(methods, 2, simplify = FALSE)
  rows    <- list()

  strata <- c(list(ALL = seq_along(A$states)),
              split(seq_along(A$states), A$states))

  for (pr in pairs) {
    for (sname in names(strata)) {
      keep <- strata[[sname]]
      if (length(keep) == 0) next
      a <- as.vector(A$flags[[pr[1]]][keep, , drop = FALSE])
      b <- as.vector(A$flags[[pr[2]]][keep, , drop = FALSE])
      ok <- !is.na(a) & !is.na(b)
      a <- a[ok]; b <- b[ok]
      if (length(a) == 0) next

      kk  <- cohens_kappa(a, b)
      nz  <- (a != 0) | (b != 0)
      both <- (a != 0) & (b != 0)

      rows[[length(rows) + 1]] <- data.frame(
        tissue = A$tissue, window = A$window,
        method.a = pr[1], method.b = pr[2],
        stratum = sname,
        n.compared = length(a),
        n.flag.a = sum(a != 0), n.flag.b = sum(b != 0),
        n.either = sum(nz), n.both = sum(both),
        n.agree.nonzero = sum(a == b & nz),
        agree.nonzero = if (sum(nz) > 0) round(sum(a == b & nz) / sum(nz), 4) else NA,
        jaccard = if (sum(nz) > 0) round(sum(both) / sum(nz), 4) else NA,
        raw.agreement = round(kk["po"], 4),
        expected.agreement = round(kk["pe"], 4),
        kappa = round(kk["kappa"], 4))
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# ------------------------------------------------
# 10. Run
# ------------------------------------------------

for (wname in names(windows)) {
  w <- windows[[wname]]

  for (tissue in c("Normal", "Tumor")) {

    cat("==================================================\n")
    cat("window ", w$label, " | tissue ", tissue, "\n", sep = "")
    cat("==================================================\n")

    A <- assemble(tissue, w)

    cat("\n-- state distribution (", tissue, " annotation) --\n", sep = "")
    print(table(factor(A$states, levels = STATE_ORDER)))

    cat("\n-- Part A: flag distribution by state --\n")
    summ <- state_summary(A)
    print(summ[, c("method", "methy.state", "n.cpgs", "n.evaluated",
                   "count.neg1", "count.0", "count.pos1",
                   "pct.neg1", "pct.0", "pct.pos1")], row.names = FALSE)
    write_out(summ, sprintf("Task20_Chr22_%s_StateSummary_%s.csv", w$label, tissue))

    cat("\n-- Part B: contrary flags (H/HM given +1, L/LM given -1) --\n")
    contra <- contrary_counts(A)
    print(contra[, c("method", "methy.state", "contrary.flag", "n.cpgs",
                     "n.evaluated", "contrary.count", "contrary.pct",
                     "n.cpgs.affected", "n.samples.affected")], row.names = FALSE)
    write_out(contra, sprintf("Task20_Chr22_%s_ContraryCounts_%s.csv", w$label, tissue))

    cat("\n-- Part C: annotated contrary-flag detail --\n")
    detail <- contrary_detail(A)
    if (nrow(detail) == 0) {
      cat("  none found\n")
    } else {
      cat("  ", nrow(detail), " contrary flags\n", sep = "")
      print(table(detail$method, detail$methy.state))
      cat("\n  of the ext-method contrary flags:\n")
      d <- detail[detail$method == "ext", ]
      if (nrow(d)) {
        cat("    in a CpG island        : ", sum(d$in.cpg.island, na.rm = TRUE),
            " / ", nrow(d), "\n", sep = "")
        cat("    with an SNV within 10bp: ", sum(d$snv.within.10bp > 0, na.rm = TRUE),
            " / ", nrow(d), "\n", sep = "")
        cat("    also a Tukey 3xIQR outlier within this cohort: ",
            sum(d$tukey.outlier, na.rm = TRUE), " / ", nrow(d), "\n", sep = "")
        cat("    isolated (no flagged neighbour within 1kb)   : ",
            sum(d$isolated.flag, na.rm = TRUE), " / ", nrow(d), "\n", sep = "")
        cat("\n  ext contrary flags, most extreme by |mad.z| (top 15):\n")
        dd <- d[order(-abs(d$mad.z)), ]
        print(utils::head(dd[, c("cgID", "pos", "methy.state", "sample", "flag",
                                 "beta", "cohort.median", "rank.in.cohort",
                                 "mad.z", "tukey.outlier", "in.cpg.island",
                                 "snv.within.10bp", "n.neighbours.1kb.flagged")], 15),
              row.names = FALSE)
      }
    }
    write_out(detail, sprintf("Task20_Chr22_%s_ContraryDetail_%s.csv", w$label, tissue))

    cat("\n-- Part D: method concordance --\n")
    cat("  agree.nonzero = of cells flagged by at least one method,\n")
    cat("  the fraction where both methods gave the same value.\n")
    conc <- concordance(A)
    print(conc[conc$stratum == "ALL",
               c("method.a", "method.b", "n.either", "n.both",
                 "agree.nonzero", "jaccard", "raw.agreement", "kappa")],
          row.names = FALSE)
    cat("\n  bio vs ext, by state:\n")
    print(conc[conc$method.a == "bio" & conc$method.b == "ext",
               c("stratum", "n.flag.a", "n.flag.b", "n.either", "n.both",
                 "agree.nonzero", "jaccard", "kappa")], row.names = FALSE)
    write_out(conc, sprintf("Task20_Chr22_%s_Concordance_%s.csv", w$label, tissue))
    cat("\n")
  }
}

cat("====================================================\n")
cat("TASK 20 COMPLETE\n")
cat("====================================================\n")
cat("Windows : index100 (Task 17/18 legacy), tight100 (min span)\n")
cat("Tissues : Normal, Tumor\n")
cat("Methods : bio, bio.loo, self, ext\n")
cat("p level :", P_LEVEL, "\n")
cat("\nRead the ext rows. bio and self are pinned by order statistics\n")
cat("at n = 53 and are included as a negative control only.\n")
cat("====================================================\n")
