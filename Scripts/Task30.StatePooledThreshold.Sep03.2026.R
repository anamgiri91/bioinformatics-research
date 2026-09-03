# ================================================================
# Task 30 - The "bio-stat" rule: one threshold per STATE, not per CpG
# Date: Sep 3, 2026
#
# THE PROPOSAL (2026-09-03 email, "Analysis to-do" item 2)
#
#   "For the L (8 CG sites), bio-stat flag based on all 'L' state,
#    that is, >96% (or lower a bit) percentile of 8*53, similarly
#    for other states."
#
# i.e. stop estimating a threshold from 53 numbers at one CpG, and
# estimate it once from every cell in the state - 8 sites x 53
# samples = 424 values in the window, 2,211 x 53 = 117,183 on all of
# chr22. This is a direct answer to the n = 53 degeneracy: an
# empirical quantile over 53 values is an order statistic and admits
# exactly one flag per CpG by construction (README, "the
# self-reference degeneracy"), while a quantile over 424 or 117,183
# values is an actual estimate.
#
# ---------------------------------------------------------------
# TWO WAYS TO POOL, AND THE DIFFERENCE MATTERS
#
#   bio.stat.abs   pool the raw beta values of the state and cut at
#                  the q-th percentile. This is the email's version
#                  literally. The threshold is ONE absolute beta
#                  value for the whole state.
#
#   bio.stat.dev   pool (beta - that site's cohort median) and cut
#                  the deviation at the q-th percentile. The
#                  threshold is one EFFECT SIZE for the whole state,
#                  re-centred at each site.
#
# They behave very differently and the difference is the point.
# States are broad: "L" spans beta 0.005 to 0.10. Under the absolute
# rule the top of the pooled tail is dominated by whichever SITES sit
# highest within the state, so the rule flags whole sites rather than
# unusual samples - every sample at a high-ish L site clears the bar
# and no sample at a low L site ever can. The deviation rule removes
# the site's own level first, so it asks the intended question: is
# this sample unusual FOR THIS SITE, by an amount large enough to
# matter across the state?
#
# bio.stat.dev is also, precisely, a data-calibrated version of the
# "|beta - median| > 0.05" floor proposed in the same email: instead
# of guessing 0.05, it reads the number off the state's own
# deviation distribution. Both are reported so the guess can be
# checked against the estimate.
#
# ---------------------------------------------------------------
# DIRECTION
#
# For each state only the informative tail is used, matching the
# existing `bio` rule:  L, LM -> upper tail (+1 only)
#                      H, HM -> lower tail (-1 only)
#                      M, R  -> both tails, q/2 in each
#
# ---------------------------------------------------------------
# SCOPES
#   window   the 100 chr22 CpGs from the email (state pools of
#            53 x {8,14,24,38,11,5} cells)
#   chr22    all 6,809 chr22 CpGs, where the pooled quantile is
#            estimated from tens of thousands of values
#
# The window numbers answer the question as asked; the chr22 numbers
# say whether the answer was an artefact of pooling only 424 values.
#
# SCORING mirrors Task 24 so the new rule can be dropped straight
# into that table: flag rate, stability under N(0, 0.01) with full
# re-estimation, magnitude, and concentration over sites and samples.
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

N_WINDOW <- 100
STATE_ORDER <- c("L", "LM", "M", "HM", "H", "R")
QGRID  <- c(0.90, 0.92, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99)
Q_MAIN <- 0.96                       # the email's number
NOISE_SD <- 0.01
N_REPS   <- 10

# informative tail per state
TAIL <- c(L = "upper", LM = "upper", M = "both",
          HM = "lower", H = "lower", R = "both")

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
gini <- function(x) {
  if (!length(x) || sum(x) == 0) return(NA_real_)
  x <- sort(as.numeric(x)); n <- length(x)
  sum((2 * seq_len(n) - n - 1) * x) / (n * sum(x))
}
jaccard <- function(a, b) {
  ok <- !is.na(a) & !is.na(b)
  A <- a[ok] != 0; B2 <- b[ok] != 0
  u <- sum(A | B2); if (!u) return(NA_real_)
  sum(A & B2 & a[ok] == b[ok]) / u
}

# ---------------------------------------------------------------
# The rule. X is the matrix the pool is taken over: raw beta for
# "abs", beta minus the row median for "dev". Thresholds are fitted
# on X and applied to X, so passing a perturbed X re-estimates
# everything - which is what the stability test needs.
# ---------------------------------------------------------------
state_pooled_flags <- function(X, st, q, states = STATE_ORDER) {
  o <- matrix(0L, nrow(X), ncol(X), dimnames = dimnames(X))
  thr <- list()
  for (s in states) {
    k <- which(st == s); if (!length(k)) next
    pool <- as.vector(X[k, , drop = FALSE]); pool <- pool[!is.na(pool)]
    if (length(pool) < 20) next
    tl <- TAIL[[s]]
    hi <- lo <- NA_real_
    if (tl == "upper") hi <- quantile(pool, q, names = FALSE)
    if (tl == "lower") lo <- quantile(pool, 1 - q, names = FALSE)
    if (tl == "both") { hi <- quantile(pool, 1 - (1 - q) / 2, names = FALSE)
                        lo <- quantile(pool, (1 - q) / 2, names = FALSE) }
    Xs <- X[k, , drop = FALSE]
    os <- matrix(0L, nrow(Xs), ncol(Xs))
    if (is.finite(hi)) os[!is.na(Xs) & Xs > hi] <-  1L
    if (is.finite(lo)) os[!is.na(Xs) & Xs < lo] <- -1L
    o[k, ] <- os
    thr[[s]] <- c(n.pool = length(pool), hi = hi, lo = lo)
  }
  o[is.na(X)] <- NA
  attr(o, "thr") <- thr
  o
}

cat("=========================================================\n")
cat("TASK 30 - STATE-POOLED (\"bio-stat\") THRESHOLDS\n")
cat("=========================================================\n\n")

thr_rows <- list(); score_rows <- list(); sweep_rows <- list()
site_rows <- list(); sample_rows <- list()

for (tissue in c("Normal", "Tumor")) {

  cat("#########################################################\n")
  cat("### ", tissue, "\n", sep = "")
  cat("#########################################################\n\n")

  cols <- paste0(substr(tissue, 1, 1), 1:53)

  nrm <- read_any(inputs$Normal, tab = TRUE)
  n22 <- nrm[nrm$Chromosome == "chr22", ]; n22 <- n22[order(n22$Start), ]
  win_ids <- n22$Composite.Element.REF[seq_len(N_WINDOW)]
  rm(nrm, n22); invisible(gc())

  d <- read_any(inputs[[tissue]], tab = TRUE)
  d <- d[d$Chromosome == "chr22", ]; d <- d[order(d$Start), ]
  ids_all <- d$Composite.Element.REF
  B_all <- as.matrix(d[, cols, drop = FALSE]); storage.mode(B_all) <- "numeric"
  rownames(B_all) <- ids_all
  st_all <- d$methy.state
  pos_all <- d$Start
  rm(d); invisible(gc())

  fx <- read_any(flag_ext[[tissue]])
  F_all <- as.matrix(fx[match(ids_all, fx$cgID), cols, drop = FALSE])
  storage.mode(F_all) <- "numeric"; rownames(F_all) <- ids_all
  rm(fx); invisible(gc())

  for (scope in c("window", "chr22")) {
    idx <- if (scope == "window") match(win_ids, ids_all) else seq_along(ids_all)
    idx <- idx[!is.na(idx)]
    B <- B_all[idx, , drop = FALSE]; F <- F_all[idx, , drop = FALSE]
    st <- st_all[idx]; pos <- pos_all[idx]
    med <- apply(B, 1, median, na.rm = TRUE)
    D <- B - med                     # signed deviation from the site median
    dev <- abs(D)

    cat("--------------------------------------------------\n")
    cat("scope = ", scope, "   sites = ", nrow(B), "\n", sep = "")
    cat("--------------------------------------------------\n")

    variants <- list(bio.stat.abs = B, bio.stat.dev = D)

    # ---- thresholds and per-state behaviour at the email's q ----
    for (vn in names(variants)) {
      O <- state_pooled_flags(variants[[vn]], st, Q_MAIN)
      th <- attr(O, "thr")
      cat("\n", vn, "  at q = ", Q_MAIN, "\n", sep = "")
      for (s in STATE_ORDER) {
        k <- which(st == s); if (!length(k)) next
        t_s <- th[[s]]
        Os <- O[k, , drop = FALSE]
        nf <- sum(Os != 0, na.rm = TRUE)
        per_site <- rowSums(Os != 0, na.rm = TRUE)
        per_samp <- colSums(Os != 0, na.rm = TRUE)
        mg <- dev[k, , drop = FALSE][!is.na(Os) & Os != 0]
        cat(sprintf("   %-3s pool=%7d  thr(lo,hi)=(%s,%s)  flags=%5d (%5.2f%%)  sites>=1=%4d/%-4d  samples>=1=%2d/53  med|db|=%.4f\n",
            s, if (is.null(t_s)) 0L else as.integer(t_s[["n.pool"]]),
            if (is.null(t_s) || !is.finite(t_s[["lo"]])) "  -  " else sprintf("%.4f", t_s[["lo"]]),
            if (is.null(t_s) || !is.finite(t_s[["hi"]])) "  -  " else sprintf("%.4f", t_s[["hi"]]),
            nf, 100 * nf / max(1, sum(!is.na(Os))),
            sum(per_site > 0), length(k), sum(per_samp > 0),
            if (length(mg)) median(mg) else NA_real_))
        thr_rows[[length(thr_rows) + 1]] <- data.frame(
          tissue = tissue, scope = scope, variant = vn, q = Q_MAIN,
          methy.state = s, tail = TAIL[[s]],
          sites = length(k),
          pool.n = if (is.null(t_s)) NA_integer_ else as.integer(t_s[["n.pool"]]),
          thr.lower = if (is.null(t_s)) NA_real_ else round(t_s[["lo"]], 4),
          thr.upper = if (is.null(t_s)) NA_real_ else round(t_s[["hi"]], 4),
          flags = nf,
          flag.rate.pct = round(100 * nf / max(1, sum(!is.na(Os))), 3),
          sites.with.flag = sum(per_site > 0),
          pct.sites.with.flag = round(100 * sum(per_site > 0) / length(k), 1),
          max.per.site = if (nf) max(per_site) else 0L,
          gini.sites = round(gini(per_site), 3),
          samples.with.flag = sum(per_samp > 0),
          max.per.sample = if (nf) max(per_samp) else 0L,
          gini.samples = round(gini(per_samp), 3),
          median.abs.dbeta = if (length(mg)) round(median(mg), 4) else NA_real_,
          pct.under.0.05 = if (length(mg)) round(100 * mean(mg < 0.05), 1) else NA_real_,
          pct.under.0.10 = if (length(mg)) round(100 * mean(mg < 0.10), 1) else NA_real_)
      }
    }

    # ---- percentile sweep: "or lower a bit" ----
    for (vn in names(variants)) for (q in QGRID) {
      O <- state_pooled_flags(variants[[vn]], st, q)
      mg <- dev[!is.na(O) & O != 0]
      sweep_rows[[length(sweep_rows) + 1]] <- data.frame(
        tissue = tissue, scope = scope, variant = vn, q = q,
        flags = sum(O != 0, na.rm = TRUE),
        flag.rate.pct = round(100 * mean(O != 0, na.rm = TRUE), 3),
        sites.with.flag = sum(rowSums(O != 0, na.rm = TRUE) > 0),
        pct.sites.with.flag = round(100 * mean(rowSums(O != 0, na.rm = TRUE) > 0), 1),
        samples.with.flag = sum(colSums(O != 0, na.rm = TRUE) > 0),
        median.abs.dbeta = if (length(mg)) round(median(mg), 4) else NA_real_,
        pct.under.0.05 = if (length(mg)) round(100 * mean(mg < 0.05), 1) else NA_real_,
        pct.under.0.10 = if (length(mg)) round(100 * mean(mg < 0.10), 1) else NA_real_,
        jaccard.vs.ext = round(jaccard(O, F), 4))
    }

    # ---- headline scoring at Q_MAIN, incl. stability ----
    for (vn in names(variants)) {
      O <- state_pooled_flags(variants[[vn]], st, Q_MAIN)
      s_acc <- 0
      for (r in seq_len(N_REPS)) {
        Bn <- B + matrix(rnorm(length(B), 0, NOISE_SD), nrow(B))
        Bn[Bn < 0] <- 0; Bn[Bn > 1] <- 1; Bn[is.na(B)] <- NA
        Xn <- if (vn == "bio.stat.abs") Bn else Bn - apply(Bn, 1, median, na.rm = TRUE)
        s_acc <- s_acc + jaccard(O, state_pooled_flags(Xn, st, Q_MAIN))
      }
      mg <- dev[!is.na(O) & O != 0]
      rate_by_state <- sapply(STATE_ORDER, function(s) {
        k <- which(st == s); if (!length(k)) return(NA_real_)
        mean(O[k, ] != 0, na.rm = TRUE) })
      mag_by_state <- sapply(STATE_ORDER, function(s) {
        k <- which(st == s); if (!length(k)) return(NA_real_)
        v <- as.vector(O[k, ]); w <- as.vector(dev[k, ])
        ok <- !is.na(v) & v != 0 & !is.na(w)
        if (!any(ok)) return(NA_real_); median(w[ok]) })
      score_rows[[length(score_rows) + 1]] <- data.frame(
        tissue = tissue, scope = scope, method = vn, q = Q_MAIN,
        flag.rate.pct = round(100 * mean(O != 0, na.rm = TRUE), 3),
        stability.jaccard = round(s_acc / N_REPS, 4),
        state.rate.ratio = round(max(rate_by_state, na.rm = TRUE) /
                                 max(1e-9, min(rate_by_state, na.rm = TRUE)), 2),
        magnitude.H.vs.R = round(mag_by_state[["H"]] / mag_by_state[["R"]], 3),
        median.abs.dbeta = if (length(mg)) round(median(mg), 4) else NA_real_,
        pct.trivial = if (length(mg)) round(100 * mean(mg < 0.10), 2) else NA_real_,
        jaccard.vs.ext = round(jaccard(O, F), 4))
    }

    # ---- the window's L sites, spelled out, since that is the ask ----
    if (scope == "window") {
      for (vn in names(variants)) {
        O <- state_pooled_flags(variants[[vn]], st, Q_MAIN)
        for (s in STATE_ORDER) {
          k <- which(st == s); if (!length(k)) next
          Os <- O[k, , drop = FALSE]
          site_rows[[length(site_rows) + 1]] <- data.frame(
            tissue = tissue, variant = vn, q = Q_MAIN, methy.state = s,
            cgID = rownames(Os), pos = pos[k],
            cohort.median.beta = round(med[k], 4),
            biostat.neg1 = rowSums(Os == -1, na.rm = TRUE),
            biostat.pos1 = rowSums(Os ==  1, na.rm = TRUE),
            ext.neg1 = rowSums(F[k, , drop = FALSE] == -1, na.rm = TRUE),
            ext.pos1 = rowSums(F[k, , drop = FALSE] ==  1, na.rm = TRUE))
          sample_rows[[length(sample_rows) + 1]] <- data.frame(
            tissue = tissue, variant = vn, q = Q_MAIN, methy.state = s,
            sample = colnames(Os),
            biostat.flags = colSums(Os != 0, na.rm = TRUE),
            ext.flags = colSums(F[k, , drop = FALSE] != 0, na.rm = TRUE))
        }
      }
    }
    cat("\n")
  }
  rm(B_all, F_all); invisible(gc())
}

write_out(do.call(rbind, thr_rows),   "Task30_StatePooled_Thresholds.csv")
write_out(do.call(rbind, sweep_rows), "Task30_StatePooled_QSweep.csv")
write_out(do.call(rbind, score_rows), "Task30_StatePooled_Scores.csv")
write_out(do.call(rbind, site_rows),  "Task30_StatePooled_Window_BySite.csv")
write_out(do.call(rbind, sample_rows),"Task30_StatePooled_Window_BySample.csv")

cat("\n---- headline scores (compare against Task24_MethodBenchmark_Scores.csv) ----\n")
print(do.call(rbind, score_rows), row.names = FALSE)

cat("\n---- percentile sweep, chr22 scope ----\n")
sw <- do.call(rbind, sweep_rows)
print(sw[sw$scope == "chr22", ], row.names = FALSE)

cat("\n=========================================================\n")
cat("TASK 30 COMPLETE\n")
cat("=========================================================\n")
