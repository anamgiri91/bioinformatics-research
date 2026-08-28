# ================================================================
# Task 23 - Tasks 10 / 13A / 14 regenerated with the R state restored
# Date: Aug 28, 2026
#
# The bug (REVIEW.md Finding 2): Task10.FlagStateTable.Aug6.2026.R,
# Task13.FlagBias.StateAndSample.NoExtRef.Aug2026.R and
# Task14.FlagBias.StateAndSample.ExtRef.Aug2026.R all declared
#
#   state_order <- c("L", "LM", "M", "HM", "H", "Rc")
#
# The data has "R". The subsequent tab[ord, ] kept only matched
# states, so every R-state CpG was dropped silently - 6.2% of Normal
# and 40.3% of Tumor. Task 14's tumour table is the strongest result
# in the project and was computed with the single largest category of
# tumour CpGs missing.
#
# Those three scripts are now fixed and carry a guard that stops on
# any state present in the data but absent from state_order, so the
# class of bug cannot recur. Rerunning them end to end, however, means
# rebuilding referenceMeth() over 380,355 x 53 (~7 min each) and
# loading tcga.rda, which is not always to hand.
#
# It is also unnecessary. Task 13's four flag matrices already store
# the finished flags for all 380,355 CpGs alongside methy.state, and
# Task13.FullFlagMatrix.Aug20.2026.R applies no state filter at all.
# Everything the three broken tables reported is recoverable from
# them by aggregation alone - no reference rebuild, no tcga.rda, and
# no chance of a different random path through the pipeline.
#
# This script therefore reproduces the corrected tables from those
# matrices, and prints the before/after delta so the size of the bug
# is on the record rather than merely asserted.
#
# Outputs (per tissue x reference):
#   Task23_StateFlagTable_<ref>_<tissue>.csv     counts + row %
#   Task23_StateBiasMeanMedian_<ref>_<tissue>.csv
#   Task23_SampleFlagSummary_<ref>_<tissue>.csv
# plus:
#   Task23_RStateImpact.csv    what the "Rc" bug removed
# ================================================================

options(stringsAsFactors = FALSE)

pick_dir <- function(...) {
  cand <- c(...)
  for (d in cand) if (dir.exists(d)) return(d)
  stop("no candidate directory exists:\n  ", paste(cand, collapse = "\n  "))
}
repo_dir <- pick_dir("/mmfs1/home/wln26/Experiments.Outlier.July31.2026",
                     path.expand("~/Desktop/bioinformatics-research"))
out_dir <- file.path(repo_dir, "Results")

STATE_ORDER <- c("L", "LM", "M", "HM", "H", "R")
N_CPG_EXPECTED <- 380355

flag_files <- list(
  selfref = list(Normal = "Task13_Normal_selfref_flags.csv",
                 Tumor  = "Task13_Tumor_selfref_flags.csv"),
  extref  = list(Normal = "Task13_Normal_extref_flags.csv",
                 Tumor  = "Task13_Tumor_extref_flags.csv"))

has_dt <- requireNamespace("data.table", quietly = TRUE)
read_big <- function(p) {
  if (has_dt) as.data.frame(data.table::fread(p, showProgress = FALSE))
  else read.csv(p, check.names = FALSE)
}
write_out <- function(tbl, fname) {
  path <- file.path(out_dir, fname)
  write.csv(tbl, path, row.names = FALSE)
  cat("  written:", path, " (", nrow(tbl), " rows )\n", sep = "")
}

cat("====================================================\n")
cat("TASK 23 - CORRECTED STATE TABLES (R RESTORED)\n")
cat("====================================================\n\n")

impact <- list()

for (ref in names(flag_files)) {
  for (tissue in c("Normal", "Tumor")) {

    cat("##################################################\n")
    cat("### ", ref, " / ", tissue, "\n", sep = "")
    cat("##################################################\n")

    fl <- read_big(file.path(out_dir, flag_files[[ref]][[tissue]]))
    cols <- paste0(substr(tissue, 1, 1), 1:53)

    # ---- the guard the three original scripts now carry ----
    present <- unique(stats::na.omit(fl$methy.state))
    missing <- setdiff(present, STATE_ORDER)
    if (length(missing))
      stop("methy.state values absent from STATE_ORDER: ",
           paste(missing, collapse = ", "))
    stopifnot(nrow(fl) == N_CPG_EXPECTED)

    st <- factor(fl$methy.state, levels = STATE_ORDER)
    cat("  CpGs:", nrow(fl), " - state distribution:\n")
    print(table(st))

    F <- as.matrix(fl[, cols, drop = FALSE]); storage.mode(F) <- "numeric"

    # ---- state x flag contingency, all six states ----
    rows <- list()
    for (s in STATE_ORDER) {
      k <- which(st == s)
      if (!length(k)) next
      v <- as.vector(F[k, ])
      nobs <- sum(!is.na(v))
      c1 <- sum(v == -1, na.rm = TRUE); c0 <- sum(v == 0, na.rm = TRUE)
      c2 <- sum(v ==  1, na.rm = TRUE)
      pc <- function(x) if (nobs > 0) round(100 * x / nobs, 3) else NA
      rows[[length(rows) + 1]] <- data.frame(
        reference = ref, tissue = tissue, methy.state = s,
        n.cpgs = length(k), n.cells = length(v),
        n.evaluated = nobs, n.NA = sum(is.na(v)),
        count.neg1 = c1, count.0 = c0, count.pos1 = c2,
        pct.neg1 = pc(c1), pct.0 = pc(c0), pct.pos1 = pc(c2))
    }
    tab <- do.call(rbind, rows)
    print(tab[, c("methy.state", "n.cpgs", "n.evaluated", "count.neg1",
                  "count.0", "count.pos1", "pct.neg1", "pct.0", "pct.pos1")],
          row.names = FALSE)
    write_out(tab, sprintf("Task23_StateFlagTable_%s_%s.csv", ref, tissue))

    # ---- what the bug removed ----
    rrow <- tab[tab$methy.state == "R", ]
    if (nrow(rrow)) {
      impact[[length(impact) + 1]] <- data.frame(
        reference = ref, tissue = tissue,
        R.cpgs = rrow$n.cpgs,
        pct.of.all.cpgs = round(100 * rrow$n.cpgs / nrow(fl), 2),
        R.flags = rrow$count.neg1 + rrow$count.pos1,
        all.flags = sum(tab$count.neg1 + tab$count.pos1),
        pct.of.all.flags = round(100 * (rrow$count.neg1 + rrow$count.pos1) /
                                 sum(tab$count.neg1 + tab$count.pos1), 2))
    }

    # ---- per-CpG rates, then mean/median by state ----
    rate.neg1 <- rowMeans(F == -1, na.rm = TRUE)
    rate.pos1 <- rowMeans(F ==  1, na.rm = TRUE)
    rate.any  <- rowMeans(F !=  0, na.rm = TRUE)
    bias <- do.call(rbind, lapply(STATE_ORDER, function(s) {
      k <- which(st == s)
      if (!length(k)) return(NULL)
      data.frame(reference = ref, tissue = tissue, methy.state = s,
        n.cpgs = length(k),
        mean.rate.neg1 = round(mean(rate.neg1[k], na.rm = TRUE), 5),
        median.rate.neg1 = round(median(rate.neg1[k], na.rm = TRUE), 5),
        mean.rate.pos1 = round(mean(rate.pos1[k], na.rm = TRUE), 5),
        median.rate.pos1 = round(median(rate.pos1[k], na.rm = TRUE), 5),
        mean.rate.any = round(mean(rate.any[k], na.rm = TRUE), 5),
        median.rate.any = round(median(rate.any[k], na.rm = TRUE), 5))
    }))
    cat("\n  per-CpG flag rate by state:\n")
    print(bias[, c("methy.state", "n.cpgs", "mean.rate.neg1",
                   "mean.rate.pos1", "mean.rate.any", "median.rate.any")],
          row.names = FALSE)
    write_out(bias, sprintf("Task23_StateBiasMeanMedian_%s_%s.csv", ref, tissue))

    # ---- per-sample summary ----
    smp <- do.call(rbind, lapply(seq_along(cols), function(j) {
      v <- F[, j]
      data.frame(reference = ref, tissue = tissue, sample = cols[j],
        n.cpgs = nrow(F),
        count.neg1 = sum(v == -1, na.rm = TRUE),
        count.0 = sum(v == 0, na.rm = TRUE),
        count.pos1 = sum(v == 1, na.rm = TRUE),
        n.NA = sum(is.na(v)))
    }))
    smp$total.flags <- smp$count.neg1 + smp$count.pos1
    smp$pct.flagged <- round(100 * smp$total.flags /
                             (smp$n.cpgs - smp$n.NA), 3)
    smp <- smp[order(-smp$total.flags), ]
    rownames(smp) <- NULL
    cat("\n  top 5 samples by flag burden:\n")
    print(head(smp[, c("sample", "count.neg1", "count.pos1",
                       "total.flags", "pct.flagged")], 5), row.names = FALSE)
    write_out(smp, sprintf("Task23_SampleFlagSummary_%s_%s.csv", ref, tissue))

    rm(F, fl); invisible(gc())
    cat("\n")
  }
}

cat("==================================================\n")
cat("WHAT THE \"Rc\" BUG REMOVED\n")
cat("==================================================\n")
imp <- do.call(rbind, impact)
print(imp, row.names = FALSE)
write_out(imp, "Task23_RStateImpact.csv")

cat("\n====================================================\n")
cat("TASK 23 COMPLETE\n")
cat("====================================================\n")
