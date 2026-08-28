# ================================================================
# Task 15 — Flag bias check: by methylation state and by sample
# Self-referential only (no external reference)
#
# Addresses To-Do items 1(1) and 1(2) from the Aug 6, 2026 meeting notes:
#   1(1) against methylation state -- do more -1/+1 flags land in the
#        L and H states? check using both MEAN and MEDIAN to see if
#        they disagree (a mean/median gap would itself suggest a
#        skewed, non-representative summary).
#   1(2) against sample -- do samples 14-17 stand out as unusual
#        relative to the rest of the 53-sample cohort? bar plots with
#        14-17 highlighted.
#
# Run separately for Normal and Tumor.
#
# NOTE: methy.state is passed through untouched -- no filtering or
# reordering -- so the "Rc" bug from Task10.FlagStateTable.Aug6.2026.R
# (which silently dropped the R state) cannot reoccur here; all 6
# states, including R, are included throughout.
#
# NOTE on row/column alignment: flag_matrix rows are in the same order
# as meta_df (BRCA.53AN / BRCA.53AT) rows -- both trace back to the
# same beta matrix built directly from that data frame -- so Part A
# uses direct positional alignment, no join. Part B operates on flag
# matrix columns (samples), which keep their original N1..N53 /
# T1..T53 names from the source file.
# ================================================================

library(OutlierMeth)
library(dplyr)
library(ggplot2)

setwd("/home/s_s355/research3.data/TCGA/filter.BRCA.UCEC.Aug17.2022/sorted.TCGA.380355cg.files")
out_dir <- "/mmfs1/home/wln26/Experiments.Outlier.July31.2026/Results"

write_out <- function(tbl, fname) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fwrite(tbl, fname)
  } else {
    write.csv(tbl, fname, row.names = FALSE)
  }
}

cat("Loading BRCA 53 Alive datasets...\n")
BRCA.53AN <- read.table("Sorted.BRCA.53Alive.Normal.380355cg.75col.May28.2026.txt", header = TRUE)
BRCA.53AT <- read.table("Sorted.BRCA.53Alive.Tumor.380355cg.75col.May28.2026.txt",  header = TRUE)
cat("Finished loading.\n")

normal.beta <- as.matrix(BRCA.53AN[, 5:57])
tumor.beta  <- as.matrix(BRCA.53AT[, 5:57])
rownames(normal.beta) <- BRCA.53AN$Composite.Element.REF
rownames(tumor.beta)  <- BRCA.53AT$Composite.Element.REF

# ---- self-referential flag matrix, padded back to the full row count
#      so it lines up positionally with meta_df ----
self_flag <- function(beta) {
  ref  <- referenceMeth(beta)
  flag <- flagMeth(beta[rownames(ref), ], reference = ref, p = 0.01)
  full_flag <- matrix(NA, nrow = nrow(beta), ncol = ncol(beta),
                       dimnames = list(rownames(beta), colnames(beta)))
  full_flag[rownames(flag), ] <- flag
  full_flag
}

cat("\nBuilding self-referential flags (Normal)...\n")
flag.normal <- self_flag(normal.beta)
cat("Building self-referential flags (Tumor)...\n")
flag.tumor  <- self_flag(tumor.beta)

# ================================================================
# PART A -- against methylation state (mean AND median)
# ================================================================
# Vectorized row-wise counts (no pivot_longer) -- the flag matrix is
# 380,355 x 53 (~20M cells), so reshaping to long format first would
# be needlessly slow/memory-heavy for what is just a per-row tally.

state_summary <- function(flag_matrix, meta_df, label) {

  cnt.neg1 <- apply(flag_matrix, 1, function(r) sum(r == -1, na.rm = TRUE))
  cnt.0    <- apply(flag_matrix, 1, function(r) sum(r ==  0, na.rm = TRUE))
  cnt.1    <- apply(flag_matrix, 1, function(r) sum(r ==  1, na.rm = TRUE))

  per_cpg <- data.frame(
    methy.state = meta_df$methy.state,
    n.neg1 = cnt.neg1,
    n.0    = cnt.0,
    n.1    = cnt.1
  )

  by_state <- per_cpg %>%
    group_by(methy.state) %>%
    summarise(
      n.cpgs      = n(),
      mean.neg1   = mean(n.neg1, na.rm = TRUE),
      median.neg1 = median(n.neg1, na.rm = TRUE),
      mean.0      = mean(n.0, na.rm = TRUE),
      median.0    = median(n.0, na.rm = TRUE),
      mean.1      = mean(n.1, na.rm = TRUE),
      median.1    = median(n.1, na.rm = TRUE),
      .groups = "drop"
    )

  fname <- file.path(out_dir, paste0("Task15_", label, "_flagbias_by_state.csv"))
  write_out(by_state, fname)
  cat(label, "by-state summary written to:", fname, "\n")
  print(as.data.frame(by_state))

  by_state
}

cat("\n=== Part A: flag bias by methylation state ===\n")
state.normal <- state_summary(flag.normal, BRCA.53AN, "Normal")
state.tumor  <- state_summary(flag.tumor,  BRCA.53AT, "Tumor")

# ================================================================
# PART B -- against sample, highlighting 14-17
# ================================================================

sample_summary <- function(flag_matrix, label) {

  n_cpgs   <- apply(flag_matrix, 2, function(col) sum(!is.na(col)))
  cnt.neg1 <- apply(flag_matrix, 2, function(col) sum(col == -1, na.rm = TRUE))
  cnt.0    <- apply(flag_matrix, 2, function(col) sum(col ==  0, na.rm = TRUE))
  cnt.1    <- apply(flag_matrix, 2, function(col) sum(col ==  1, na.rm = TRUE))

  out <- data.frame(
    sample     = colnames(flag_matrix),
    n.cpgs     = n_cpgs,
    count.neg1 = cnt.neg1,
    count.0    = cnt.0,
    count.1    = cnt.1,
    pct.neg1   = round(cnt.neg1 / n_cpgs * 100, 3),
    pct.0      = round(cnt.0    / n_cpgs * 100, 3),
    pct.1      = round(cnt.1    / n_cpgs * 100, 3),
    stringsAsFactors = FALSE
  )

  out$sample.num <- as.integer(gsub("[^0-9]", "", out$sample))
  out$highlight  <- out$sample.num %in% 14:17

  fname <- file.path(out_dir, paste0("Task15_", label, "_flagbias_by_sample.csv"))
  write_out(out, fname)
  cat(label, "by-sample summary written to:", fname, "\n")

  out
}

cat("\n=== Part B: flag bias by sample (highlighting 14-17) ===\n")
sample.normal <- sample_summary(flag.normal, "Normal")
sample.tumor  <- sample_summary(flag.tumor,  "Tumor")

# ---- bar plots: total flagged CpGs (|count.neg1| + |count.1|) per
#      sample, samples 14-17 highlighted ----
make_sample_barplot <- function(sample_df, label) {
  sample_df <- sample_df %>%
    mutate(total.flagged = count.neg1 + count.1,
           sample = factor(sample, levels = sample[order(sample.num)]))

  ggplot(sample_df, aes(x = sample, y = total.flagged, fill = highlight)) +
    geom_col() +
    scale_fill_manual(values = c(`FALSE` = "grey70", `TRUE` = "firebrick"),
                       labels = c("Other samples", "Samples 14-17"),
                       name = NULL) +
    labs(title = paste0(label, " -- total flagged CpGs per sample (self-referential)"),
         x = "Sample", y = "Total flagged CpGs (count.neg1 + count.1)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7))
}

p.normal <- make_sample_barplot(sample.normal, "Normal")
p.tumor  <- make_sample_barplot(sample.tumor,  "Tumor")

ggsave(file.path(out_dir, "Task15_Normal_sample_barplot.pdf"), p.normal, width = 11, height = 6)
ggsave(file.path(out_dir, "Task15_Tumor_sample_barplot.pdf"),  p.tumor,  width = 11, height = 6)

cat("\nDONE.\n")
cat("Part A CSVs:  Task15_Normal_flagbias_by_state.csv, Task15_Tumor_flagbias_by_state.csv\n")
cat("Part B CSVs:  Task15_Normal_flagbias_by_sample.csv, Task15_Tumor_flagbias_by_sample.csv\n")
cat("Part B plots: Task15_Normal_sample_barplot.pdf, Task15_Tumor_sample_barplot.pdf\n")
