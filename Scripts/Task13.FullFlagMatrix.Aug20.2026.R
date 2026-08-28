# ================================================================
# Task 13 — Full outlier flag matrix (380,355 CpGs x 65 columns)
# Self-referential AND external-reference (tcga.rda) versions,
# Normal and Tumor run separately -> 4 output CSVs total.
# ================================================================
#
# Output columns (65 total = 12 metadata + 53 samples):
#   cgID, chr, pos, methy.state, outliers.coef2, outliers.coef3,
#   flag.count.neg1, flag.count.0, flag.count.1,
#   flag.pct.neg1,   flag.pct.0,   flag.pct.1,
#   [53 individual sample flag columns: N1..N53 or T1..T53]
#
# ASSUMPTIONS (flag if wrong, easy to change):
#  - "pos" = Start coordinate only. If you want End too, say so and
#    this becomes a 66th column.
#  - "without using references" = self-referential (own 53 samples
#    build their own threshold table via referenceMeth()).
#    "with reference" = external tcga.rda panel (~2,000 independent
#    samples), loaded pre-built, not re-derived from our data.
#  - Row order is preserved throughout via positional alignment, not
#    a join: beta's rownames are set directly from meta_df's
#    Composite.Element.REF in original file order, and both flagMeth()
#    and the self-ref padding step below preserve that same row order.
#    No merge/join is used anywhere, so there's no risk of silent
#    row-mismatch.
#
# NOTE on external-reference rows: tcga.rda covers 370,201 of our
# 380,355 CpGs (per Task 12). For the ~10,154 unmatched CpGs, flagMeth()
# cannot compute a flag -> all 53 sample columns are NA for that row,
# and flag.count.*/flag.pct.* are also left NA (not silently 0) so an
# unmatched CpG can never be misread as "flagged normal for everyone."
#
# NOTE on the "Rc" bug from Task10.FlagStateTable.Aug6.2026.R:
# that script's state_order vector had "Rc" instead of "R", so the R
# methylation-state (the single largest category in Tumor, ~40% of
# CpGs) was silently dropped from that table. This script does not
# filter or reorder methy.state at all -- all 6 states pass through
# untouched -- so that bug cannot reoccur here.

library(OutlierMeth)

setwd("/home/s_s355/research3.data/TCGA/filter.BRCA.UCEC.Aug17.2022/sorted.TCGA.380355cg.files")

out_dir <- "/mmfs1/home/wln26/Experiments.Outlier.July31.2026/Results"

# ---- optional fast CSV writer; falls back to base write.csv if
#      data.table isn't installed (avoids a repeat of the earlier
#      stringi/farver install headaches for something this script
#      doesn't strictly need) ----
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

cat("Loading external TCGA reference panel...\n")
load("/mmfs1/home/wln26/OutlierMeth_tarball_check/OutlierMeth/data/tcga.rda")

# ---- build the 65-column table for one flag matrix ----
build_flag_table <- function(flag_matrix, meta_df) {

  n_samples <- ncol(flag_matrix)

  all_na   <- apply(flag_matrix, 1, function(r) all(is.na(r)))
  cnt.neg1 <- apply(flag_matrix, 1, function(r) sum(r == -1, na.rm = TRUE))
  cnt.0    <- apply(flag_matrix, 1, function(r) sum(r ==  0, na.rm = TRUE))
  cnt.1    <- apply(flag_matrix, 1, function(r) sum(r ==  1, na.rm = TRUE))
  cnt.neg1[all_na] <- NA
  cnt.0[all_na]    <- NA
  cnt.1[all_na]    <- NA

  pct.neg1 <- round(cnt.neg1 / n_samples * 100, 3)
  pct.0    <- round(cnt.0    / n_samples * 100, 3)
  pct.1    <- round(cnt.1    / n_samples * 100, 3)

  out <- data.frame(
    cgID            = meta_df$Composite.Element.REF,
    chr             = meta_df$Chromosome,
    pos             = meta_df$Start,
    methy.state     = meta_df$methy.state,
    outliers.coef2  = meta_df$outliers.coef2,
    outliers.coef3  = meta_df$outliers.coef3,
    flag.count.neg1 = cnt.neg1,
    flag.count.0    = cnt.0,
    flag.count.1    = cnt.1,
    flag.pct.neg1   = pct.neg1,
    flag.pct.0      = pct.0,
    flag.pct.1      = pct.1,
    stringsAsFactors = FALSE
  )

  cbind(out, as.data.frame(flag_matrix))
}

run_and_write <- function(beta, meta_df, label, ref_mode) {

  cat("\n===", label, "-", ref_mode, "reference ===\n")
  t0 <- Sys.time()

  if (ref_mode == "self") {
    ref  <- referenceMeth(beta)
    flag <- flagMeth(beta[rownames(ref), ], reference = ref, p = 0.01)
    # referenceMeth()'s na.omit() can drop a handful of CpG rows; pad
    # back to the full 380,355 so every output file has the same row
    # count/order as the other three
    full_flag <- matrix(NA, nrow = nrow(beta), ncol = ncol(beta),
                         dimnames = list(rownames(beta), colnames(beta)))
    full_flag[rownames(flag), ] <- flag
    flag <- full_flag
  } else {
    # tcga is already a pre-built ~2,000-sample reference; flagMeth()
    # returns a matrix with the same dimnames as beta, in beta's order
    flag <- flagMeth(beta, reference = tcga, p = 0.01)
  }

  tbl <- build_flag_table(flag, meta_df)

  fname <- file.path(out_dir, paste0("Task13_", label, "_", ref_mode, "ref_flags.csv"))
  write_out(tbl, fname)

  cat(label, ref_mode, "reference done in",
      round(difftime(Sys.time(), t0, units = "mins"), 2), "min. Written to:\n ", fname, "\n")
  invisible(tbl)
}

# ---- Normal ----
run_and_write(normal.beta, BRCA.53AN, "Normal", "self")
run_and_write(normal.beta, BRCA.53AN, "Normal", "ext")

# ---- Tumor ----
run_and_write(tumor.beta,  BRCA.53AT, "Tumor",  "self")
run_and_write(tumor.beta,  BRCA.53AT, "Tumor",  "ext")

cat("\nDONE. 4 files written to", out_dir, ":\n")
cat(" Task13_Normal_selfref_flags.csv\n")
cat(" Task13_Normal_extref_flags.csv\n")
cat(" Task13_Tumor_selfref_flags.csv\n")
cat(" Task13_Tumor_extref_flags.csv\n")
