library(OutlierMeth)
 
cat("Loading BRCA 53 Alive datasets...\n")
BRCA.53AN <- read.table(
  "Sorted.BRCA.53Alive.Normal.380355cg.75col.May28.2026.txt",
  header = TRUE
)
BRCA.53AT <- read.table(
  "Sorted.BRCA.53Alive.Tumor.380355cg.75col.May28.2026.txt",
  header = TRUE
)
cat("Finished loading.\n")
 
normal.beta <- as.matrix(BRCA.53AN[, 5:57])
tumor.beta  <- as.matrix(BRCA.53AT[, 5:57])
rownames(normal.beta) <- BRCA.53AN$Composite.Element.REF
rownames(tumor.beta)  <- BRCA.53AT$Composite.Element.REF
 
run_full <- function(beta, label) {
  t0 <- Sys.time()
  ref  <- referenceMeth(beta)
  flag <- flagMeth(beta[rownames(ref), ], reference = ref, p = 0.01)
  cat(label, "done in", round(difftime(Sys.time(), t0, units = "mins"), 2), "min\n")
  list(per.cpg = rowSums(abs(flag), na.rm = TRUE),
       per.sample = colSums(abs(flag), na.rm = TRUE))
}
 
cat("\n=== BEFORE removing samples 14-17 ===\n")
before.normal <- run_full(normal.beta, "Normal (before)")
before.tumor  <- run_full(tumor.beta,  "Tumor (before)")
 
# ---- remove sample columns N14/N15/N16/N17 and T14/T15/T16/T17 ----
# NOTE: adjust the "N"/"T" prefixes below if your column names differ
remove.ids          <- 14:17
remove.cols.normal  <- paste0("N", remove.ids)
remove.cols.tumor   <- paste0("T", remove.ids)
 
cat("\nColumns targeted for removal:\n")
cat("  Normal:", paste(remove.cols.normal, collapse = ", "), "\n")
cat("  Tumor :", paste(remove.cols.tumor,  collapse = ", "), "\n")
 
missing.normal <- setdiff(remove.cols.normal, colnames(normal.beta))
missing.tumor  <- setdiff(remove.cols.tumor,  colnames(tumor.beta))
if (length(missing.normal) > 0) cat("WARNING - not found in Normal columns:", missing.normal, "\n")
if (length(missing.tumor)  > 0) cat("WARNING - not found in Tumor columns:",  missing.tumor,  "\n")
 
normal.clean <- normal.beta[, !(colnames(normal.beta) %in% remove.cols.normal)]
tumor.clean  <- tumor.beta[,  !(colnames(tumor.beta)  %in% remove.cols.tumor)]
 
cat("\nDimensions before removal -> Normal:", ncol(normal.beta), " Tumor:", ncol(tumor.beta), "\n")
cat("Dimensions after removal  -> Normal:", ncol(normal.clean), " Tumor:", ncol(tumor.clean), "\n")
 
cat("\n=== AFTER removing samples 14-17 ===\n")
after.normal <- run_full(normal.clean, "Normal (after)")
after.tumor  <- run_full(tumor.clean,  "Tumor (after)")
 
cat("\n===== Mean per-CpG outlier count, before vs after =====\n")
cat("Normal before:", mean(before.normal$per.cpg), " after:", mean(after.normal$per.cpg), "\n")
cat("Tumor before:",  mean(before.tumor$per.cpg),  " after:", mean(after.tumor$per.cpg),  "\n")
 
cat("\nDONE.\n")
