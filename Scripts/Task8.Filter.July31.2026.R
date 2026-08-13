# task8_filter.R
# Remove low-variance H/L methylation-state CpGs from BRCA 53 Alive, rerun real OutlierMeth
# (compares full-genome vs H/L-filtered, per Dr. Sun's note)

library(OutlierMeth)

cat("Loading BRCA 53 Alive datasets...\n")
AN <- read.table("Sorted.BRCA.53Alive.Normal.380355cg.75col.May28.2026.txt", header = TRUE)
AT <- read.table("Sorted.BRCA.53Alive.Tumor.380355cg.75col.May28.2026.txt", header = TRUE)
cat("Finished loading.\n")

normal.beta <- as.matrix(AN[, 5:57])
tumor.beta  <- as.matrix(AT[, 5:57])
rownames(normal.beta) <- AN$Composite.Element.REF
rownames(tumor.beta)  <- AT$Composite.Element.REF

cat("\nOriginal rows:", nrow(normal.beta), "\n")

run_full <- function(beta, label) {
  t0 <- Sys.time()
  ref  <- referenceMeth(beta)
  flag <- flagMeth(beta[rownames(ref), ], reference = ref, p = 0.01)
  cat(label, "done in", round(difftime(Sys.time(), t0, units = "mins"), 2), "min\n")
  list(per.cpg = rowSums(abs(flag), na.rm = TRUE),
       per.sample = colSums(abs(flag), na.rm = TRUE))
}

cat("\n=== BEFORE filtering (full genome) ===\n")
before.normal <- run_full(normal.beta, "Normal (before)")
before.tumor  <- run_full(tumor.beta,  "Tumor (before)")

# ---- filter out H/L methylation-state sites using existing methy.state column ----
keep.normal <- !(AN$methy.state %in% c("H", "L"))
keep.tumor  <- !(AT$methy.state %in% c("H", "L"))

normal.filtered <- normal.beta[keep.normal[match(rownames(normal.beta), AN$Composite.Element.REF)], ]
tumor.filtered  <- tumor.beta[keep.tumor[match(rownames(tumor.beta), AT$Composite.Element.REF)], ]

cat("\nFiltered rows -> Normal:", nrow(normal.filtered), " Tumor:", nrow(tumor.filtered), "\n")

cat("\n=== AFTER filtering (H/L removed) ===\n")
after.normal <- run_full(normal.filtered, "Normal (after)")
after.tumor  <- run_full(tumor.filtered,  "Tumor (after)")

cat("\n===== Mean per-CpG outlier count, before vs after =====\n")
cat("Normal before:", mean(before.normal$per.cpg), " after:", mean(after.normal$per.cpg), "\n")
cat("Tumor before:",  mean(before.tumor$per.cpg),  " after:", mean(after.tumor$per.cpg),  "\n")

cat("\nDONE.\n")
