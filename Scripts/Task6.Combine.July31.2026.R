# task6_combine.R
# Combine BRCA 53 Alive + 32 Dead cohorts (Normal, Tumor) and run real OutlierMeth

library(OutlierMeth)

cat("Loading all four datasets...\n")
AN <- read.table("Sorted.BRCA.53Alive.Normal.380355cg.75col.May28.2026.txt", header = TRUE)
AT <- read.table("Sorted.BRCA.53Alive.Tumor.380355cg.75col.May28.2026.txt", header = TRUE)
DN <- read.table("Sorted.BRCA.32Dead.Normal.380355cg.54col.May28.2026.txt", header = TRUE)
DT <- read.table("Sorted.BRCA.32Dead.Tumor.380355cg.54col.May28.2026.txt", header = TRUE)
cat("Finished loading.\n")

alive.normal <- as.matrix(AN[, 5:57])   # 53 samples
alive.tumor  <- as.matrix(AT[, 5:57])
dead.normal  <- as.matrix(DN[, 5:36])   # 32 samples
dead.tumor   <- as.matrix(DT[, 5:36])

rownames(alive.normal) <- AN$Composite.Element.REF
rownames(alive.tumor)  <- AT$Composite.Element.REF
rownames(dead.normal)  <- DN$Composite.Element.REF
rownames(dead.tumor)   <- DT$Composite.Element.REF

cat("Row order check (should all be TRUE):\n")
cat(identical(rownames(alive.normal), rownames(dead.normal)), "\n")
cat(identical(rownames(alive.tumor),  rownames(dead.tumor)),  "\n")

combined.normal <- cbind(alive.normal, dead.normal)   # 85 samples
combined.tumor  <- cbind(alive.tumor,  dead.tumor)    # 85 samples

cat("\nCombined dimensions:\n")
print(dim(combined.normal))
print(dim(combined.tumor))

run_full <- function(beta, label) {
  t0 <- Sys.time()
  ref  <- referenceMeth(beta)
  flag <- flagMeth(beta[rownames(ref), ], reference = ref, p = 0.01)
  cat(label, "referenceMeth+flagMeth done in",
      round(difftime(Sys.time(), t0, units = "mins"), 2), "min\n")
  list(per.cpg = rowSums(abs(flag), na.rm = TRUE),
       per.sample = colSums(abs(flag), na.rm = TRUE))
}

cat("\nRunning combined Normal (85 samples, full genome) ...\n")
res.normal <- run_full(combined.normal, "Combined Normal")

cat("\nRunning combined Tumor (85 samples, full genome) ...\n")
res.tumor <- run_full(combined.tumor, "Combined Tumor")

cat("\n===== Per-sample outlier CpG counts: Combined Normal =====\n")
print(sort(res.normal$per.sample, decreasing = TRUE))

cat("\n===== Per-sample outlier CpG counts: Combined Tumor =====\n")
print(sort(res.tumor$per.sample, decreasing = TRUE))

cat("\nMean per-CpG outlier count -> Normal:", mean(res.normal$per.cpg),
    " Tumor:", mean(res.tumor$per.cpg), "\n")

cat("\nDONE.\n")
