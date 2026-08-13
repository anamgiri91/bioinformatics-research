# experiment7_remove_outlier_samples.R
# Uses the real OutlierMeth package (now installed) to:
#  1. Build a reference table per CpG (referenceMeth) for Normal and Tumor
#  2. Flag each CpG x sample cell as hyper/hypo outlier (flagMeth)
#  3. Sum flags per sample -> outlier CpG count per sample
#  4. Flag and remove outlier SAMPLES (mean + 2*SD rule)
#  5. Recompute numeric summary before/after removal

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

cat("\nBuilding reference table (Normal) ... this is the slow step, please wait\n")
t0 <- Sys.time()
ref.normal <- referenceMeth(normal.beta)
cat("Done in", round(difftime(Sys.time(), t0, units = "mins"), 2), "minutes\n")
saveRDS(ref.normal, "~/experiment7_ref_normal.rds")

cat("\nBuilding reference table (Tumor) ...\n")
t0 <- Sys.time()
ref.tumor <- referenceMeth(tumor.beta)
cat("Done in", round(difftime(Sys.time(), t0, units = "mins"), 2), "minutes\n")
saveRDS(ref.tumor, "~/experiment7_ref_tumor.rds")

cat("\nFlagging outlier CpGs per sample (Normal) ...\n")
flag.normal <- flagMeth(normal.beta[rownames(ref.normal), ], reference = ref.normal, p = 0.01)

cat("Flagging outlier CpGs per sample (Tumor) ...\n")
flag.tumor <- flagMeth(tumor.beta[rownames(ref.tumor), ], reference = ref.tumor, p = 0.01)

# per-sample outlier CpG counts: sum of |flag| (1 = hyper or hypo) down each column
normal.sample.counts <- colSums(abs(flag.normal), na.rm = TRUE)
tumor.sample.counts  <- colSums(abs(flag.tumor), na.rm = TRUE)

cat("\n===== Per-sample outlier CpG counts: Normal =====\n")
print(sort(normal.sample.counts, decreasing = TRUE))

cat("\n===== Per-sample outlier CpG counts: Tumor =====\n")
print(sort(tumor.sample.counts, decreasing = TRUE))

# flag outlier samples: more than mean + 2*SD outlier CpGs
flag_outlier_samples <- function(counts, k = 2) {
  m <- mean(counts); s <- sd(counts)
  names(counts)[counts > m + k * s]
}

normal.outlier.samples <- flag_outlier_samples(normal.sample.counts)
tumor.outlier.samples  <- flag_outlier_samples(tumor.sample.counts)

cat("\nNormal samples flagged as outliers (mean + 2*SD):\n")
print(normal.outlier.samples)
cat("\nTumor samples flagged as outliers (mean + 2*SD):\n")
print(tumor.outlier.samples)

# remove flagged samples and recompute numeric summary
normal.beta.clean <- normal.beta[, !(colnames(normal.beta) %in% normal.outlier.samples)]
tumor.beta.clean  <- tumor.beta[,  !(colnames(tumor.beta)  %in% tumor.outlier.samples)]

cat("\nDimensions before removal:\n")
print(dim(normal.beta)); print(dim(tumor.beta))
cat("Dimensions after removal:\n")
print(dim(normal.beta.clean)); print(dim(tumor.beta.clean))

numeric.summary.clean <- data.frame(
  Dataset = c("Normal (samples removed)", "Tumor (samples removed)"),
  Mean   = c(mean(normal.beta.clean), mean(tumor.beta.clean)),
  Median = c(median(normal.beta.clean), median(tumor.beta.clean)),
  SD     = c(sd(normal.beta.clean), sd(tumor.beta.clean)),
  Min    = c(min(normal.beta.clean), min(tumor.beta.clean)),
  Max    = c(max(normal.beta.clean), max(tumor.beta.clean))
)

cat("\n===== Numeric Summary After Removing Outlier Samples =====\n")
print(numeric.summary.clean)

cat("\nDONE.\n")
