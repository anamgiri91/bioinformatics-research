library(OutlierMeth)

cat("Loading BRCA 53 Alive datasets...\n")
BRCA.53AN <- read.table("Sorted.BRCA.53Alive.Normal.380355cg.75col.May28.2026.txt", header = TRUE)
BRCA.53AT <- read.table("Sorted.BRCA.53Alive.Tumor.380355cg.75col.May28.2026.txt", header = TRUE)

normal.beta <- as.matrix(BRCA.53AN[, 5:57])
tumor.beta  <- as.matrix(BRCA.53AT[, 5:57])
rownames(normal.beta) <- BRCA.53AN$Composite.Element.REF
rownames(tumor.beta)  <- BRCA.53AT$Composite.Element.REF

cat("Loading external TCGA reference panel...\n")
load("/mmfs1/home/wln26/OutlierMeth_tarball_check/OutlierMeth/data/tcga.rda")

cat("\n=== Self-referential (own 53 samples as reference) ===\n")
ref.self.normal  <- referenceMeth(normal.beta)
flag.self.normal <- flagMeth(normal.beta[rownames(ref.self.normal), ], reference = ref.self.normal, p = 0.01)

ref.self.tumor  <- referenceMeth(tumor.beta)
flag.self.tumor <- flagMeth(tumor.beta[rownames(ref.self.tumor), ], reference = ref.self.tumor, p = 0.01)

cat("Self-referential mean per-CpG count -> Normal:", mean(rowSums(abs(flag.self.normal), na.rm = TRUE)),
    " Tumor:", mean(rowSums(abs(flag.self.tumor), na.rm = TRUE)), "\n")

cat("\n=== External reference (tcga.rda, ~2,000 independent samples) ===\n")
flag.ext.normal <- flagMeth(normal.beta, reference = tcga, p = 0.01)
flag.ext.tumor  <- flagMeth(tumor.beta,  reference = tcga, p = 0.01)

cat("External-reference mean per-CpG count -> Normal:", mean(rowSums(abs(flag.ext.normal), na.rm = TRUE)),
    " Tumor:", mean(rowSums(abs(flag.ext.tumor), na.rm = TRUE)), "\n")

cat("\nDONE.\n")
