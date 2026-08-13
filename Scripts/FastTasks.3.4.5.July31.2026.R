# fast_tasks_3_4_5.R — chrX/chr22 (Alive + Dead) and 5-CpG printout, all via row subsetting

library(OutlierMeth)

run_chrom <- function(beta, refnames, label) {
  t0 <- Sys.time()
  ref  <- referenceMeth(beta)
  flag <- flagMeth(beta[rownames(ref), ], reference = ref, p = 0.01)
  cat(label, "done in", round(difftime(Sys.time(), t0, units = "secs"), 1), "sec\n")
  list(per.cpg = rowSums(abs(flag), na.rm = TRUE),
       per.sample = colSums(abs(flag), na.rm = TRUE))
}

# ---- Task 3: BRCA 53 Alive chrX/chr22 ----
cat("Loading BRCA 53 Alive...\n")
AN <- read.table("Sorted.BRCA.53Alive.Normal.380355cg.75col.May28.2026.txt", header = TRUE)
AT <- read.table("Sorted.BRCA.53Alive.Tumor.380355cg.75col.May28.2026.txt", header = TRUE)

mk <- function(df, chrom, cols) {
  m <- as.matrix(df[df$Chromosome == chrom, cols])
  rownames(m) <- df$Composite.Element.REF[df$Chromosome == chrom]
  m
}

cat("\n=== Task 3: Alive chrX/chr22 ===\n")
r1 <- run_chrom(mk(AN, "chrX",  5:57), NULL, "Alive chrX Normal")
r2 <- run_chrom(mk(AN, "chr22", 5:57), NULL, "Alive chr22 Normal")
r3 <- run_chrom(mk(AT, "chrX",  5:57), NULL, "Alive chrX Tumor")
r4 <- run_chrom(mk(AT, "chr22", 5:57), NULL, "Alive chr22 Tumor")
cat("Mean per-CpG outlier count -> chrX N:", mean(r1$per.cpg), " chr22 N:", mean(r2$per.cpg),
    " chrX T:", mean(r3$per.cpg), " chr22 T:", mean(r4$per.cpg), "\n")

# ---- Task 5: BRCA 32 Dead chrX/chr22 ----
cat("\nLoading BRCA 32 Dead...\n")
DN <- read.table("Sorted.BRCA.32Dead.Normal.380355cg.54col.May28.2026.txt", header = TRUE)
DT <- read.table("Sorted.BRCA.32Dead.Tumor.380355cg.54col.May28.2026.txt", header = TRUE)

cat("\n=== Task 5: Dead chrX/chr22 ===\n")
d1 <- run_chrom(mk(DN, "chrX",  5:36), NULL, "Dead chrX Normal")
d2 <- run_chrom(mk(DN, "chr22", 5:36), NULL, "Dead chr22 Normal")
d3 <- run_chrom(mk(DT, "chrX",  5:36), NULL, "Dead chrX Tumor")
d4 <- run_chrom(mk(DT, "chr22", 5:36), NULL, "Dead chr22 Tumor")
cat("Mean per-CpG outlier count -> chrX N:", mean(d1$per.cpg), " chr22 N:", mean(d2$per.cpg),
    " chrX T:", mean(d3$per.cpg), " chr22 T:", mean(d4$per.cpg), "\n")

# ---- Task 4: 5 specific CpGs, real OutlierMeth flags ----
cat("\n=== Task 4: 5-CpG printout (real OutlierMeth) ===\n")
five_cpgs <- c("cg27487046","cg06838427","cg01150641","cg23625715","cg06973463")

tumor.beta.full <- as.matrix(AT[, 5:57])
rownames(tumor.beta.full) <- AT$Composite.Element.REF
sub5 <- tumor.beta.full[five_cpgs, ]

ref5  <- referenceMeth(sub5)
flag5 <- flagMeth(sub5[rownames(ref5), ], reference = ref5, p = 0.01)

cat("Per-CpG outlier sample counts (Tumor):\n")
print(rowSums(abs(flag5), na.rm = TRUE))
cat("\nRaw flags (1=high outlier, -1=low outlier, 0=normal):\n")
print(flag5)

cat("\nDONE.\n")
