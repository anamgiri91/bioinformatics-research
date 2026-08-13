cat("Loading BRCA 32 Dead datasets...\n")

setwd("/home/s_s355/research3.data/TCGA/filter.BRCA.UCEC.Aug17.2022/sorted.TCGA.380355cg.files")

BRCA.32DN <- read.table(
"Sorted.BRCA.32Dead.Normal.380355cg.54col.May28.2026.txt",
header=TRUE)

BRCA.32DT <- read.table(
"Sorted.BRCA.32Dead.Tumor.380355cg.54col.May28.2026.txt",
header=TRUE)

cat("Finished loading.\n")

dim(BRCA.32DN)
dim(BRCA.32DT)
cat("\n===== Numeric Summary =====\n")

normal.beta <- as.matrix(BRCA.32DN[,5:36])
tumor.beta  <- as.matrix(BRCA.32DT[,5:36])

numeric.summary <- data.frame(
  Dataset=c("Dead Normal","Dead Tumor"),
  Mean=c(mean(normal.beta), mean(tumor.beta)),
  Median=c(median(normal.beta), median(tumor.beta)),
  SD=c(sd(normal.beta), sd(tumor.beta)),
  Min=c(min(normal.beta), min(tumor.beta)),
  Max=c(max(normal.beta), max(tumor.beta))
)

print(numeric.summary)
cat("\n===== coef2 =====\n")

print(summary(BRCA.32DN$outliers.coef2))
print(summary(BRCA.32DT$outliers.coef2))

cat("\n===== coef3 =====\n")

print(summary(BRCA.32DN$outliers.coef3))
print(summary(BRCA.32DT$outliers.coef3))

cat("\n===== coef2 counts =====\n")

print(table(BRCA.32DN$outliers.coef2))
print(table(BRCA.32DT$outliers.coef2))

cat("\n===== coef3 counts =====\n")

print(table(BRCA.32DN$outliers.coef3))
print(table(BRCA.32DT$outliers.coef3))
cat("\n===== chrX vs chr22 =====\n")

chrX.normal <- subset(BRCA.32DN, Chromosome=="chrX")
chr22.normal <- subset(BRCA.32DN, Chromosome=="chr22")

chrX.tumor <- subset(BRCA.32DT, Chromosome=="chrX")
chr22.tumor <- subset(BRCA.32DT, Chromosome=="chr22")

cat("Rows\n")
cat("chrX Normal:", nrow(chrX.normal), "\n")
cat("chr22 Normal:", nrow(chr22.normal), "\n")
cat("chrX Tumor :", nrow(chrX.tumor), "\n")
cat("chr22 Tumor:", nrow(chr22.tumor), "\n")

cat("\nMean coef2\n")
print(c(
  chrX.normal=mean(chrX.normal$outliers.coef2),
  chr22.normal=mean(chr22.normal$outliers.coef2),
  chrX.tumor=mean(chrX.tumor$outliers.coef2),
  chr22.tumor=mean(chr22.tumor$outliers.coef2)
))

cat("\nMean coef3\n")
print(c(
  chrX.normal=mean(chrX.normal$outliers.coef3),
  chr22.normal=mean(chr22.normal$outliers.coef3),
  chrX.tumor=mean(chrX.tumor$outliers.coef3),
  chr22.tumor=mean(chr22.tumor$outliers.coef3)
))

cat("\nMethylation States\n")
print(table(chrX.normal$methy.state))
print(table(chr22.normal$methy.state))
print(table(chrX.tumor$methy.state))
print(table(chr22.tumor$methy.state))
