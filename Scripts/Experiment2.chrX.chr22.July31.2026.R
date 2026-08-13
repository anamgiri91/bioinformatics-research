####################################################
# Experiment 2
# chrX vs chr22
####################################################

setwd("/home/s_s355/research3.data/TCGA/filter.BRCA.UCEC.Aug17.2022/sorted.TCGA.380355cg.files")

BRCA.53AN <- read.table(
"Sorted.BRCA.53Alive.Normal.380355cg.75col.May28.2026.txt",
header=TRUE)

BRCA.53AT <- read.table(
"Sorted.BRCA.53Alive.Tumor.380355cg.75col.May28.2026.txt",
header=TRUE)

chrX.normal <- subset(BRCA.53AN,
Chromosome=="chrX")

chr22.normal <- subset(BRCA.53AN,
Chromosome=="chr22")

chrX.tumor <- subset(BRCA.53AT,
Chromosome=="chrX")

chr22.tumor <- subset(BRCA.53AT,
Chromosome=="chr22")

cat("\nRows\n")

cat("chrX Normal:",nrow(chrX.normal),"\n")
cat("chr22 Normal:",nrow(chr22.normal),"\n")
cat("chrX Tumor:",nrow(chrX.tumor),"\n")
cat("chr22 Tumor:",nrow(chr22.tumor),"\n")

cat("\nSummary coef2\n")

summary(chrX.normal$outliers.coef2)
summary(chr22.normal$outliers.coef2)
summary(chrX.tumor$outliers.coef2)
summary(chr22.tumor$outliers.coef2)

cat("\nSummary coef3\n")

summary(chrX.normal$outliers.coef3)
summary(chr22.normal$outliers.coef3)
summary(chrX.tumor$outliers.coef3)
summary(chr22.tumor$outliers.coef3)

cat("\nMethy.state\n")

table(chrX.normal$methy.state)
table(chr22.normal$methy.state)
table(chrX.tumor$methy.state)
table(chr22.tumor$methy.state)
