cat("Loading datasets...\n")

setwd("/home/s_s355/research3.data/TCGA/filter.BRCA.UCEC.Aug17.2022/sorted.TCGA.380355cg.files")

normal <- read.table(
"Sorted.BRCA.53Alive.Normal.380355cg.75col.May28.2026.txt",
header=TRUE)

tumor <- read.table(
"Sorted.BRCA.53Alive.Tumor.380355cg.75col.May28.2026.txt",
header=TRUE)

cat("Original rows\n")
print(nrow(normal))

###################################################
# Remove H and L
###################################################

normal.filtered <-
subset(normal,
!(methy.state %in% c("H","L")))

tumor.filtered <-
subset(tumor,
!(methy.state %in% c("H","L")))

cat("\nFiltered rows\n")
print(nrow(normal.filtered))
print(nrow(tumor.filtered))

###################################################
# Compare coef2
###################################################

cat("\nMean coef2 Before\n")
print(mean(normal$outliers.coef2))
print(mean(tumor$outliers.coef2))

cat("\nMean coef2 After\n")
print(mean(normal.filtered$outliers.coef2))
print(mean(tumor.filtered$outliers.coef2))

###################################################
# Compare coef3
###################################################

cat("\nMean coef3 Before\n")
print(mean(normal$outliers.coef3))
print(mean(tumor$outliers.coef3))

cat("\nMean coef3 After\n")
print(mean(normal.filtered$outliers.coef3))
print(mean(tumor.filtered$outliers.coef3))

###################################################
# States remaining
###################################################

cat("\nRemaining states\n")
print(table(normal.filtered$methy.state))
print(table(tumor.filtered$methy.state))
