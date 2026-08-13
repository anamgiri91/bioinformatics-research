cat("Loading datasets...\n")

setwd("/home/s_s355/research3.data/TCGA/filter.BRCA.UCEC.Aug17.2022/sorted.TCGA.380355cg.files")

alive.normal <- read.table(
"Sorted.BRCA.53Alive.Normal.380355cg.75col.May28.2026.txt",
header=TRUE)

alive.tumor <- read.table(
"Sorted.BRCA.53Alive.Tumor.380355cg.75col.May28.2026.txt",
header=TRUE)

dead.normal <- read.table(
"Sorted.BRCA.32Dead.Normal.380355cg.54col.May28.2026.txt",
header=TRUE)

dead.tumor <- read.table(
"Sorted.BRCA.32Dead.Tumor.380355cg.54col.May28.2026.txt",
header=TRUE)

cat("Finished loading.\n\n")

##########################################################
# Combine beta values
##########################################################

alive.normal.beta <- alive.normal[,5:57]
alive.tumor.beta  <- alive.tumor[,5:57]

dead.normal.beta <- dead.normal[,5:36]
dead.tumor.beta  <- dead.tumor[,5:36]

combined.normal <- cbind(alive.normal.beta,
                         dead.normal.beta)

combined.tumor <- cbind(alive.tumor.beta,
                        dead.tumor.beta)

cat("Dimensions\n")
print(dim(combined.normal))
print(dim(combined.tumor))

##########################################################
# Numeric summary
##########################################################

summary.table <- data.frame(

Dataset=c("Combined Normal",
          "Combined Tumor"),

Mean=c(mean(as.matrix(combined.normal)),
       mean(as.matrix(combined.tumor))),

Median=c(median(as.matrix(combined.normal)),
         median(as.matrix(combined.tumor))),

SD=c(sd(as.matrix(combined.normal)),
     sd(as.matrix(combined.tumor)))
)

cat("\nNumeric Summary\n")
print(summary.table)

##########################################################
# Existing Outlier coefficients
##########################################################

combined.coef2.normal <-
c(alive.normal$outliers.coef2,
dead.normal$outliers.coef2)

combined.coef2.tumor <-
c(alive.tumor$outliers.coef2,
dead.tumor$outliers.coef2)

combined.coef3.normal <-
c(alive.normal$outliers.coef3,
dead.normal$outliers.coef3)

combined.coef3.tumor <-
c(alive.tumor$outliers.coef3,
dead.tumor$outliers.coef3)

cat("\nMean coef2\n")
print(mean(combined.coef2.normal))
print(mean(combined.coef2.tumor))

cat("\nMean coef3\n")
print(mean(combined.coef3.normal))
print(mean(combined.coef3.tumor))
