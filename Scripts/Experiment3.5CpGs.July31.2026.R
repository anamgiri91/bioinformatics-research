setwd("/home/s_s355/research3.data/TCGA/filter.BRCA.UCEC.Aug17.2022/sorted.TCGA.380355cg.files")

BRCA.53AN <- read.table(
"Sorted.BRCA.53Alive.Normal.380355cg.75col.May28.2026.txt",
header=TRUE)

BRCA.53AT <- read.table(
"Sorted.BRCA.53Alive.Tumor.380355cg.75col.May28.2026.txt",
header=TRUE)

outliers <- subset(BRCA.53AT,
outliers.coef3>=8)

head(outliers,5)

outliers[1:5,c(1:57)]

BRCA.53AN[
BRCA.53AN$Composite.Element.REF %in%
outliers$Composite.Element.REF[1:5],
c(1:57)]
