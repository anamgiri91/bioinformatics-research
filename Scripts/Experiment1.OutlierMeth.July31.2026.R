setwd("/home/s_s355/research3.data/TCGA/filter.BRCA.UCEC.Aug17.2022/sorted.TCGA.380355cg.files")

####################################################
# Experiment 1
# BRCA 53 Alive
# Normal vs Tumor
####################################################

cat("Loading BRCA datasets...\n")

BRCA.53AN <- read.table(
    "Sorted.BRCA.53Alive.Normal.380355cg.75col.May28.2026.txt",
    header = TRUE
)

BRCA.53AT <- read.table(
    "Sorted.BRCA.53Alive.Tumor.380355cg.75col.May28.2026.txt",
    header = TRUE
)

cat("Finished loading.\n")

dim(BRCA.53AN)
dim(BRCA.53AT)
