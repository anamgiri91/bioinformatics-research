base_dir <- "/home/s_s355/research3.data/TCGA/filter.BRCA.UCEC.Aug17.2022/sorted.TCGA.380355cg.files"

BRCA.53AN <- read.table(
    file.path(base_dir, "Sorted.BRCA.53Alive.Normal.380355cg.75col.May28.2026.txt"),
    header = TRUE
)

BRCA.53AT <- read.table(
    file.path(base_dir, "Sorted.BRCA.53Alive.Tumor.380355cg.75col.May28.2026.txt"),
    header = TRUE
)
