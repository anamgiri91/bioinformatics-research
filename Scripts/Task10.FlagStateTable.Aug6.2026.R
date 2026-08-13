
library(OutlierMeth)
library(dplyr)
library(tidyr)
 
cat("Loading BRCA 53 Alive datasets...\n")
BRCA.53AN <- read.table(
  "Sorted.BRCA.53Alive.Normal.380355cg.75col.May28.2026.txt",
  header = TRUE
)
BRCA.53AT <- read.table(
  "Sorted.BRCA.53Alive.Tumor.380355cg.75col.May28.2026.txt",
  header = TRUE
)
 
normal.beta <- as.matrix(BRCA.53AN[, 5:57])
tumor.beta  <- as.matrix(BRCA.53AT[, 5:57])
rownames(normal.beta) <- BRCA.53AN$Composite.Element.REF
rownames(tumor.beta)  <- BRCA.53AT$Composite.Element.REF
 
cat("Building reference + flags (Normal) ...\n")
ref.normal  <- referenceMeth(normal.beta)
flag.normal <- flagMeth(normal.beta[rownames(ref.normal), ], reference = ref.normal, p = 0.01)
 
cat("Building reference + flags (Tumor) ...\n")
ref.tumor  <- referenceMeth(tumor.beta)
flag.tumor <- flagMeth(tumor.beta[rownames(ref.tumor), ], reference = ref.tumor, p = 0.01)
 
# ---- helper: cross-tabulate flag value against methy.state ----
build_state_flag_table <- function(flag_matrix, methy_state_lookup) {
  flag_df <- as.data.frame(flag_matrix)
  flag_df$cpg <- rownames(flag_matrix)
 
  flag_long <- flag_df %>%
    pivot_longer(-cpg, names_to = "sample", values_to = "flag")
 
  flag_long$methy.state <- methy_state_lookup[flag_long$cpg]
 
  tab <- table(flag_long$methy.state, flag_long$flag)
 
  # enforce a consistent row order; drop any state not present in the data
  state_order <- c("L", "LM", "M", "HM", "H", "Rc")
  state_order <- state_order[state_order %in% rownames(tab)]
  tab <- tab[state_order, , drop = FALSE]
 
  pct <- round(prop.table(tab, margin = 1) * 100, 2)  # row % (within each state)
 
  combined <- matrix(
    paste0(tab, " (", pct, "%)"),
    nrow = nrow(tab), dimnames = dimnames(tab)
  )
 
  list(counts = tab, percent = pct, combined = combined)
}
 
methy_lookup_normal <- setNames(BRCA.53AN$methy.state, BRCA.53AN$Composite.Element.REF)
methy_lookup_tumor  <- setNames(BRCA.53AT$methy.state, BRCA.53AT$Composite.Element.REF)
 
cat("\n===== Flag x Methylation State — Normal =====\n")
result.normal <- build_state_flag_table(flag.normal, methy_lookup_normal)
print(result.normal$combined)
 
cat("\n===== Flag x Methylation State — Tumor =====\n")
result.tumor <- build_state_flag_table(flag.tumor, methy_lookup_tumor)
print(result.tumor$combined)
 
write.csv(result.normal$combined, "/mmfs1/home/wln26/Experiments.Outlier.July31.2026/Results/task10_flag_by_state_normal.csv")
write.csv(result.tumor$combined,  "/mmfs1/home/wln26/Experiments.Outlier.July31.2026/Results/task10_flag_by_state_tumor.csv")
 
cat("\nDONE. CSVs written to working directory.\n")
