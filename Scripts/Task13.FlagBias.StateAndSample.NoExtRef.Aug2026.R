## ============================================================================
## Task 13 -- Outlier flag (-1,0,1) bias analysis, NO external reference
##
## Follows up on Aug 6, 2026 notes:
##   1) L-state tends to show more "outliers" under boxplot/normality-based
##      methods (outlier2) -- check whether the same L/H bias shows up in the
##      OutlierMeth flag (-1/0/1) system.
##   2) OutlierMeth does not use local spatial correlation.
##   3) Top/bottom 1% flagging may not be "proper" across all states.
##   4) For L-state CpGs low methylation is EXPECTED (not a true outlier);
##      for H-state CpGs high methylation is EXPECTED. Flagging the "-1" at
##      L-state sites and "+1" at H-state sites may therefore be mislabeling
##      biology as an outlier.
##
## This script:
##   Part A - Flag vs methylation state
##            - contingency table of flag (-1/0/1) x state  (counts + row %)
##            - per-CpG "outlier rate" (fraction of samples flagged, any
##              direction) summarized by state using BOTH mean and median,
##              to see if mean/median diverge (a sign the state's
##              distribution is skewed / driven by a handful of CpGs)
##            - same summary split by direction (-1 only, +1 only) since the
##              L/H "expected extreme" issue is direction-specific
##   Part B - Flag vs sample (are samples 14-17 unusual?)
##            - per-sample total |flag| count (Normal, Tumor)
##            - z-score and rank of samples 14-17 against the other 49
##              samples
##            - explicit flag: TRUE if sample falls outside mean +/- 2 SD
##              (and separately, outside 1.5*IQR) of the *other* samples
##
## Data loading / column layout kept identical to Task10 / Task12.
## ============================================================================

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
cat("Finished loading.\n")

normal.beta <- as.matrix(BRCA.53AN[, 5:57])
tumor.beta  <- as.matrix(BRCA.53AT[, 5:57])
rownames(normal.beta) <- BRCA.53AN$Composite.Element.REF
rownames(tumor.beta)  <- BRCA.53AT$Composite.Element.REF

methy_lookup_normal <- setNames(BRCA.53AN$methy.state, BRCA.53AN$Composite.Element.REF)
methy_lookup_tumor  <- setNames(BRCA.53AT$methy.state, BRCA.53AT$Composite.Element.REF)

state_order <- c("L", "LM", "M", "HM", "H", "Rc")

cat("\nBuilding self-reference + flags (Normal) ...\n")
ref.normal  <- referenceMeth(normal.beta)
flag.normal <- flagMeth(normal.beta[rownames(ref.normal), ], reference = ref.normal, p = 0.01)

cat("Building self-reference + flags (Tumor) ...\n")
ref.tumor  <- referenceMeth(tumor.beta)
flag.tumor <- flagMeth(tumor.beta[rownames(ref.tumor), ], reference = ref.tumor, p = 0.01)

## ---------------------------------------------------------------------
## Part A: Flag vs methylation state
## ---------------------------------------------------------------------

# A1: contingency table of flag x state (counts + row %) -- same as Task10
build_state_flag_table <- function(flag_matrix, methy_state_lookup) {
  flag_df <- as.data.frame(flag_matrix)
  flag_df$cpg <- rownames(flag_matrix)

  flag_long <- flag_df %>%
    pivot_longer(-cpg, names_to = "sample", values_to = "flag")
  flag_long$methy.state <- methy_state_lookup[flag_long$cpg]

  tab <- table(flag_long$methy.state, flag_long$flag)
  ord <- state_order[state_order %in% rownames(tab)]
  tab <- tab[ord, , drop = FALSE]
  pct <- round(prop.table(tab, margin = 1) * 100, 2)

  list(counts = tab, percent = pct, flag_long = flag_long)
}

# A2: per-CpG outlier rate (fraction of samples flagged), by state,
#     summarized with mean AND median -- overall and by direction
state_bias_summary <- function(flag_matrix, methy_state_lookup) {
  n.samples <- ncol(flag_matrix)

  rate.any   <- rowMeans(abs(flag_matrix) == 1, na.rm = TRUE)   # any flag
  rate.neg1  <- rowMeans(flag_matrix == -1, na.rm = TRUE)       # "-1" only
  rate.pos1  <- rowMeans(flag_matrix ==  1, na.rm = TRUE)       # "+1" only

  df <- data.frame(
    cpg   = rownames(flag_matrix),
    state = methy_state_lookup[rownames(flag_matrix)],
    rate.any  = rate.any,
    rate.neg1 = rate.neg1,
    rate.pos1 = rate.pos1
  )

  summ <- df %>%
    group_by(state) %>%
    summarise(
      n.cpg          = n(),
      mean.any       = mean(rate.any,  na.rm = TRUE),
      median.any     = median(rate.any, na.rm = TRUE),
      mean.neg1      = mean(rate.neg1, na.rm = TRUE),
      median.neg1    = median(rate.neg1, na.rm = TRUE),
      mean.pos1      = mean(rate.pos1, na.rm = TRUE),
      median.pos1    = median(rate.pos1, na.rm = TRUE),
      .groups = "drop"
    )

  ord <- state_order[state_order %in% summ$state]
  summ <- summ[match(ord, summ$state), ]
  summ
}

cat("\n===== Part A: Flag x State contingency -- Normal =====\n")
tabA.normal <- build_state_flag_table(flag.normal, methy_lookup_normal)
print(tabA.normal$percent)

cat("\n===== Part A: Flag x State contingency -- Tumor =====\n")
tabA.tumor <- build_state_flag_table(flag.tumor, methy_lookup_tumor)
print(tabA.tumor$percent)

cat("\n===== Part A: mean vs median outlier RATE by state -- Normal =====\n")
bias.normal <- state_bias_summary(flag.normal, methy_lookup_normal)
print(bias.normal)

cat("\n===== Part A: mean vs median outlier RATE by state -- Tumor =====\n")
bias.tumor <- state_bias_summary(flag.tumor, methy_lookup_tumor)
print(bias.tumor)

write.csv(as.data.frame.matrix(tabA.normal$percent), "task13_stateFlagPct_normal_noExtRef.csv")
write.csv(as.data.frame.matrix(tabA.tumor$percent),  "task13_stateFlagPct_tumor_noExtRef.csv")
write.csv(bias.normal, "task13_stateBias_meanMedian_normal_noExtRef.csv", row.names = FALSE)
write.csv(bias.tumor,  "task13_stateBias_meanMedian_tumor_noExtRef.csv", row.names = FALSE)

## ---------------------------------------------------------------------
## Part B: Flag vs sample -- are samples 14-17 unusual?
## ---------------------------------------------------------------------

sample_outlier_check <- function(flag_matrix, target_ids, label) {
  per.sample.total <- colSums(abs(flag_matrix), na.rm = TRUE)
  per.sample.neg1  <- colSums(flag_matrix == -1, na.rm = TRUE)
  per.sample.pos1  <- colSums(flag_matrix ==  1, na.rm = TRUE)
  per.sample.zero  <- colSums(flag_matrix ==  0, na.rm = TRUE)

  targets <- intersect(target_ids, names(per.sample.total))
  others  <- setdiff(names(per.sample.total), targets)

  mu   <- mean(per.sample.total[others])
  sdv  <- sd(per.sample.total[others])
  q1   <- quantile(per.sample.total[others], 0.25)
  q3   <- quantile(per.sample.total[others], 0.75)
  iqr  <- q3 - q1
  lo.iqr <- q1 - 1.5 * iqr
  hi.iqr <- q3 + 1.5 * iqr

  out <- data.frame(
    sample     = names(per.sample.total),
    total.flag = as.numeric(per.sample.total),
    neg1.count = as.numeric(per.sample.neg1),
    zero.count = as.numeric(per.sample.zero),
    pos1.count = as.numeric(per.sample.pos1),
    z.vs.others = (as.numeric(per.sample.total) - mu) / sdv,
    rank.desc   = rank(-as.numeric(per.sample.total)),
    flag.2SD    = abs((as.numeric(per.sample.total) - mu) / sdv) > 2,
    flag.1.5IQR = as.numeric(per.sample.total) < lo.iqr | as.numeric(per.sample.total) > hi.iqr,
    is.target14_17 = names(per.sample.total) %in% targets
  )
  out <- out[order(out$rank.desc), ]

  cat("\n--", label, "-- reference stats from the other", length(others), "samples:\n")
  cat("   mean =", round(mu, 3), " SD =", round(sdv, 3),
      " 1.5xIQR bounds = [", round(lo.iqr, 3), ",", round(hi.iqr, 3), "]\n")
  cat("   Samples 14-17 rows:\n")
  print(out[out$is.target14_17, ])

  out
}

target.normal <- paste0("N", 14:17)
target.tumor  <- paste0("T", 14:17)

cat("\n===== Part B: sample-level flag summary -- Normal =====\n")
sampB.normal <- sample_outlier_check(flag.normal, target.normal, "Normal")

cat("\n===== Part B: sample-level flag summary -- Tumor =====\n")
sampB.tumor <- sample_outlier_check(flag.tumor, target.tumor, "Tumor")

write.csv(sampB.normal, "task13_sampleFlagSummary_normal_noExtRef.csv", row.names = FALSE)
write.csv(sampB.tumor,  "task13_sampleFlagSummary_tumor_noExtRef.csv", row.names = FALSE)

cat("\nDONE. Task13 (no external reference) CSVs written to working directory.\n")