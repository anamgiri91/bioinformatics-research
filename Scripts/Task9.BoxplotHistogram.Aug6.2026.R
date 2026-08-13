library(dplyr)
library(tidyr)
library(ggplot2)
 
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
 
# ---- helper: reshape one dataset (optionally filtered) into long form ----
make_long <- function(df, states_to_keep = NULL) {
  d <- df %>% select(methy.state, outliers.coef2, outliers.coef3)
  if (!is.null(states_to_keep)) {
    d <- d %>% filter(methy.state %in% states_to_keep)
  }
  d %>%
    pivot_longer(cols = c(outliers.coef2, outliers.coef3),
                 names_to = "metric", values_to = "count") %>%
    mutate(metric = recode(metric,
                            outliers.coef2 = "outlier2",
                            outliers.coef3 = "outlier3"))
}
 
# ---- helper: side-by-side boxplot (outlier2 vs outlier3 per state) ----
make_boxplot <- function(d_long, title_txt) {
  ggplot(d_long, aes(x = methy.state, y = count, fill = metric)) +
    geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.5) +
    labs(title = title_txt, x = "Methylation state",
         y = "Outlier count", fill = "Metric") +
    theme_minimal()
}
 
# ---- helper: histograms, one panel per state, SAME xlim across all panels ----
# global_xlim is passed in explicitly so before/after plots are comparable too,
# not just panels within a single plot.
make_histogram <- function(d_long, title_txt, global_xlim) {
  ggplot(d_long, aes(x = count, fill = metric)) +
    geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
    facet_wrap(~ methy.state, scales = "free_y") +
    coord_cartesian(xlim = global_xlim) +
    labs(title = title_txt, x = "Outlier count",
         y = "Frequency", fill = "Metric") +
    theme_minimal()
}
 
run_task9 <- function(df, dataset_label) {
 
  states.before <- unique(df$methy.state)
  states.after  <- setdiff(states.before, c("H", "L"))
 
  long.before <- make_long(df)
  long.after  <- make_long(df, states_to_keep = states.after)
 
  # one global x-range so before/after histograms are visually comparable
  global_xlim <- range(long.before$count, na.rm = TRUE)
 
  box.before <- make_boxplot(long.before,
                 paste0(dataset_label, " — BEFORE removing H/L: outlier2 vs outlier3 by state"))
  box.after  <- make_boxplot(long.after,
                 paste0(dataset_label, " — AFTER removing H/L: outlier2 vs outlier3 by state"))
 
  hist.before <- make_histogram(long.before,
                 paste0(dataset_label, " — BEFORE removing H/L: distribution by state"),
                 global_xlim)
  hist.after  <- make_histogram(long.after,
                 paste0(dataset_label, " — AFTER removing H/L: distribution by state"),
                 global_xlim)
 
  list(box.before = box.before, box.after = box.after,
       hist.before = hist.before, hist.after = hist.after)
}
 
cat("\nRunning Task 9 for Normal...\n")
plots.normal <- run_task9(BRCA.53AN, "Normal")
 
cat("Running Task 9 for Tumor...\n")
plots.tumor <- run_task9(BRCA.53AT, "Tumor")
 
ggsave("/mmfs1/home/wln26/Experiments.Outlier.July31.2026/Results/task9_normal_box_before.png",  plots.normal$box.before,  width = 8, height = 5)
ggsave("/mmfs1/home/wln26/Experiments.Outlier.July31.2026/Results/task9_normal_box_after.png",   plots.normal$box.after,   width = 8, height = 5)
ggsave("/mmfs1/home/wln26/Experiments.Outlier.July31.2026/Results/task9_normal_hist_before.png", plots.normal$hist.before, width = 9, height = 6)
ggsave("/mmfs1/home/wln26/Experiments.Outlier.July31.2026/Results/task9_normal_hist_after.png",  plots.normal$hist.after,  width = 9, height = 6)
 
ggsave("/mmfs1/home/wln26/Experiments.Outlier.July31.2026/Results/task9_tumor_box_before.png",   plots.tumor$box.before,   width = 8, height = 5)
ggsave("/mmfs1/home/wln26/Experiments.Outlier.July31.2026/Results/task9_tumor_box_after.png",    plots.tumor$box.after,    width = 8, height = 5)
ggsave("/mmfs1/home/wln26/Experiments.Outlier.July31.2026/Results/task9_tumor_hist_before.png",  plots.tumor$hist.before,  width = 9, height = 6)
ggsave("/mmfs1/home/wln26/Experiments.Outlier.July31.2026/Results/task9_tumor_hist_after.png",   plots.tumor$hist.after,   width = 9, height = 6)
 
cat("\nDONE. 8 PNGs written to working directory.\n")
