# ================================================================
# Task 52 - Section 6 of the frozen plan: do close CpGs share reads?
# Date: Sep 29, 2026
#
# Runs section 6 of Results/Task47_FrozenProtocol.md (commit a832294)
# on the GTEx colon RRBS cohort (GSE233417, 29 donors), using the
# read-level linkage built by Task 50 and the held-out outcomes saved
# by Task 51.
#
# LINKAGE (frozen)
#   For adjacent CpGs (genome CpG index i and i + 1) within 200 bp:
#   among reads that call both sites, the share with the same state at
#   both, minus p_i p_j + (1 - p_i)(1 - p_j) from the same reads. It is
#   computed within each sample (at least 10 such reads) and then
#   averaged over samples, so differences between people cannot
#   create it.
#
# HYPOTHESES (frozen)
#   H8  mean linkage is above 0, and it falls with distance
#   H9  among low-spread pairs of the section 5 cohort (one site has no
#       sample more than 0.1 from its median; the pairs the spread
#       guard penalised), pairs with above-median linkage have higher
#       held-out prediction gain than pairs below the median, in both
#       halves.
#   The verdict uses the point estimates, as the plan states them.
#   95% intervals from a bootstrap of 1-Mb blocks (1,000 resamples,
#   seed 20260929) are added for information.
#
# OUTPUTS (Results/)
#   Task52_H8_ByDistance.csv, Task52_Hypotheses.csv,
#   Fig38_SharedReads.png
# ================================================================

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(jsonlite) })

pick_dir <- function(...) { for (d in c(...)) if (dir.exists(d)) return(d)
  stop("no candidate directory exists") }
repo_dir <- pick_dir("/mmfs1/home/wln26/Experiments.Outlier.July31.2026",
                     path.expand("~/Desktop/bioinformatics-research"))
res <- file.path(repo_dir, "Results"); dat <- file.path(repo_dir, "Data", "gtex_colon_rrbs")
r4 <- function(x) round(x, 4)
proto <- fromJSON(file.path(res, "Task47_FrozenProtocol.json"))
for (f in names(proto$files_sha256))
  stopifnot(digest::digest(file = file.path(repo_dir, f), algo = "sha256") == proto$files_sha256[[f]])
cat("freeze verified\n")
B <- 1000L; SEED <- 20260929L
source(file.path(repo_dir,"Scripts","GTExColon.BlockBootstrap.R"))
check_block_bootstrap()

cc0 <- readRDS(file.path(dat, "cache_cohort.rds"))
mp <- cc0$map; setkey(mp, chrn, idx)
lk <- copy(cc0$linkage)
rm(cc0); invisible(gc())
lk[, pos := mp[.(lk$chrn, lk$idx), pos]]
lk[, pos2 := mp[.(lk$chrn, lk$idx + 1L), pos]]
lk <- lk[!is.na(pos) & !is.na(pos2)][, gap := pos2 - pos][gap > 0 & gap <= 200]
lk[, half := fifelse(pos < median(pos), "first", "second"), by = chrn]
lk[, block := paste(chrn, pos %/% 1e6L)]
cat("adjacent pairs within 200 bp with read-level linkage:", nrow(lk), "\n")

block_boot <- function(d, kind) block_boot_exact(d,kind,B=B,seed=SEED)

# ----------------------------------------------------------------
# H8
# ----------------------------------------------------------------
SUB <- list(all = rep(TRUE, nrow(lk)), first = lk$half == "first", second = lk$half == "second")
h8 <- rbindlist(lapply(names(SUB), function(sn) {
  d <- lk[SUB[[sn]]]
  m <- mean(d$linkage); rho <- cor(d$linkage, d$gap, method = "spearman")
  ci_m <- block_boot(d, "mean"); ci_r <- block_boot(d, "spearman")
  data.table(hypothesis = "H8", subset = sn, pairs = nrow(d), mean_linkage = r4(m), lo95_mean = r4(ci_m[1]), hi95_mean = r4(ci_m[2]),
             share_positive = r4(mean(d$linkage > 0)), spearman_with_distance = r4(rho),
             lo95_rho = r4(ci_r[1]), hi95_rho = r4(ci_r[2]), holds = m > 0 & rho < 0)
}))
bins <- lk[, .(pairs = .N, mean_linkage = r4(mean(linkage)), median_gap = as.numeric(median(gap))),
           by = .(bin = cut(gap, c(0, 10, 20, 40, 60, 100, 150, 200)))][order(bin)]
fwrite(bins, file.path(res, "Task52_H8_ByDistance.csv"))
cat("\nH8: read-level linkage of adjacent CpGs within 200 bp\n"); print(h8); print(bins)

# ----------------------------------------------------------------
# H9
# ----------------------------------------------------------------
po <- readRDS(file.path(dat, "pair_outcomes.rds"))
h9d <- merge(po[low_spread == TRUE & gap <= 200], lk[, .(chrn, idx, linkage, block)], by = c("chrn", "idx"))
cat("\nlow-spread cohort pairs within 200 bp with linkage:", nrow(h9d), "of", sum(po$low_spread & po$gap <= 200), "\n")
diff_above <- function(x) { md <- median(x$linkage)
  mean(x$heldout_gain[x$linkage > md]) - mean(x$heldout_gain[x$linkage <= md]) }
SUB9 <- list(all = rep(TRUE, nrow(h9d)), first = h9d$half == "first", second = h9d$half == "second")
h9 <- rbindlist(lapply(names(SUB9), function(sn) {
  d <- h9d[SUB9[[sn]]]; md <- median(d$linkage)
  dd <- diff_above(d); ci <- block_boot(d, "gain_difference")
  data.table(hypothesis = "H9", subset = sn, pairs = nrow(d), median_linkage = r4(md),
             gain_above = round(mean(d$heldout_gain[d$linkage > md]), 5),
             gain_below = round(mean(d$heldout_gain[d$linkage <= md]), 5),
             difference = round(dd, 5), lo95 = round(ci[1], 5), hi95 = round(ci[2], 5), holds = dd > 0)
}))
cat("\nH9: held-out prediction gain of low-spread pairs, above vs below median linkage\n"); print(h9)
fwrite(h8,file.path(res,"Task52_H8_Detail.csv"))
fwrite(h9,file.path(res,"Task52_H9_Detail.csv"))
# exploratory: the same across all pairs, not only low-spread
ex <- merge(po[gap <= 200], lk[, .(chrn, idx, linkage)], by = c("chrn", "idx"))
cat("\nExploratory, all cohort pairs within 200 bp: Spearman(linkage, held-out gain) =",
    r4(cor(ex$linkage, ex$heldout_gain, method = "spearman")), "over", nrow(ex), "pairs\n")

hyps <- rbind(h8[, .(hypothesis, subset, pairs, estimate = mean_linkage, lo95 = lo95_mean, hi95 = hi95_mean,
                     second_estimate = spearman_with_distance, holds)],
              h9[, .(hypothesis, subset, pairs, estimate = difference, lo95, hi95, second_estimate = NA_real_, holds)])
fwrite(hyps, file.path(res, "Task52_Hypotheses.csv"))
cat("\nVerdicts:\n"); print(hyps[, .(holds_everywhere = all(holds)), by = hypothesis])

# ----------------------------------------------------------------
# figure
# ----------------------------------------------------------------
SURF <- "#fcfcfb"; INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"; GRID <- "#e1e0d9"; BLUE <- "#2a78d6"
bci <- rbindlist(lapply(levels(lk$bin <- cut(lk$gap, c(0, 10, 20, 40, 60, 100, 150, 200))), function(bn) {
  d <- lk[bin == bn]; ci <- block_boot(d, "mean")
  data.table(bin = bn, mean = mean(d$linkage), lo = ci[1], hi = ci[2], pairs = nrow(d))
}))
bci[, bin := factor(bin, levels = bin)]
pA <- ggplot(bci, aes(bin, mean)) +
  geom_hline(yintercept = 0, colour = MUTED, linewidth = 0.4) +
  geom_linerange(aes(ymin = lo, ymax = hi), colour = "#c3c2b7", linewidth = 1) +
  geom_point(colour = BLUE, size = 3) +
  labs(title = "A. Read-level linkage of adjacent CpGs, by distance",
       x = "Distance between the two CpGs (bp)", y = "Linkage (same-state share above chance)") +
  theme_minimal(base_size = 10.5)
h9p <- h9[, .(group = c("Below-median linkage", "Above-median linkage"), gain = c(gain_below, gain_above)), by = subset]
h9p[, subset := factor(subset, levels = c("all", "first", "second"), labels = c("All pairs", "First halves", "Second halves"))]
pB <- ggplot(h9p, aes(gain, subset, fill = group)) +
  geom_point(shape = 21, colour = SURF, stroke = 0.6, size = 3.4) +
  scale_fill_manual(values = c(`Below-median linkage` = MUTED, `Above-median linkage` = BLUE), name = NULL) +
  labs(title = "B. Low-spread pairs: held-out prediction gain", x = "Held-out prediction gain (higher is better)", y = NULL) +
  theme_minimal(base_size = 10.5) + theme(legend.position = "top", legend.justification = "left")
th <- theme(plot.background = element_rect(fill = SURF, colour = NA), panel.grid.minor = element_blank(),
            panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
            plot.title = element_text(colour = INK, face = "bold", size = 11.5), axis.text = element_text(colour = INK2),
            axis.title = element_text(colour = INK2))
png(file.path(res, "Fig38_SharedReads.png"), width = 11, height = 4.6, units = "in", res = 200, bg = SURF)
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 2, heights = grid::unit(c(0.12, 0.88), "npc"))))
grid::grid.text("Do close CpGs share reads? GTEx colon RRBS, 29 donors (frozen section 6)",
                x = 0.01, hjust = 0, gp = grid::gpar(fontface = "bold", fontsize = 13, col = INK),
                vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1:2))
print(pA + th, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
print(pB + th, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 2))
invisible(dev.off())
cat("\ndone\n")
