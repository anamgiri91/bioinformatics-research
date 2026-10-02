# ================================================================
# Task 39 - RRBS 18-sample lung data (chr22): site summaries,
#           methylation states, and similarity of neighbouring CpGs
# Date: Sep 28, 2026  (revised the same day after an audit of the
#                      first run; see REVISIONS below)
#
# THE REQUEST (supervisor, items 4-6)
#
#   4. (1) 6-number summary + outlier.IQR2 + outlier.IQR3, NA.count;
#          how many CG sites have 0 NA.
#      (2) methylation states for those sites; check the patterns of
#          CG sites close to each other.
#   5. Discretise beta into 10 levels, [0,0.1) -> 1 ... [0.9,1] -> 10,
#      and score the similarity of two consecutive CG sites; compare the
#      distances used in genetics / epigenetics, including Chatterjee's
#      xi (JASA 2021). Run on this RRBS data (confirmed 2026-09-28).
#   6. States as defined in the TSG paper: L, H, M, LM, HM, R.
#
# INPUT
#
#   Data/NN.hg38.18P.forw.chr22.w.header.txt   (gitignored)
#     chr, pos, P1NF..P18NF. One row per forward-strand CpG on hg38
#     chr22, pos = the C, 1-based (section 2 checks this against the
#     hg38 windows in Data/hg38_seq). Values are methylation ratios
#     (methylated / total reads). Read depth is not in the file.
#   TCGA-BRCA 53 Alive Normal / Tumour source matrices - used only in
#     section 1, to check the definitions below against the lab's own
#     columns. Skipped if absent.
#
# DEFINITIONS
#
#   Six-number summary: Min, Q1, Median, Mean, Q3, Max over the
#   non-missing samples; quantiles by R's type 7. outliers.coef2 / coef3
#   count samples below Q1 - k*IQR or above Q3 + k*IQR (k = 2, 3).
#
#   States - the TSG-paper rule (Sun, Pritchard, McFall & Tian 2025,
#   Epigenetics Insights 18:e004, "as defined by Tian et al. [13]" =
#   Tian, Bertelsmann, Yu & Sun 2016, Cancer Informatics 15:1-9), as
#   coded in meth.class() of
#   Oct11.methylation.status.identification.no.bimod.R:
#     L   all m < 0.2          H   all m > 0.8
#     M   no m < 0.2 and no m > 0.8
#     HM  all m > 0.4          LM  all m < 0.6
#     R   none of the above
#   HM and LM cannot both hold once L, H and M have failed (HM rules out
#   m < 0.2, so a non-M site must then have some m > 0.8, which rules
#   out LM), so the order of the last two tests does not matter.
#   Tian et al. also had a 'bimodal' state (over 40% of samples below
#   0.2 AND over 40% above 0.8) that the TSG paper folded into R; it is
#   counted here to show whether dropping it matters.
#
#   Discrete level: floor(10 m) + 1, with m = 1 put in level 10. The
#   levels line up with the state cut-offs except at exactly 0.4 and
#   0.8: 0.8 is level 9 but not "> 0.8", 0.4 is level 5 but not "> 0.4".
#   Ratios such as 2/5 and 4/5 land there often, so they are counted.
#
# SIMILARITY MEASURES, per pair of consecutive complete CpGs (x, y are
# the 18-sample vectors)
#
#   agree        share of samples with identical level (1 - Hamming/18);
#                genetics: identity-by-state between genotype vectors
#   agree1       share of samples within one level
#   L1           mean |level_x - level_y|        (Manhattan)
#   L2beta       root mean squared difference of the RAW betas
#   kappa_w      quadratic-weighted Cohen's kappa on the levels; it is
#                algebraically Lin's concordance correlation,
#                2 cov / (var_x + var_y + (mean_x - mean_y)^2)
#   pearson_beta Pearson r on the RAW betas, not on the levels
#   spearman     Spearman rho on the levels. This is the lab's own
#                co-methylation measure: the TSG paper calls |rho| >= 0.8
#                "highly correlated"
#   xi           Chatterjee's xi_n on the levels, made symmetric as
#                max(xi(x,y), xi(y,x)). Ties are broken at random, as the
#                paper prescribes; with 18 samples on ten levels ties are
#                everywhere and one draw moves xi a lot, so xi is the mean
#                over XI_DRAWS draws and the draw-to-draw SD is kept.
#                With n = 18 its maximum is 1 - 3/(n+1) = 0.842, not 1,
#                so it is not on the scale of the correlations.
#
#   When a measure is undefined:
#     spearman, xi   either level vector constant
#     pearson_beta   either beta vector constant
#     kappa_w        0 when exactly one vector is constant; undefined
#                    only when both are constant and equal
#   Only L, H and M sites can be constant (LM, HM and R need values on
#   both sides of a cut-off), so the correlation-type measures break at
#   L and H sites; the tables report how often.
#
# LITERATURE (checked 2026-09-28)
#
#   Affinito et al. 2020, Genomics 112:144-150: co-methylation of nearby
#     CpGs falls as the distance between them grows.
#   Guo et al. 2017, Nat Genet 49:635-642: methylation haplotype blocks,
#     adjacent CpGs with LD r^2 >= 0.5. There r^2 is taken across the
#     sequencing reads of one sample; it needs read-level data this file
#     does not have, and it is not a similarity across samples.
#
# REVISIONS (2026-09-28, after auditing the first run)
#
#   - Dropped the AUC of each measure for "same state". The target does
#     not fit: same-state R-R neighbours are dissimilar by construction
#     and L-L alike, and each measure had been scored on its own subset.
#   - Chance agreement of states is now also matched to each distance
#     bin's state mix; the single chromosome-wide figure understated it.
#   - Overall IQR-2/3 outlier counts are reported, next to the rate that
#     continuous data would give with 18 samples.
#   - pearson renamed pearson_beta; xi averaged over tie-breaks; kappa_w
#     computed in centred form, so a constant vector gives exactly 0.
#   - Added: section 1's TCGA check, the hg38 position check, the
#     bimodal count, value granularity, undefined rates by state pair,
#     distance decay within one state, kinds of state change, and the
#     email's own cg1 / cg2 example.
# ================================================================

suppressPackageStartupMessages({ library(data.table); library(matrixStats) })

pick_dir <- function(...) { for (d in c(...)) if (dir.exists(d)) return(d)
  stop("no candidate directory exists") }
repo_dir <- pick_dir("/mmfs1/home/wln26/Experiments.Outlier.July31.2026",
                     path.expand("~/Desktop/bioinformatics-research"))
out_dir  <- file.path(repo_dir, "Results")
in_file  <- file.path(repo_dir, "Data", "NN.hg38.18P.forw.chr22.w.header.txt")
seq_dir  <- file.path(repo_dir, "Data", "hg38_seq")
tcga_dir <- Filter(dir.exists, c(
  "/home/s_s355/research3.data/TCGA/filter.BRCA.UCEC.Aug17.2022/sorted.TCGA.380355cg.files",
  path.expand("~/Downloads/tcga_brca_380355cg")))
r4 <- function(x) round(x, 4)
set.seed(20260928)

XI_DRAWS <- 20     # random tie-breaks averaged into each xi
N_SIM    <- 1e5    # simulated continuous sites for the IQR-outlier yardstick
states   <- c("L", "LM", "M", "HM", "H", "R")

# ----------------------------------------------------------------
# 0. definitions shared by the RRBS data, the TCGA check and the
#    simulation, so all three go through the same code
# ----------------------------------------------------------------
six_num <- function(m) {
  q <- rowQuantiles(m, probs = c(0, .25, .5, .75, 1), na.rm = TRUE, type = 7L)
  iqr <- q[, 4] - q[, 2]
  out_k <- function(k) rowSums(m < q[, 2] - k * iqr | m > q[, 4] + k * iqr, na.rm = TRUE)
  data.table(Min = q[, 1], Q1 = q[, 2], Median = q[, 3], Mean = rowMeans2(m, na.rm = TRUE),
             Q3 = q[, 4], Max = q[, 5], NAcount = rowCounts(is.na(m)),
             st.dev = rowSds(m, na.rm = TRUE), diff.Q3Q1 = iqr,
             outliers.coef2 = out_k(2), outliers.coef3 = out_k(3))
}

meth_class <- function(mm) {
  n   <- rowSums(!is.na(mm))
  nL2 <- rowSums(mm < 0.2, na.rm = TRUE); nG8 <- rowSums(mm > 0.8, na.rm = TRUE)
  nG4 <- rowSums(mm > 0.4, na.rm = TRUE); nL6 <- rowSums(mm < 0.6, na.rm = TRUE)
  fifelse(nL2 == n, "L", fifelse(nG8 == n, "H",
  fifelse(nL2 == 0 & nG8 == 0, "M", fifelse(nG4 == n, "HM",
  fifelse(nL6 == n, "LM", "R")))))
}

is_bimodal <- function(mm) {                     # Tian et al. 2016
  n <- rowSums(!is.na(mm))
  rowSums(mm < 0.2, na.rm = TRUE) / n > 0.4 & rowSums(mm > 0.8, na.rm = TRUE) / n > 0.4
}

# ----------------------------------------------------------------
# 1. do these definitions reproduce the lab's own TCGA columns?
# ----------------------------------------------------------------
if (length(tcga_dir)) {
  rule_chk <- rbindlist(lapply(c("Normal", "Tumor"), function(tissue) {
    tg <- fread(file.path(tcga_dir[1],
      sprintf("Sorted.BRCA.53Alive.%s.380355cg.75col.May28.2026.txt", tissue)))
    sc <- 5:(which(names(tg) == "Min") - 1)
    tm <- as.matrix(tg[, ..sc])
    mine <- six_num(tm)
    qc <- c("Min", "Q1", "Median", "Q3", "Max")
    res <- data.table(tissue, sites = nrow(tm), samples = ncol(tm),
      max_abs_diff_quantiles = max(abs(as.matrix(mine[, ..qc]) - as.matrix(tg[, ..qc])), na.rm = TRUE),
      coef2_match = mean(mine$outliers.coef2 == tg$outliers.coef2),
      coef2_mismatch_sites = sum(mine$outliers.coef2 != tg$outliers.coef2),
      coef3_match = mean(mine$outliers.coef3 == tg$outliers.coef3),
      methy_state_match = mean(meth_class(tm) == tg$methy.state),
      bimodal_sites = sum(is_bimodal(tm)))
    rm(tg, tm, mine); invisible(gc()); res
  }))
  fwrite(rule_chk, file.path(out_dir, "Task39_TCGA_RuleCheck.csv"))
  cat("Definitions against the TCGA 53 Alive source columns:\n"); print(rule_chk)
  stopifnot(rule_chk$max_abs_diff_quantiles < 1e-6, rule_chk$coef2_match > 0.9999,
            rule_chk$coef3_match > 0.9999, rule_chk$methy_state_match == 1)
} else cat("TCGA source matrices not found; definition check skipped\n")

# ----------------------------------------------------------------
# 2. read, and check that pos is the C of every forward-strand CpG
# ----------------------------------------------------------------
d <- fread(in_file, na.strings = "NA")
samp <- grep("^P[0-9]+NF$", names(d), value = TRUE)
stopifnot(length(samp) == 18, !is.unsorted(d$pos, strictly = TRUE))
m <- as.matrix(d[, ..samp])
n_site <- nrow(m)
cat("\nsites:", n_site, " samples:", ncol(m), "\n")
cat("value range:", range(m, na.rm = TRUE), "\n")

if (dir.exists(seq_dir)) {
  wins <- list.files(seq_dir, pattern = "_pm500\\.txt$", full.names = TRUE)
  win_ok <- vapply(wins, function(f) {
    start <- as.integer(sub(".*_(\\d+)_pm500\\.txt$", "\\1", f)) - 500
    ln <- readLines(f, warn = FALSE)
    sq <- toupper(gsub("[^A-Za-z]", "", paste(ln[!startsWith(ln, ">")], collapse = "")))
    p <- d$pos[d$pos >= start & d$pos < start + nchar(sq) - 1]
    cg <- gregexpr("CG", sq, fixed = TRUE)[[1]] + start - 1
    all(substring(sq, p - start + 1, p - start + 2) == "CG") && setequal(p, cg)
  }, logical(1))
  cat("hg38 windows where every CpG is listed, at its C:", sum(win_ok), "of", length(win_ok), "\n")
  stopifnot(all(win_ok))
}

# ----------------------------------------------------------------
# 3. six-number summary, IQR outliers, NA count (all sites)
# ----------------------------------------------------------------
summ <- cbind(data.table(chr = d$chr, pos = d$pos), six_num(m))
summ[NAcount == 18, c("Mean", "outliers.coef2", "outliers.coef3") := NA]
num <- setdiff(names(summ), c("chr", "pos", "NAcount", "outliers.coef2", "outliers.coef3"))
summ[, (num) := lapply(.SD, function(v) signif(v, 6)), .SDcols = num]
fwrite(summ, file.path(out_dir, "Task39_RRBS_SiteSummary.csv.gz"))

na_tab <- summ[, .N, by = NAcount][order(NAcount)]
na_tab[, share := r4(N / n_site)]
fwrite(na_tab, file.path(out_dir, "Task39_RRBS_NAcount.csv"))
cat("\nNA count distribution:\n"); print(na_tab)
complete <- which(summ$NAcount == 0)
cat("\nCG sites with 0 NA:", length(complete), "of", n_site,
    sprintf("(%.2f%%)\n", 100 * length(complete) / n_site))

# ----------------------------------------------------------------
# 4. methylation states (complete sites)
# ----------------------------------------------------------------
# cross-check the vectorised rule against the original per-row function
meth.class.orig <- function(x) { x <- x[!is.na(x)]; n <- length(x)
  if (sum(x < .2) == n) "L" else if (sum(x > .8) == n) "H"
  else if (sum(x < .2) == 0 & sum(x > .8) == 0) "M"
  else if (sum(x > .4) == n) "HM" else if (sum(x < .6) == n) "LM" else "R" }
chk <- sample(complete, 5000)
stopifnot(identical(meth_class(m[chk, ]), apply(m[chk, ], 1, meth.class.orig)))

cm <- m[complete, ]
cs <- data.table(chr = d$chr[complete], pos = d$pos[complete],
                 state = factor(meth_class(cm), levels = states))
cs <- cbind(cs, summ[complete, .(Median, Mean, diff.Q3Q1, outliers.coef2, outliers.coef3)])
fwrite(cs, file.path(out_dir, "Task39_RRBS_CompleteSiteStates.csv.gz"))

bim <- is_bimodal(cm)
st_tab <- cs[, .(N = .N, median_IQR = r4(median(diff.Q3Q1)),
                 bimodal_Tian2016 = sum(bim[.I])), keyby = state]
st_tab[, share := r4(N / sum(N))]
fwrite(st_tab, file.path(out_dir, "Task39_RRBS_StateCounts.csv"))
cat("\nStates of complete sites:\n"); print(st_tab)

# ----------------------------------------------------------------
# 5. how many IQR-2 / IQR-3 outliers, against what continuous data give
# ----------------------------------------------------------------
oc_row <- function(scope, x) x[, .(scope = scope, sites = .N,
  IQR_zero = r4(mean(diff.Q3Q1 == 0)),
  any_coef2 = r4(mean(outliers.coef2 > 0)), any_coef3 = r4(mean(outliers.coef3 > 0)),
  any_coef2_if_IQR_pos = r4(mean(outliers.coef2[diff.Q3Q1 > 0] > 0)),
  flagged_with_IQR_zero = r4(mean(diff.Q3Q1[outliers.coef2 > 0] == 0)),
  sites_coef2 = sum(outliers.coef2 > 0), sites_coef3 = sum(outliers.coef3 > 0))]
sim_norm <- six_num(matrix(rnorm(18 * N_SIM), ncol = 18))
oc <- rbindlist(c(
  lapply(states, function(s) oc_row(paste("complete,", s), cs[state == s])),
  list(oc_row("complete, all states", cs),
       oc_row("all sites with >= 4 values", summ[NAcount <= 14]),
       oc_row("continuous Normal data, 18 samples (simulated)", sim_norm))))
fwrite(oc, file.path(out_dir, "Task39_RRBS_OutlierCounts.csv"))
cat("\nIQR-2 / IQR-3 outliers (share of sites with at least one):\n"); print(oc)

# ----------------------------------------------------------------
# 6. value granularity: exact 0 / 1, cut-off values, implied depth
# ----------------------------------------------------------------
v <- as.vector(cm)
mid <- v[v > 0 & v < 1]
# smallest denominator k with v = j/k (values carry 6 significant
# digits); a LOWER bound on read depth, e.g. 0.5 could be 1/2 or 20/40
den <- rep(NA_integer_, length(mid))
for (kk in 2:1000) {
  todo <- which(is.na(den)); if (!length(todo)) break
  hit <- abs(mid[todo] * kk - round(mid[todo] * kk)) < 1e-6 * kk
  den[todo[hit]] <- kk
}
edge_site <- rowSums(abs(cm - 0.4) < 1e-9 | abs(cm - 0.8) < 1e-9) > 0
dq <- quantile(den, c(.1, .25, .5, .75, .9), na.rm = TRUE)
gran <- data.table(
  metric = c("values (complete sites)", "share exactly 0", "share exactly 1",
             "share exactly 0.4", "share exactly 0.8",
             "sites holding a 0.4 or 0.8 (level and state cut-offs disagree)",
             "share of such sites",
             "0<v<1 values with a denominator <= 1000",
             paste("implied depth lower bound,", names(dq))),
  value = c(length(v), r4(mean(v == 0)), r4(mean(v == 1)),
            r4(mean(abs(v - 0.4) < 1e-9)), r4(mean(abs(v - 0.8) < 1e-9)),
            sum(edge_site), r4(mean(edge_site)),
            r4(mean(!is.na(den))), as.numeric(dq)))
fwrite(gran, file.path(out_dir, "Task39_RRBS_ValueGranularity.csv"))
cat("\nValue granularity, complete sites:\n"); print(gran)

# ----------------------------------------------------------------
# 7. neighbouring complete CpGs: state patterns
# ----------------------------------------------------------------
k <- nrow(cs)
pr <- data.table(pos1 = cs$pos[-k], pos2 = cs$pos[-1],
                 s1 = cs$state[-k], s2 = cs$state[-1])
pr[, gap := pos2 - pos1]
# were the two complete CpGs also adjacent among ALL chr22 CpGs?
pr[, adjacent_in_file := (match(pos2, d$pos) - match(pos1, d$pos)) == 1L]
cat("\nconsecutive complete pairs:", nrow(pr), " also adjacent among all chr22 CpGs:",
    sprintf("%.1f%%\n", 100 * mean(pr$adjacent_in_file)))
gap_breaks <- c(0, 10, 50, 100, 200, 500, 1000, 10000, Inf)
pr[, gap_bin := cut(gap, gap_breaks, right = TRUE, dig.lab = 6)]

# chance agreement two ways: chromosome-wide state mix, and the state
# mix of the pairs in that distance bin (near pairs are mostly L, far
# pairs mostly R / HM, so the two differ)
chance_global <- sum(prop.table(table(cs$state))^2)
by_gap <- pr[, {
  p1 <- prop.table(table(s1)); p2 <- prop.table(table(s2))
  pool <- prop.table(table(factor(c(as.character(s1), as.character(s2)), states)))
  c(list(pairs = .N, same_state = r4(mean(s1 == s2)),
         chance_global = r4(chance_global), chance_matched = r4(sum(p1 * p2))),
    as.list(setNames(r4(as.numeric(pool)), paste0("share_", states))))
}, keyby = gap_bin]
fwrite(by_gap, file.path(out_dir, "Task39_RRBS_StateAgreementByGap.csv"))
cat("\nSame-state rate of consecutive complete CpGs by distance (bp), with chance\n",
    "and the state mix of the pairs in each bin:\n"); print(by_gap)

trans <- pr[gap <= 1000, as.data.table(table(from = s1, to = s2))]
trans[, `:=`(from = factor(from, states), to = factor(to, states))]
trans[, row_share := r4(N / sum(N)), by = from]
fwrite(trans, file.path(out_dir, "Task39_RRBS_StateTransitions_1kb.csv"))
cat("\nTransition probabilities, neighbours within 1 kb (row = first CpG):\n")
print(dcast(trans, from ~ to, value.var = "row_share"))
cat("\nShare of neighbours within 1 kb that keep the first CpG's state:\n")
print(trans[from == to, .(state = from, stay = row_share)])

ord <- c(L = 1, LM = 2, M = 3, HM = 4, H = 5)
chg <- pr[gap <= 1000 & s1 != s2]
chg[, kind := fifelse(s1 == "R" | s2 == "R", "to or from R",
              fifelse(abs(ord[as.character(s1)] - ord[as.character(s2)]) == 1,
                      "one step (L-LM, LM-M, M-HM, HM-H)", "jump of two or more steps"))]
kinds <- chg[, .(pairs = .N), keyby = kind][, share := r4(pairs / sum(pairs))]
fwrite(kinds, file.path(out_dir, "Task39_RRBS_StateChangeKinds_1kb.csv"))
cat("\nNeighbours within 1 kb whose states differ, by kind of change:\n"); print(kinds)

# runs: maximal stretches of consecutive complete CpGs, each within 1 kb
# of the next, that share one state
cs[, brk := c(TRUE, diff(pos) > 1000 | state[-1] != state[-.N])]
cs[, run := cumsum(brk)]
runs <- cs[, .(state = state[1], start = pos[1], end = pos[.N], n_cpg = .N), by = run]
run_tab <- runs[, .(runs = .N, mean_len = r4(mean(n_cpg)), max_len = max(n_cpg),
                    runs_ge5 = sum(n_cpg >= 5)), keyby = state]
fwrite(run_tab, file.path(out_dir, "Task39_RRBS_StateRuns.csv"))
fwrite(runs[order(-n_cpg)][1:50], file.path(out_dir, "Task39_RRBS_LongestRuns.csv"))
cat("\nSame-state runs (<= 1 kb spacing):\n"); print(run_tab)

# ----------------------------------------------------------------
# 8. 10-level discretisation and pairwise similarity
# ----------------------------------------------------------------
lev <- pmin(floor(round(cm * 10, 8)) + 1L, 10L)
storage.mode(lev) <- "integer"
X <- lev[-k, , drop = FALSE]; Y <- lev[-1, , drop = FALSE]
BX <- cm[-k, , drop = FALSE]; BY <- cm[-1, , drop = FALSE]

rowcor <- function(a, b) {
  a <- a - rowMeans(a); b <- b - rowMeans(b)
  den <- sqrt(rowSums(a^2) * rowSums(b^2))
  ifelse(den > 0, rowSums(a * b) / den, NA_real_)
}
rowrank <- function(a) t(apply(a, 1, rank))            # average ranks

# quadratic-weighted kappa, in its centred (Lin's concordance) form:
# 1 - mean((a-b)^2) / E_indep[(a-b)^2] = 2 cov / (var_a + var_b + (ma-mb)^2)
kappa_w <- function(a, b) {
  ma <- rowMeans(a); mb <- rowMeans(b)
  ac <- a - ma; bc <- b - mb
  den <- rowMeans(ac^2) + rowMeans(bc^2) + (ma - mb)^2
  ifelse(den > 0, 2 * rowMeans(ac * bc) / den, NA_real_)
}

# Chatterjee xi_n(x -> y); ties in x broken at random
xi_one <- function(x, y) {
  n <- length(y)
  l <- n - rank(y, ties.method = "min") + 1               # #{j : y_j >= y_i}
  den <- 2 * sum(l * (n - l))
  if (den == 0) return(NA_real_)
  r <- rank(y, ties.method = "max")[order(x, runif(n))]  # #{j : y_j <= y_i}, x order
  1 - n * sum(abs(diff(r))) / den
}
stopifnot(abs(xi_one(1:18, 1:18) - (1 - 3 / 19)) < 1e-12)  # the n = 18 maximum
xi_sym_avg <- function(x, y, draws = XI_DRAWS) {
  if (sd(x) == 0 || sd(y) == 0) return(c(NA_real_, NA_real_))
  v <- replicate(draws, max(xi_one(x, y), xi_one(y, x)))
  c(mean(v), sd(v))
}

sim <- data.table(pos1 = pr$pos1, pos2 = pr$pos2, gap = pr$gap, gap_bin = pr$gap_bin,
                  s1 = pr$s1, s2 = pr$s2)
sim[, agree        := rowMeans(X == Y)]
sim[, agree1       := rowMeans(abs(X - Y) <= 1)]
sim[, L1           := rowMeans(abs(X - Y))]
sim[, L2beta       := sqrt(rowMeans((BX - BY)^2))]
sim[, kappa_w      := kappa_w(X, Y)]
sim[, pearson_beta := rowcor(BX, BY)]
sim[, spearman     := rowcor(rowrank(X), rowrank(Y))]
t0 <- Sys.time()
xs <- vapply(seq_len(nrow(X)), function(i) xi_sym_avg(X[i, ], Y[i, ]), numeric(2))
sim[, `:=`(xi = xs[1, ], xi_tiebreak_sd = xs[2, ])]
cat("\nxi (mean of", XI_DRAWS, "tie-breaks) computed in",
    round(as.numeric(Sys.time() - t0, units = "secs")), "s\n")
cat("SD of xi between single tie-break draws: median", r4(median(sim$xi_tiebreak_sd, na.rm = TRUE)),
    " 90th percentile", r4(quantile(sim$xi_tiebreak_sd, .9, na.rm = TRUE)),
    "; the averaged xi has about 1/sqrt(", XI_DRAWS, ") of that\n")
sim[, constant_any := rowSds(X) == 0 | rowSds(Y) == 0]
sim[, pair_type := factor(fifelse(s1 == s2, paste0(s1, "-", s2), "different"),
                          levels = c(paste0(states, "-", states), "different"))]

meas <- c("agree", "agree1", "L1", "L2beta", "kappa_w", "pearson_beta", "spearman", "xi")
sim_out <- copy(sim)
sim_out[, c(meas, "xi_tiebreak_sd") := lapply(.SD, r4), .SDcols = c(meas, "xi_tiebreak_sd")]
fwrite(sim_out, file.path(out_dir, "Task39_RRBS_PairSimilarity.csv.gz"))

# ----------------------------------------------------------------
# 9. summaries of the similarity measures
# ----------------------------------------------------------------
can_be_na <- c("kappa_w", "pearson_beta", "spearman", "xi")
pair_summary <- function(dt, by) {
  med <- dt[, c(list(pairs = .N), lapply(.SD, function(v) r4(median(v, na.rm = TRUE)))),
            keyby = by, .SDcols = meas]
  und <- dt[, c(lapply(.SD, function(v) r4(mean(is.na(v)))),
                list(kappa_zero = r4(mean(abs(kappa_w) < 1e-12, na.rm = TRUE)))),
            keyby = by, .SDcols = can_be_na]
  setnames(und, can_be_na, paste0("undef_", can_be_na))
  med[und, on = by]
}

sim_gap <- pair_summary(sim, "gap_bin")
fwrite(sim_gap, file.path(out_dir, "Task39_RRBS_SimilarityByGap.csv"))
cat("\nMedian similarity of consecutive complete CpGs, by distance, and the share of\n",
    "pairs where each correlation-type measure is undefined:\n"); print(sim_gap)

sim_state <- pair_summary(sim[gap <= 1000], "pair_type")
fwrite(sim_state, file.path(out_dir, "Task39_RRBS_SimilarityByStatePair_1kb.csv"))
cat("\nWithin 1 kb, by state pair (kappa_zero = share of defined kappa_w that are 0):\n")
print(sim_state)

# does similarity fall with distance WITHIN one state, or only because
# the state mix changes with distance?
by_state_gap <- pair_summary(sim[s1 == s2 & s1 != "M"], c("pair_type", "gap_bin"))
by_state_gap[, pair_type := droplevels(pair_type)]
fwrite(by_state_gap, file.path(out_dir, "Task39_RRBS_SimilarityByGapWithinState.csv"))
cat("\nMedian agreement within one state, by distance:\n")
print(dcast(by_state_gap, gap_bin ~ pair_type, value.var = "agree"))
cat("\nMedian pearson_beta within one state, by distance:\n")
print(dcast(by_state_gap, gap_bin ~ pair_type, value.var = "pearson_beta"))

near <- sim[gap <= 1000]
all_def <- near[complete.cases(near[, ..meas])]
mc <- round(cor(all_def[, ..meas], method = "spearman"), 2)
fwrite(data.table(measure = rownames(mc), mc), file.path(out_dir, "Task39_RRBS_MeasureCorrelation_1kb.csv"))
cat("\nCorrelation between measures (Spearman), on the", nrow(all_def), "of", nrow(near),
    "pairs within 1 kb where all eight are defined:\n"); print(mc)
cat("state-pair mix of those pairs:\n"); print(round(prop.table(table(droplevels(all_def$pair_type))), 3))

# the email's own illustration, through the same functions; it has no
# betas, so Pearson is taken on its levels here
ex_measures <- function(x, y, draws = 2000) {
  a <- matrix(x, 1); b <- matrix(y, 1)
  xd <- if (sd(x) == 0 || sd(y) == 0) NA_real_ else replicate(draws, max(xi_one(x, y), xi_one(y, x)))
  data.table(agree = mean(x == y), agree1 = mean(abs(x - y) <= 1), L1 = mean(abs(x - y)),
             kappa_w = kappa_w(a, b), pearson_levels = rowcor(a, b),
             spearman = rowcor(rowrank(a), rowrank(b)),
             xi_mean = mean(xd), xi_min = min(xd), xi_max = max(xd))
}
cg2 <- c(1, 2, 1, 1, 2, 1, 1)
ex <- rbind(cbind(case = "email: cg1 = 1 1 1 1 2 1 1", ex_measures(c(1, 1, 1, 1, 2, 1, 1), cg2)),
            cbind(case = "same cg2, cg1 = 1 1 1 1 1 1 1", ex_measures(rep(1, 7), cg2)))
num_ex <- setdiff(names(ex), "case")
ex[, (num_ex) := lapply(.SD, function(v) round(v, 3)), .SDcols = num_ex]
cat("\nThe email's example (cg2 = 1 2 1 1 2 1 1); xi over 2000 tie-breaks:\n"); print(ex)

cat("\nData examples, neighbours within 50 bp differing in exactly one sample:\n")
one_diff <- which(sim$gap <= 50 & rowSums(X != Y) == 1)
for (i in c(one_diff[!sim$constant_any[one_diff]][1], one_diff[sim$constant_any[one_diff]][1])) {
  cat(sprintf("  %d (%s): %s\n  %d (%s): %s\n    agree=%.3f L1=%.3f kappa_w=%s pearson_beta=%s spearman=%s xi=%s\n",
    sim$pos1[i], as.character(sim$s1[i]), paste(X[i, ], collapse = " "),
    sim$pos2[i], as.character(sim$s2[i]), paste(Y[i, ], collapse = " "),
    sim$agree[i], sim$L1[i], format(r4(sim$kappa_w[i])), format(r4(sim$pearson_beta[i])),
    format(r4(sim$spearman[i])), format(r4(sim$xi[i]))))
}
cat("\ndone\n")
