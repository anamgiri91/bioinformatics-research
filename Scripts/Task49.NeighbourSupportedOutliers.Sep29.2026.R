# ================================================================
# Task 49 - Step 3 of the frozen plan: do outlier flags with a
#           similar flagged neighbour replicate better?
# Date: Sep 29, 2026
#
# Runs section 4 of Results/Task47_FrozenProtocol.md (frozen in
# commit a832294) exactly. The script first checks the freeze.
#
# THE QUESTION (professor's item 1(2))
#   Can neighbours help separate biological outliers from technical
#   ones? A flag that belongs to the person should also show in the
#   same person's tumour. A one-probe technical failure should not.
#
# DATA
#   Flags: Task 13 external-reference OutlierMeth flags, 53 Alive
#   normal samples, autosomal probes whose flags are not missing.
#   Primary set: flags with |beta - median of the 53 normals| >= 0.10.
#   Neighbours: the previous and next probe with flags, hg38 order,
#   same chromosome. A neighbour counts if flagged in the same sample
#   and direction (no floor on the neighbour).
#   Similarity: frozen F+ for each neighbour pair, computed on the 32
#   Dead normal samples only (other people).
#
# OUTCOME
#   Replicated = the same patient's tumour has a flag at the probe in
#   the same direction (Task 13 tumour flags). Chance = the share of
#   the other 52 tumours with a same-direction flag there. Excess =
#   replicated - chance.
#
# RULES
#   run        >= 3 consecutive flagged probes, same sample and
#              direction, each step <= 1,000 bp (epimutacions)
#   adjacency  >= 1 flagged neighbour within 1,000 bp
#   similarity >= 1 flagged neighbour with F+ >= 0.5, any distance
#   support    = the highest F+ among flagged neighbours (0 if none)
#
# PRIMARY TEST (frozen)
#   K = the number of run-rule flags. Mean excess of the top K flags
#   by support (ties split evenly) minus that of the run-rule flags.
#   1,000 bootstrap resamples of patients, seed 20260929. Success if
#   the whole 95% interval is above 0.
#
# SECONDARY (reported, the claim does not rest on them)
#   rules at their own sizes; top K by size |beta - median|; all flags
#   without the floor; excluding Zhou M_general probes; tumour flags
#   replicated in the matched normal; similarity cut-offs 0.3 and 0.7.
#
# OUTPUTS (Results/; aggregates only, no sample-by-probe rows)
#   Task49_Primary.csv        the frozen primary test
#   Task49_Sets.csv           every set in every analysis, with intervals
#   Task49_BySupport.csv      excess replication by support band and run length
#   Fig36_OutlierReplication.png
# ================================================================

suppressPackageStartupMessages({ library(data.table); library(matrixStats); library(ggplot2); library(jsonlite) })

pick_dir <- function(...) { for (d in c(...)) if (dir.exists(d)) return(d)
  stop("no candidate directory exists") }
repo_dir <- pick_dir("/mmfs1/home/wln26/Experiments.Outlier.July31.2026",
                     path.expand("~/Desktop/bioinformatics-research"))
tcga_dir <- pick_dir(path.expand("~/Downloads/tcga_brca_380355cg"),
                     "/home/s_s355/research3.data/TCGA/filter.BRCA.UCEC.Aug17.2022")
res <- file.path(repo_dir, "Results"); r4 <- function(x) round(x, 4)

# ----------------------------------------------------------------
# 0. the freeze is intact
# ----------------------------------------------------------------
proto <- fromJSON(file.path(res, "Task47_FrozenProtocol.json"))
for (f in names(proto$files_sha256))
  stopifnot(digest::digest(file = file.path(repo_dir, f), algo = "sha256") == proto$files_sha256[[f]])
source(file.path(repo_dir, "Scripts", "Task47.FrozenScore.Sep29.2026.R"))
P3 <- proto$step3_outlier_replication
B_BOOT <- P3$uncertainty$B; SEED <- P3$uncertainty$seed
FLOOR <- 0.10; THETA <- 0.5; RUN_GAP <- 1000
stopifnot(B_BOOT == 1000, SEED == 20260929)
cat("freeze verified\n")

# ----------------------------------------------------------------
# 1. data
# ----------------------------------------------------------------
fn <- fread(file.path(res, "Task13_Normal_extref_flags.csv"))
ft <- fread(file.path(res, "Task13_Tumor_extref_flags.csv"))
stopifnot(identical(fn$cgID, ft$cgID))
beta_of <- function(file, prefix, n) {
  d <- fread(file.path(tcga_dir, file), showProgress = FALSE)
  stopifnot(identical(d$Composite.Element.REF, fn$cgID))
  as.matrix(d[, paste0(prefix, seq_len(n)), with = FALSE])
}
bN <- beta_of("Sorted.BRCA.53Alive.Normal.380355cg.75col.May28.2026.txt", "N", 53)
bT <- beta_of("Sorted.BRCA.53Alive.Tumor.380355cg.75col.May28.2026.txt", "T", 53)
bD <- beta_of("Sorted.BRCA.32Dead.Normal.380355cg.54col.May28.2026.txt", "ND", 32)
Fn <- as.matrix(fn[, paste0("N", 1:53), with = FALSE]); Ft <- as.matrix(ft[, paste0("T", 1:53), with = FALSE])

man <- fread(file.path(repo_dir, "Data", "zhou_HM450", "HM450.hg38.manifest.tsv.gz"),
             select = c("Probe_ID", "CpG_chrm", "CpG_beg"))
mk <- fread(file.path(repo_dir, "Data", "zhou_HM450", "HM450.hg38.mask.tsv.gz"), select = c("Probe_ID", "M_general"))
pr <- data.table(row = seq_along(fn$cgID), id = fn$cgID)
pr <- mk[man[pr, on = c(Probe_ID = "id")], on = "Probe_ID"]   # left join: the mask file lists only
pr <- pr[!is.na(CpG_beg) & CpG_chrm %in% paste0("chr", 1:22) & rowSums(is.na(Fn[row, , drop = FALSE])) == 0]
pr[, `:=`(chrn = as.integer(sub("chr", "", CpG_chrm)), pos = CpG_beg + 1L,
          general = M_general %in% TRUE)]                   # masked probes; unlisted = no mask
setorder(pr, chrn, pos)
n <- nrow(pr); rw <- pr$row
same_prev <- c(FALSE, pr$chrn[-1] == pr$chrn[-n]); gap_prev <- c(NA, diff(pr$pos)); gap_prev[!same_prev] <- NA
same_next <- c(same_prev[-1], FALSE); gap_next <- c(gap_prev[-1], NA)
cat("evaluable autosomal probes:", n, "\n")

# F+ between each probe and the next, from the 32 Dead normals
ip <- which(same_next)
sc <- score_F_plus(bD[rw[ip], ], bD[rw[ip + 1], ], gap_next[ip])
F_next <- numeric(n); F_next[ip] <- sc$F_plus
F_prev <- c(0, F_next[-n])
cat("F+ for", length(ip), "neighbour pairs; median", r4(median(sc$F_plus)), "; share >= 0.5:",
    r4(mean(sc$F_plus >= THETA)), "\n")

# ----------------------------------------------------------------
# 2. build the flag table for one direction of the question
# ----------------------------------------------------------------
run_len <- function(M, link) {                # run length of each flagged cell (0 if not flagged)
  out <- matrix(0L, nrow(M), ncol(M))
  for (j in seq_len(ncol(M))) {
    x <- M[, j]; grp <- cumsum(!(x & c(FALSE, x[-length(x)]) & link))
    rl <- tabulate(grp[x], nbins = max(grp)); out[x, j] <- rl[grp[x]]
  }
  out
}
flag_table <- function(Fsrc, Frep, bsrc) {
  med <- rowMedians(bsrc[rw, ]); dev <- abs(bsrc[rw, ] - med)
  link <- same_prev & !is.na(gap_prev) & gap_prev <= RUN_GAP
  rbindlist(lapply(c(1L, -1L), function(d) {
    M <- Fsrc[rw, ] == d; R <- Frep[rw, ] == d
    Mp <- rbind(FALSE, M[-n, , drop = FALSE]) & same_prev
    Mn <- rbind(M[-1, , drop = FALSE], FALSE) & same_next
    support <- pmax(F_prev * Mp, F_next * Mn)
    adj <- (Mp & !is.na(gap_prev) & gap_prev <= RUN_GAP) | (Mn & !is.na(gap_next) & gap_next <= RUN_GAP)
    rl <- run_len(M, link)
    chance_all <- rowSums(R)
    w <- which(M, arr.ind = TRUE)
    rep <- R[w]; chance <- (chance_all[w[, 1]] - rep) / (ncol(R) - 1)
    data.table(patient = w[, 2], probe = w[, 1], direction = d, magnitude = dev[w], support = support[w],
               adj = adj[w], run = rl[w], replicated = rep, chance = chance, excess = rep - chance,
               general = pr$general[w[, 1]])
  }))
}

# ----------------------------------------------------------------
# 3. the sets and the bootstrap
# ----------------------------------------------------------------
top_mean <- function(s, y, K) {
  t <- -sort(-s, partial = K)[K]
  above <- s > t; tied <- s == t
  (sum(y[above]) + (K - sum(above)) * mean(y[tied])) / K
}
stats_of <- function(ft) {                    # every set, from one table of flags
  K <- sum(ft$run >= 3)
  c(all = mean(ft$excess), run = mean(ft$excess[ft$run >= 3]), adjacency = mean(ft$excess[ft$adj]),
    similarity = mean(ft$excess[ft$support >= THETA]),
    similarity_0.3 = mean(ft$excess[ft$support >= 0.3]), similarity_0.7 = mean(ft$excess[ft$support >= 0.7]),
    topK_support = top_mean(ft$support, ft$excess, K), topK_size = top_mean(ft$magnitude, ft$excess, K),
    diff_primary = top_mean(ft$support, ft$excess, K) - mean(ft$excess[ft$run >= 3]),
    diff_size = top_mean(ft$support, ft$excess, K) - top_mean(ft$magnitude, ft$excess, K))
}
SIZES <- function(ft) c(all = nrow(ft), run = sum(ft$run >= 3), adjacency = sum(ft$adj),
                        similarity = sum(ft$support >= THETA), similarity_0.3 = sum(ft$support >= 0.3),
                        similarity_0.7 = sum(ft$support >= 0.7), topK_support = sum(ft$run >= 3),
                        topK_size = sum(ft$run >= 3), diff_primary = NA, diff_size = NA)
boot <- function(ft, npt = 53) {
  set.seed(SEED)
  by_pt <- split(seq_len(nrow(ft)), factor(ft$patient, levels = seq_len(npt)))
  t(vapply(seq_len(B_BOOT), function(b) stats_of(ft[unlist(by_pt[sample.int(npt, npt, replace = TRUE)], use.names = FALSE)]),
           stats_of(ft)))
}
analyse <- function(ft, label) {
  est <- stats_of(ft); bt <- boot(ft)
  data.table(analysis = label, set = names(est), flags = SIZES(ft)[names(est)],
             replicated = c(all = mean(ft$replicated), run = mean(ft$replicated[ft$run >= 3]),
                            adjacency = mean(ft$replicated[ft$adj]), similarity = mean(ft$replicated[ft$support >= THETA]),
                            similarity_0.3 = mean(ft$replicated[ft$support >= 0.3]),
                            similarity_0.7 = mean(ft$replicated[ft$support >= 0.7]),
                            topK_support = NA, topK_size = NA, diff_primary = NA, diff_size = NA)[names(est)],
             excess = est, lo95 = apply(bt, 2, quantile, 0.025, na.rm = TRUE),
             hi95 = apply(bt, 2, quantile, 0.975, na.rm = TRUE))
}

t0 <- Sys.time()
fN <- flag_table(Fn, Ft, bN)                 # normal flags, replicated in the matched tumour
fT <- flag_table(Ft, Fn, bT)                 # secondary: tumour flags, replicated in the matched normal
cat("normal flags (autosomes):", nrow(fN), " passing the floor:", sum(fN$magnitude >= FLOOR),
    " | tumour flags:", nrow(fT), "\n")
prim <- fN[magnitude >= FLOOR]
sets <- rbindlist(list(
  analyse(prim, "primary: normal flags >= 0.10"),
  analyse(fN, "secondary: all normal flags, no floor"),
  analyse(prim[general == FALSE], "secondary: primary set without Zhou M_general probes"),
  analyse(fT[magnitude >= FLOOR], "secondary: tumour flags >= 0.10, replicated in the matched normal")))
for (v in c("replicated", "excess", "lo95", "hi95")) set(sets, j = v, value = r4(sets[[v]]))
fwrite(sets, file.path(res, "Task49_Sets.csv"))
cat("bootstrap done in", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")

primary <- sets[analysis == "primary: normal flags >= 0.10" & set %in% c("run", "topK_support", "diff_primary")]
pt <- data.table(K = primary[set == "run", flags],
                 run_rule_excess = primary[set == "run", excess],
                 topK_support_excess = primary[set == "topK_support", excess],
                 difference = primary[set == "diff_primary", excess],
                 lo95 = primary[set == "diff_primary", lo95], hi95 = primary[set == "diff_primary", hi95])
pt[, success := lo95 > 0]
fwrite(pt, file.path(res, "Task49_Primary.csv"))
cat("\nPRIMARY (frozen): top K by similarity support vs the run rule, excess replication in the matched tumour\n")
print(pt)
cat("\nAll sets (excess = replicated minus chance; 95% patient-bootstrap interval):\n")
print(sets[, .(analysis = substr(analysis, 1, 38), set, flags, replicated, excess, lo95, hi95)], nrows = 100)

# ----------------------------------------------------------------
# 4. descriptive: excess by support band and by run length
# ----------------------------------------------------------------
prim[, `:=`(sband = cut(support, c(-Inf, 0, 0.25, 0.5, 0.75, 1),
                        labels = c("0 (no flagged neighbour)", "(0, 0.25]", "(0.25, 0.5]", "(0.5, 0.75]", "(0.75, 1]")),
            rband = factor(pmin(run, 4L), levels = 1:4, labels = c("1", "2", "3", "4 or more")))]
summ_by <- function(col, what) prim[, .(flags = .N, replicated = r4(mean(replicated)), chance = r4(mean(chance)),
                                        excess = r4(mean(excess))), by = .(band = get(col))][, what := what]
bands <- rbind(summ_by("sband", "support"), summ_by("rband", "run length"))
setcolorder(bands, "what"); setorder(bands, what, band)
fwrite(bands, file.path(res, "Task49_BySupport.csv"))
cat("\nPrimary set by support band and by run length:\n"); print(bands)

# EXPLORATORY (not in the frozen plan): why is the patient interval wide?
# How concentrated are the flags, and how much does one patient move the
# primary difference (leave one patient out)?
pt_n <- prim[, .N, by = patient][order(-N)]
cat("\nExploratory: the 4 patients with the most primary-set flags hold", r4(sum(pt_n$N[1:4]) / nrow(prim)),
    "of them; the top patient alone holds", r4(pt_n$N[1] / nrow(prim)), "\n")
jk <- vapply(seq_len(53), function(k) stats_of(prim[patient != k])[["diff_primary"]], 0)
cat("Exploratory: primary difference with one patient left out ranges from", r4(min(jk)), "to", r4(max(jk)),
    "; below 0 when leaving out", sum(jk < 0), "of 53 patients\n")
fwrite(data.table(analysis = "exploratory: leave one patient out, primary difference",
                  min = r4(min(jk)), median = r4(median(jk)), max = r4(max(jk)), n_below_0 = sum(jk < 0),
                  top4_patient_share_of_flags = r4(sum(pt_n$N[1:4]) / nrow(prim))),
       file.path(res, "Task49_Exploratory_Concentration.csv"))

# ----------------------------------------------------------------
# 5. figure
# ----------------------------------------------------------------
SURF <- "#fcfcfb"; INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"; GRID <- "#e1e0d9"; BLUE <- "#2a78d6"
LAB <- c(all = "All flags", run = "Run rule (3+ in a row, epimutacions)", adjacency = "Adjacency rule (1 neighbour, 1 kb)",
         similarity = "Similarity rule (neighbour F+ >= 0.5)", topK_support = "Top K by similarity support",
         topK_size = "Top K by size alone")
f36 <- sets[set %in% names(LAB) & analysis %in% c("primary: normal flags >= 0.10", "secondary: all normal flags, no floor")]
f36[, `:=`(name = factor(LAB[set], levels = rev(LAB)),
           panel = factor(fifelse(grepl("^primary", analysis), "Flags at least 0.10 from the median (primary)",
                                  "All flags (no floor)"),
                          levels = c("Flags at least 0.10 from the median (primary)", "All flags (no floor)")),
           ours = set %in% c("similarity", "topK_support"))]
p36 <- ggplot(f36, aes(excess, name)) +
  geom_vline(xintercept = 0, colour = MUTED, linewidth = 0.4) +
  geom_segment(aes(x = lo95, xend = hi95, yend = name), colour = "#c3c2b7", linewidth = 1, lineend = "round") +
  geom_point(aes(fill = ours), shape = 21, colour = SURF, stroke = 0.6, size = 3.2) +
  geom_text(aes(label = paste0("n = ", format(flags, big.mark = ",", trim = TRUE))), nudge_y = 0.32, size = 2.8, colour = INK2) +
  scale_fill_manual(values = c(`TRUE` = BLUE, `FALSE` = MUTED), labels = c(`TRUE` = "Uses the similarity score",
                                                                          `FALSE` = "Established or baseline"), name = NULL) +
  facet_wrap(~panel, nrow = 1, scales = "free_x") +
  labs(title = "Do normal-tissue outlier flags show again in the same patient's tumour?",
       subtitle = paste0("Excess replication = share of flags repeated in the matched tumour, minus the share of other ",
                         "patients' tumours flagged there.\n53 TCGA-BRCA patients, OutlierMeth external-reference flags. ",
                         "Lines: 95% intervals from resampling patients.\nThe frozen test is the paired difference, top K by ",
                         "similarity support minus the run rule: ", sprintf("%+.3f", pt$difference), " (95% interval ",
                         sprintf("%.3f", pt$lo95), " to ", sprintf("%.3f", pt$hi95), "), so it ",
                         if (pt$success) "passes." else "narrowly fails."),
       x = "Excess replication in the matched tumour (higher = more person-level)", y = NULL,
       caption = "Source: Results/Task49_Sets.csv (Task 49, frozen plan of Task 47)") +
  theme_minimal(base_size = 10.5) +
  theme(plot.background = element_rect(fill = SURF, colour = NA),
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(colour = GRID, linewidth = 0.3), panel.spacing.x = unit(16, "pt"),
        strip.text = element_text(colour = INK, face = "bold", hjust = 0),
        plot.title = element_text(colour = INK, face = "bold", size = 13), plot.title.position = "plot",
        plot.subtitle = element_text(colour = INK2, size = 9.5), plot.caption = element_text(colour = MUTED, hjust = 0),
        plot.caption.position = "plot", legend.position = "top", legend.justification = "left",
        axis.text = element_text(colour = INK2))
ggsave(file.path(res, "Fig36_OutlierReplication.png"), p36, width = 11, height = 5.5, dpi = 200, bg = SURF)
cat("\ndone\n")
