# ================================================================
# Task 66 pilot, step 3b - what did the mHap mate rule hide, for one
#           plasma donor (GSM4502069 / SRR11615795)?
# Date: Oct 1, 2026
#
# This is a post-result sensitivity analysis of the Task 65 limitation.
# It is not a new test, and it computes no across-donor covariance.
#
# INPUT (from Task66.PlasmaFragments.Extract.Oct01.2026.py, gitignored)
#   Data/gse149438_fragment_pilot/SRR11615795_linked.pat.gz
#   Data/gse149438_fragment_pilot/SRR11615795_mhaprule.pat.gz
#   Data/gse149438_fragment_pilot/SRR11615795_fragment_stats.json
#   Data/gse149438_plasma_mhap/SRX8181977.mhap.gz  (mHapBrowser's file,
#     for a pipeline diagnostic only)
#
# WHAT IS MEASURED (rules fixed before the comparison was run)
#   For adjacent hg19 CpGs at most 200 bp apart, the joint counts of
#   SharedReadNoise.Core.R from the two record sets:
#   - Check: the site counts N and M must be identical, because both
#     sets carry the same calls and differ only in linkage. The script
#     stops otherwise.
#   - Shared calls: r = n00 + n01 + n10 + n11. The mHap rule can only
#     lose shared calls, never add them. Hidden share = 1 - sum(r_mhap)
#     / sum(r_linked), by distance bin.
#   - The per-donor noise term c_k (Task 53 formula) under each record
#     set, for pairs with N_i, N_j >= 6 (the Task 59/65 universe floor).
#     c_k = 0 when r = 0. The comparison covers pairs where r >= 2 in
#     the linked set, and pairs with r = 1 are counted separately. The
#     captured share is mean(c_mhap) / mean(c_linked), by distance bin.
#   - Fragment structure (from the extraction JSON): the share of
#     fragments whose mates' CpG ranges are disjoint, by insert size.
#   - Diagnostic against mHapBrowser's records for the same donor:
#     correlation of site depth on shared CpGs, and shared calls by
#     distance. Trimming, aligner (BSMAP versus Bismark) and duplicate
#     marking differ, so close agreement is expected, not exact.
#
# OUTPUTS
#   Results/Task66_Pilot_FragmentStructure.csv
#   Results/Task66_Pilot_HiddenSharing.csv
#   Results/Task66_Pilot_vsMHapBrowser.csv
#   Results/Fig47_PlasmaHiddenSharing.png
# ================================================================

suppressPackageStartupMessages({ library(data.table); library(jsonlite); library(ggplot2) })
pick_dir <- function(...) { for (d in c(...)) if (dir.exists(d)) return(d)
  stop("no candidate directory exists") }
repo_dir <- pick_dir("/mmfs1/home/wln26/Experiments.Outlier.July31.2026",
                     path.expand("~/Desktop/bioinformatics-research"))
res <- file.path(repo_dir, "Results"); pil <- file.path(repo_dir, "Data", "gse149438_fragment_pilot")
source(file.path(repo_dir, "Scripts", "SharedReadNoise.Core.R"))
AUTO <- paste0("chr", 1:22); MINN <- 6L
cpg <- readRDS(file.path(repo_dir, "Data", "hg19_seq", "cpg_positions_hg19.rds")); I0 <- c(0L, cumsum(lengths(cpg))[-22])
BINS <- c(0, 10, 20, 40, 60, 100, 150, 200)

read_patlike <- function(f) fread(cmd = sprintf("gzip -dc '%s'", f), header = FALSE, sep = "\t",
                                  col.names = c("chr", "idx", "pat", "n"), colClasses = c("character", "integer", "character", "integer"))
read_mhap <- function(f) {                       # Task 62 / 65 reader
  m <- fread(cmd = sprintf("gzip -dc '%s'", f), header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "hap", "n", "strand"),
             colClasses = c("character", "integer", "integer", "character", "integer", "character"))[chr %in% AUTO]
  m[, k := match(chr, AUTO)]
  m[, li := { p <- cpg[[k[1]]]; match(start, p) }, by = k]; m[, lj := { p <- cpg[[k[1]]]; match(end, p) }, by = k]
  stopifnot(m[is.na(li) | is.na(lj) | (lj - li + 1L) != nchar(hap), .N] == 0)
  m[, .(chr, idx = I0[k] + li, pat = chartr("01", "TC", hap), n)]
}
pairs_of <- function(p) rbindlist(lapply(1:22, function(k) {
  pk <- p[chr == AUTO[k], .(idx, pat, n)]; if (!nrow(pk)) return(NULL)
  pr <- srn_pat_counts(pk, I0[k], length(cpg[[k]]), cpg[[k]], max_gap = 200L)
  if (!nrow(pr)) return(NULL)
  pr[, `:=`(chrn = k, gap = cpg[[k]][idx - I0[k] + 1L] - cpg[[k]][idx - I0[k]])][]
}))
noise <- function(d) { r <- d$n00 + d$n01 + d$n10 + d$n11
  fifelse(r >= 2 & d$Ni >= 2 & d$Nj >= 2, (r * d$n11 - (d$n10 + d$n11) * (d$n01 + d$n11)) / (d$Ni * d$Nj * pmax(r - 1, 1)), 0) }

# ----------------------------------------------------------------
# 1. joint counts under both record sets, and the identity check
# ----------------------------------------------------------------
L <- pairs_of(read_patlike(file.path(pil, "SRR11615795_linked.pat.gz")))
M <- pairs_of(read_patlike(file.path(pil, "SRR11615795_mhaprule.pat.gz")))
J <- merge(L, M, by = c("chrn", "idx", "gap"), suffixes = c(".l", ".m"))
stopifnot(nrow(J) == nrow(L), nrow(J) == nrow(M),
          all(J$Ni.l == J$Ni.m), all(J$Nj.l == J$Nj.m), all(J$Mi.l == J$Mi.m), all(J$Mj.l == J$Mj.m))
J[, `:=`(r_l = n00.l + n01.l + n10.l + n11.l, r_m = n00.m + n01.m + n10.m + n11.m)]
stopifnot(all(J$r_m <= J$r_l))
cat("checks passed: identical site counts in both record sets for", nrow(J), "pairs; the mHap rule only removes shared calls\n")
J[, `:=`(c_l = noise(.SD[, .(Ni = Ni.l, Nj = Nj.l, n00 = n00.l, n01 = n01.l, n10 = n10.l, n11 = n11.l)]),
         c_m = noise(.SD[, .(Ni = Ni.m, Nj = Nj.m, n00 = n00.m, n01 = n01.m, n10 = n10.m, n11 = n11.m)]),
         bin = cut(gap, BINS))]

# ----------------------------------------------------------------
# 2. hidden shared calls and the missed noise term, by distance
# ----------------------------------------------------------------
summ <- function(d) d[, .(pairs = .N, shared_linked = sum(r_l), shared_mhap = sum(r_m),
                          hidden_share = round(1 - sum(r_m) / max(sum(r_l), 1), 4),
                          pairs_losing_some = round(mean(r_m < r_l), 4), pairs_losing_all = round(mean(r_m == 0 & r_l > 0), 4),
                          pairs_r1_linked = sum(r_l == 1),
                          mean_noise_linked = signif(mean(c_l), 4), mean_noise_mhap = signif(mean(c_m), 4),
                          noise_captured_share = round(mean(c_m) / mean(c_l), 4))]
deep <- J[Ni.l >= MINN & Nj.l >= MINN]
hid <- rbind(cbind(set = "all pairs (both sites called)", bin = "all", summ(J)),
             cbind(set = "all pairs (both sites called)", J[, summ(.SD), by = bin][order(bin)][, bin := as.character(bin)]),
             cbind(set = "both sites >= 6 reads", bin = "all", summ(deep)),
             cbind(set = "both sites >= 6 reads", deep[, summ(.SD), by = bin][order(bin)][, bin := as.character(bin)]), fill = TRUE)
fwrite(hid, file.path(res, "Task66_Pilot_HiddenSharing.csv"))
cat("\nShared calls hidden by the mHap mate rule, and the noise term it captures:\n"); print(hid)

# ----------------------------------------------------------------
# 3. fragment structure
# ----------------------------------------------------------------
js <- fromJSON(file.path(pil, "SRR11615795_fragment_stats.json"))
kinds <- names(js$insert_size_bins_50bp)
fs <- rbindlist(lapply(kinds, function(k) { v <- unlist(js$insert_size_bins_50bp[[k]])
  data.table(kind = k, insert_bin_bp = as.integer(names(v)), fragments = as.integer(v)) }))
fs[, share_of_bin := round(fragments / sum(fragments), 4), by = insert_bin_bp]
tot <- fs[, .(fragments = sum(fragments)), by = kind][, share := round(fragments / sum(fragments), 4)]
fwrite(rbind(cbind(tot, insert_bin_bp = NA_integer_, share_of_bin = NA_real_)[, .(kind, insert_bin_bp, fragments, share_of_bin = share)],
             fs), file.path(res, "Task66_Pilot_FragmentStructure.csv"))
cat("\nFragments by mate structure:\n"); print(tot); cat("\nprocessing counts:\n"); print(unlist(js$counts))

# ----------------------------------------------------------------
# 4. pipeline diagnostic against mHapBrowser's file for this donor
# ----------------------------------------------------------------
B <- pairs_of(read_mhap(file.path(repo_dir, "Data", "gse149438_plasma_mhap", "SRX8181977.mhap.gz")))
site <- function(d) unique(rbind(d[, .(chrn, idx, N = Ni)], d[, .(chrn, idx = idx + 1L, N = Nj)]), by = c("chrn", "idx"))
sm <- merge(site(M), site(B), by = c("chrn", "idx"), suffixes = c("_ours", "_mhapbrowser"))
D <- merge(M[, .(chrn, idx, gap, r_ours = n00 + n01 + n10 + n11)], B[, .(chrn, idx, r_browser = n00 + n01 + n10 + n11)], by = c("chrn", "idx"))
D[, bin := cut(gap, BINS)]
diag <- rbind(data.table(measure = "CpG sites called in our mHap-rule records", value = nrow(site(M))),
              data.table(measure = "CpG sites called in mHapBrowser's records", value = nrow(site(B))),
              data.table(measure = "sites in both", value = nrow(sm)),
              data.table(measure = "Pearson r of site depth, sites in both", value = round(cor(sm$N_ours, sm$N_mhapbrowser), 4)),
              data.table(measure = "median depth ratio, ours / mHapBrowser", value = round(median(sm$N_ours / sm$N_mhapbrowser), 4)),
              D[, .(measure = paste("shared calls ratio ours / mHapBrowser, gap", bin), value = round(sum(r_ours) / sum(r_browser), 4)), by = bin][order(bin)][, .(measure, value)])
fwrite(diag, file.path(res, "Task66_Pilot_vsMHapBrowser.csv"))
cat("\nDiagnostic against mHapBrowser's records:\n"); print(diag)

# ----------------------------------------------------------------
# 5. figure
# ----------------------------------------------------------------
SURF <- "#fcfcfb"; INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"; GRID <- "#e1e0d9"; BLUE <- "#2a78d6"; ORANGE <- "#eb6834"
f47 <- hid[set == "both sites >= 6 reads" & bin != "all"][, bin := factor(bin, levels = levels(cut(1, BINS)))]
f47 <- melt(f47, id.vars = "bin", measure.vars = c("hidden_share", "noise_captured_share"))
f47[, variable := factor(variable, levels = c("hidden_share", "noise_captured_share"),
                         labels = c("Share of shared calls hidden by the mHap rule", "Share of the noise term the mHap records capture"))]
p47 <- ggplot(f47, aes(bin, value, colour = variable, group = variable)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2.5) +
  scale_colour_manual(values = c(ORANGE, BLUE), name = NULL) + scale_y_continuous(labels = scales::percent) +
  labs(title = "One plasma donor rebuilt from raw reads: what did the mHap records hide?",
       subtitle = "GSM4502069 (SRR11615795), Bismark alignment, mates linked by read name. Pairs with both sites at >= 6 reads.",
       x = "Distance between the two CpGs (bp)", y = NULL,
       caption = "Source: Results/Task66_Pilot_HiddenSharing.csv (Task 66 pilot)") +
  theme_minimal(base_size = 10) +
  theme(plot.background = element_rect(fill = SURF, colour = NA), panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = GRID, linewidth = 0.3), legend.position = "top",
        plot.title = element_text(colour = INK, face = "bold", size = 12), plot.title.position = "plot",
        plot.subtitle = element_text(colour = INK2, size = 9), plot.caption = element_text(colour = MUTED, hjust = 0),
        plot.caption.position = "plot", axis.text = element_text(colour = INK2))
ggsave(file.path(res, "Fig47_PlasmaHiddenSharing.png"), p47, width = 9, height = 5, dpi = 200, bg = SURF)
cat("\ndone\n")
