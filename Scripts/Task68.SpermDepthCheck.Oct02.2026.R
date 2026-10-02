# ================================================================
# Task 68 - Is the sperm cohort deep enough for the three-way split?
#           A counts-only check from the study's own CpG reports.
# Date: Oct 2, 2026
#
# Cluster access is blocked, so this uses the Bismark CpG reports that
# GEO publishes for the two pilot donors Codex chose from metadata
# (GSM5058020 / SRR13602252, low exposure; GSM5058029 / SRR13602261,
# high exposure). Only depth (methylated + unmethylated counts) is
# read. No methylation level, covariance or pair statistic is computed.
#
# RULES (fixed before the data were read)
#   Depth of a CpG: both strands combined (+ at the C, - at the C + 1).
#     Autosomes only.
#   Pairs: adjacent CpGs in the report at most 200 bp apart.
#   A donor is eligible for a pair with probability e(N), where
#     N = min(depth_i, depth_j). e(N) is the chance that a 1/3 : 1/3 : 1/3
#     split of N shared fragments gives every part at least 2 reads, as in
#     the Task 59 / 65 rule. It assumes that close CpGs share their reads.
#   Projection to 52 donors: the pair's eligibility probability is the
#     mean of the two donors' e(N). The pair counts as usable with
#     probability P(Binomial(52, p) >= 20). This is a rough feasibility
#     estimate: it assumes the other 50 donors look like these two.
#
# INPUT  Data/gse165915_depth_check/*.CpG.report.txt.gz (gitignored)
# OUTPUT Results/Task68_SpermDepthCheck.csv
# ================================================================
suppressPackageStartupMessages(library(data.table))
d0 <- "Data/gse165915_depth_check"
files <- c(GSM5058020 = "GSM5058020_18098D-02-11_R_1.CpG.report.txt.gz", GSM5058029 = "GSM5058029_18098D-02-08_R_1.CpG.report.txt.gz")
AUTO <- paste0("chr", 1:22)
e_given_n <- sapply(0:200, function(n) { if (n < 6) return(0); s <- 0   # P(all three parts >= 2 | n)
  for (a in 2:(n - 4)) for (b in 2:(n - a - 2)) { c <- n - a - b; if (c >= 2) s <- s + dmultinom(c(a, b, c), prob = rep(1/3, 3)) }; s })
eN <- function(n) e_given_n[pmin(n, 200) + 1]

depth <- lapply(names(files), function(g) {
  x <- fread(cmd = sprintf("gzip -dc '%s'", file.path(d0, files[[g]])), header = FALSE, sep = "\t",
             select = 1:5, col.names = c("chr", "pos", "strand", "M", "U"))[chr %in% AUTO]
  x[strand == "-", pos := pos - 1L]
  x[, .(N = sum(M + U)), keyby = .(chr, pos)]
})
names(depth) <- names(files)
stopifnot(identical(depth[[1]]$pos, depth[[2]]$pos))          # same CpG list in both reports
D <- depth[[1]][, .(chr, pos, N1 = N)][, N2 := depth[[2]]$N]

site <- rbindlist(lapply(c("N1", "N2"), function(v) data.table(donor = names(files)[match(v, c("N1", "N2"))],
  CpGs = nrow(D), mean_depth = round(mean(D[[v]]), 2), median_depth = median(D[[v]]),
  share_ge6 = round(mean(D[[v]] >= 6), 4), share_ge10 = round(mean(D[[v]] >= 10), 4), share_zero = round(mean(D[[v]] == 0), 4))))
print(site)

D[, nxt := shift(pos, -1L), by = chr]
P <- D[!is.na(nxt) & nxt - pos <= 200L]
P[, `:=`(N1j = D$N1[match(paste(chr, nxt), paste(D$chr, D$pos))], N2j = D$N2[match(paste(chr, nxt), paste(D$chr, D$pos))])]
stopifnot(!anyNA(P$N1j), !anyNA(P$N2j))
P[, `:=`(e1 = eN(pmin(N1, N1j)), e2 = eN(pmin(N2, N2j)))]
stopifnot(!anyNA(P$e1), !anyNA(P$e2))
P[, p := (e1 + e2) / 2]
P[, p := pmin(pmax(p, 0), 1)]                                   # guard rounding just above 1
P[, usable := pbinom(19, 52, p, lower.tail = FALSE)]
pairs <- data.table(pairs_within_200bp = nrow(P),
                    mean_pair_eligibility_GSM5058020 = round(mean(P$e1), 4), mean_pair_eligibility_GSM5058029 = round(mean(P$e2), 4),
                    share_pairs_both_sites_ge6_both_donors = round(mean(pmin(P$N1, P$N1j) >= 6 & pmin(P$N2, P$N2j) >= 6), 4),
                    expected_usable_pairs = round(sum(P$usable)), expected_usable_share = signif(mean(P$usable), 3))
print(pairs)
out <- cbind(site, pairs[rep(1, nrow(site))])
fwrite(out, "Results/Task68_SpermDepthCheck.csv")
cat("done\n")

# ----------------------------------------------------------------
# Where are the projected-usable pairs? (added after the first run
# showed 2.6% of pairs projected usable; still depth only)
#   - depth of those pairs against the cohort mean (about 3x)
#   - overlap with UCSC RepeatMasker hg38 classes, where excess depth
#     often means mis-mapped reads rather than real coverage
# ----------------------------------------------------------------
rm <- fread(cmd = "gzip -dc Data/annotation/rmsk_hg38.txt.gz", header = FALSE, select = c(6, 7, 8, 12),
            col.names = c("chr", "start", "end", "repClass"))[chr %in% AUTO]
rm[, start := start + 1L]                                          # 0-based half-open to 1-based closed
setkey(rm, chr, start, end)
U <- P[usable > 0.5, .(chr, pos, nxt, N1, N1j, N2, N2j)]
U[, `:=`(start = pos, end = nxt)]
ov <- foverlaps(U, rm, by.x = c("chr", "start", "end"), type = "any", mult = "first", nomatch = NA)
U[, repClass := ov$repClass]
U[, mean_min_depth := (pmin(N1, N1j) + pmin(N2, N2j)) / 2]
where <- U[, .(pairs = .N, share = round(.N / nrow(U), 4), median_min_depth = median(mean_min_depth),
               q90_min_depth = quantile(mean_min_depth, 0.9)), by = .(repeat_class = fifelse(is.na(repClass), "no repeat", repClass))][order(-pairs)]
U[, block := paste(chr, pos %/% 1e6)]
cat("\nPairs with projected usable probability > 0.5:", nrow(U), " in", uniqueN(U$block), "1-Mb blocks\n")
print(head(where, 12))
cat("min-site depth of these pairs (mean of the two donors): median", median(U$mean_min_depth), " 90th percentile", quantile(U$mean_min_depth, 0.9), "\n")
fwrite(where, "Results/Task68_SpermDepthCheck_UsablePairs.csv")
cat("done (where)\n")
