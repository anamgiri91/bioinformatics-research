# ================================================================
# Task 51 - Section 5 of the frozen plan: the Task 46 findings on a
#           public RRBS cohort with read counts (29 GTEx colon samples)
# Date: Sep 29, 2026
#
# Runs section 5 of Results/Task47_FrozenProtocol.md (commit a832294)
# on the cohort recorded in commit 4481952, built by Task 50.
#
# DATA
#   GSE233417 colon RRBS, 29 donors, hg19. Autosomal CpGs with at
#   least 10 reads in every sample. Pairs: truly consecutive CpGs
#   (genome CpG index i and i + 1), both kept.
#   Halves: each chromosome split at the median position of its pairs.
#
# TEST (as section 3, frozen)
#   10 patient folds, seed 20260929; outcomes on held-out patients;
#   formula families over the frozen length-scale grid; top 10%; six
#   error levels over the middle 80% of each overlap.
#   H1-H5 as in Task 48. H6 uses the sequencing flags, taken on hg38
#   exactly as frozen: the hg19 sites are lifted to hg38 (UCSC chain);
#   sites that do not lift to one place are left out of H6.
#   H7 (new, frozen): F+ with S weighted per sample by the smaller of
#   the two sites' read counts beats F+ at the same agreement in more
#   than half the cells, genome-wide and in both halves.
#
# OUTPUTS (Results/)
#   Task51_Hypotheses.csv, Task51_Matched.csv, Task51_Frontier.csv,
#   Task51_Frontier_Folds.csv, Task51_Reliability_SpreadMatched.csv,
#   Task51_DistanceDecay.csv (exploratory), Fig37_GTExColon_Frontier.png
#   Data/gtex_colon_rrbs/pair_outcomes.rds (held-out outcomes per pair,
#   for section 6)
# ================================================================

suppressPackageStartupMessages({ library(data.table); library(matrixStats); library(ggplot2); library(jsonlite)
  library(rtracklayer) })

pick_dir <- function(...) { for (d in c(...)) if (dir.exists(d)) return(d)
  stop("no candidate directory exists") }
repo_dir <- pick_dir("/mmfs1/home/wln26/Experiments.Outlier.July31.2026",
                     path.expand("~/Desktop/bioinformatics-research"))
res <- file.path(repo_dir, "Results"); ann_dir <- file.path(repo_dir, "Data", "annotation")
dat <- file.path(repo_dir, "Data", "gtex_colon_rrbs"); r4 <- function(x) round(x, 4)

# ----------------------------------------------------------------
# 0. the freeze is intact
# ----------------------------------------------------------------
proto <- fromJSON(file.path(res, "Task47_FrozenProtocol.json"))
for (f in names(proto$files_sha256))
  stopifnot(digest::digest(file = file.path(repo_dir, f), algo = "sha256") == proto$files_sha256[[f]])
source(file.path(repo_dir, "Scripts", "Task47.FrozenScore.Sep29.2026.R"))
P <- proto$shared_test_design
ELL <- as.numeric(P$ell_grid_bp); TOP <- P$top_fraction; SEED <- 20260929L; NFOLD <- 10L
cat("freeze verified\n")

# ----------------------------------------------------------------
# 1. cohort, pairs, halves
# ----------------------------------------------------------------
cc0 <- readRDS(file.path(dat, "cache_cohort.rds"))
st <- cc0$sites; Bm <- cc0$M / cc0$N; Nr <- cc0$N
rm(cc0); invisible(gc())
stopifnot(!anyNA(Bm), all(Nr >= 10), ncol(Bm) == 29)
nx <- which(st$chrn[-1] == st$chrn[-nrow(st)] & st$idx[-1] == st$idx[-nrow(st)] + 1L)
pairs <- data.table(i = nx, j = nx + 1L, chrn = st$chrn[nx], pos = st$pos[nx], gap = st$pos[nx + 1L] - st$pos[nx])
stopifnot(all(pairs$gap > 0))
pairs[, half := fifelse(pos < median(pos), "first", "second"), by = chrn]
gap <- as.numeric(pairs$gap)
SUBSETS <- list(all = rep(TRUE, nrow(pairs)), first = pairs$half == "first", second = pairs$half == "second")
cat("sites:", nrow(st), " truly consecutive pairs:", nrow(pairs), " gap range:", range(gap),
    "bp, median", median(gap), "\n")

# ----------------------------------------------------------------
# 2. reliability flags on hg38 (lifted), as frozen
# ----------------------------------------------------------------
chain_file <- file.path(tempdir(), "hg19ToHg38.over.chain")
system2("gzip", c("-dc", shQuote(file.path(ann_dir, "hg38_hg19ToHg38.over.chain.gz"))), stdout = chain_file)
# Map the full CG dyad. Its leftmost mapped base remains the target C on
# reverse-chain mappings; lifting the source C alone can point to the target G.
lo <- liftOver(GRanges(paste0("chr", st$chrn), IRanges(st$pos, width = 2), strand = "+"), import.chain(chain_file))
one <- lengths(lo) == 1
g38 <- unlist(lo[one])
site38 <- data.table(k = which(one), chr = as.character(seqnames(g38)), pos = start(g38),
                     width = width(g38), strand = as.character(strand(g38)))
site38 <- site38[chr %in% paste0("chr", 1:22) & width == 2]
cat("sites lifted to one hg38 place:", nrow(site38), "of", nrow(st), "\n")
cat("reverse-orientation CpG dyads:", sum(site38$strand == "-"), "\n")
q38 <- GRanges(site38$chr, IRanges(site38$pos, width = 1))
# One chromosome at a time: rtracklayer's bigBed import with a `which` that
# spans several chromosomes returns features for the first chromosome only
# (checked 2026-09-29: chr1 alone came back), which marked every site outside
# chr1 as not unique in the first run of this script.
site38[, unique50 := FALSE]
for (ch in unique(site38$chr)) {
  ix <- which(site38$chr == ch)
  u <- import(BigBedFile(file.path(ann_dir, "hg38_k50.C2T-Converted.bb")),
              which = GRanges(ch, IRanges(1L, max(site38$pos[ix]) + 1L)))
  set(site38, ix, "unique50", overlapsAny(q38[ix], u))
}
rm(u); invisible(gc())
uq <- site38[, .(share = mean(unique50)), by = chr]
print(uq[order(share)][1:3])
stopifnot(all(uq$share > 0.5))              # a whole chromosome under 50% unique means a lookup failure
bw <- BigWigFile(file.path(ann_dir, "hg38_k50.Bismap.MultiTrackMappability.bw"))
mp <- import(bw, which = q38)
hit <- findOverlaps(q38, mp)
site38[, map50 := 0][queryHits(hit), map50 := mp$score[subjectHits(hit)]]; rm(mp); invisible(gc())
rmk <- fread(file.path(ann_dir, "rmsk_hg38.txt.gz"), header = FALSE, select = c(6, 7, 8, 12),
             col.names = c("chrom", "start", "end", "repClass"))[repClass == "SINE"]
site38[, SINE := overlapsAny(q38, GRanges(rmk$chrom, IRanges(rmk$start + 1L, rmk$end)))]; rm(rmk)
fl <- data.table(k = seq_len(nrow(st)))[site38, on = "k", `:=`(not_unique = !i.unique50, poor_mapping = i.map50 < 1, SINE = i.SINE)]
for (v in c("not_unique", "poor_mapping", "SINE")) pairs[, (paste0("flag_", v)) := fl[[v]][i] | fl[[v]][j]]
invalid_flag <- is.na(fl$not_unique[pairs$i]) | is.na(fl$not_unique[pairs$j])
for (v in c("not_unique", "poor_mapping", "SINE")) set(pairs, which(invalid_flag), paste0("flag_",v), NA)
cat("flagged pairs  not unique:", sum(pairs$flag_not_unique, na.rm = TRUE), " poor mapping:",
    sum(pairs$flag_poor_mapping, na.rm = TRUE), " SINE:", sum(pairs$flag_SINE, na.rm = TRUE),
    " no flag (not lifted):", sum(is.na(pairs$flag_SINE)), "\n")

# ----------------------------------------------------------------
# 3. the test (Task 48 machinery, plus H7)
# ----------------------------------------------------------------
top_mean <- function(s, y, K) {
  t <- -sort(-s, partial = K)[K]
  above <- s > t; tied <- s == t
  (sum(y[above]) + (K - sum(above)) * mean(y[tied])) / K
}
# Reuse the same top-K boundary for all endpoints. This preserves the original
# tie weighting while avoiding three sorts of the same score vector.
top_endpoints <- function(s, A, B, B_legacy, K) {
  t <- -sort(-s, partial = K)[K]
  above <- s > t; tied <- s == t; nt <- K-sum(above)
  avg <- function(y) (sum(y[above])+nt*mean(y[tied]))/K
  list(error=avg(A),gain=avg(B),gain_legacy=avg(B_legacy))
}
outcomes <- function(Xt, Yt, Xk, Yk) {
  mx <- rowMeans(Xt); my <- rowMeans(Yt)
  vx <- rowMeans((Xt - mx)^2); vy <- rowMeans((Yt - my)^2); cxy <- rowMeans((Xt - mx) * (Yt - my))
  bxy <- ifelse(vy > 0, cxy / vy, 0); byx <- ifelse(vx > 0, cxy / vx, 0)
  xhat_old <- pmin(1, pmax(0, mx + bxy * (Yk - my))); yhat_old <- pmin(1, pmax(0, my + byx * (Xk - mx)))
  # Exact constants can have tiny positive floating-point variance. A straight
  # line with a constant predictor must reduce to the training intercept.
  constant_x <- rowMins(Xt) == rowMaxs(Xt); constant_y <- rowMins(Yt) == rowMaxs(Yt)
  vx[constant_x] <- 0; vy[constant_y] <- 0; cxy[constant_x | constant_y] <- 0
  bxy <- ifelse(vy > 0, cxy / vy, 0); byx <- ifelse(vx > 0, cxy / vx, 0)
  xhat <- pmin(1, pmax(0, mx + bxy * (Yk - my))); yhat <- pmin(1, pmax(0, my + byx * (Xk - mx)))
  list(A = rowSums(abs(Xk - Yk)),
       B = rowSums(((abs(Xk - mx) - abs(Xk - xhat)) + (abs(Yk - my) - abs(Yk - yhat))) / 2),
       B_legacy = rowSums(((abs(Xk - mx) - abs(Xk - xhat_old)) + (abs(Yk - my) - abs(Yk - yhat_old))) / 2),
       n = ncol(Xk))
}
FAMILIES <- c("Codex F", "F no distance", "F+", "F guarded", "New C", "New C with dCor", "F+ depth-weighted")
X <- Bm[pairs$i, ]; Y <- Bm[pairs$j, ]; CW <- pmin(Nr[pairs$i, ], Nr[pairs$j, ])
set.seed(SEED); fold_of <- sample(rep(seq_len(NFOLD), length.out = ncol(X)))
A_tot <- B_tot <- B_legacy_tot <- numeric(nrow(X))
t0 <- Sys.time()
fr_folds <- rbindlist(lapply(seq_len(NFOLD), function(f) {
  tr <- which(fold_of != f); te <- which(fold_of == f)
  Xt <- X[, tr]; Yt <- Y[, tr]
  sc <- score_F_plus(Xt, Yt, gap)
  S <- sc$S; E <- sc$E; Ep <- sc$E_plus
  Sw <- 1 - rowSums(CW[, tr] * abs(Xt - Yt)) / rowSums(CW[, tr])
  S_asin <- 1 - rowMeans(abs(asin(sqrt(Xt)) - asin(sqrt(Yt)))) / (pi / 2)
  m <- pmin(rowSums(abs(Xt - rowMedians(Xt)) > 0.1), rowSums(abs(Yt - rowMedians(Yt)) > 0.1))
  D <- pmax(sc$spearman, 0); D[is.na(D)] <- 0; D <- D * m / (m + 3)
  cc <- rowMaxs(Xt) == rowMins(Xt) & rowMaxs(Yt) == rowMins(Yt)
  o <- outcomes(Xt, Yt, X[, te, drop = FALSE], Y[, te, drop = FALSE])
  A_tot <<- A_tot + o$A; B_tot <<- B_tot + o$B
  B_legacy_tot <<- B_legacy_tot + o$B_legacy
  Af <- o$A / o$n; Bf <- o$B / o$n; B_old <- o$B_legacy / o$n
  out <- rbindlist(lapply(ELL, function(l) {
    w <- if (l == 0) as.numeric(gap == 0) else exp(-gap / l); cbar <- mean(w)
    rule <- function(Sx, v) fifelse(cc, Sx, Sx * v)
    fam <- list(`Codex F` = S * (w + (1 - w) * E), `F no distance` = S * (cbar + (1 - cbar) * E),
                `F+` = S * (w + (1 - w) * Ep), `F guarded` = S * (w + (1 - w) * Ep * m / (m + 3)),
                `New C` = rule(S_asin, w + (1 - w) * D), `New C with dCor` = rule(S_asin, w + (1 - w) * E),
                `F+ depth-weighted` = Sw * (w + (1 - w) * Ep))
    rbindlist(lapply(names(SUBSETS), function(sn) {
      ix <- SUBSETS[[sn]]; K <- ceiling(TOP * sum(ix))
      rbindlist(lapply(FAMILIES, function(fm) {
        v <- top_endpoints(fam[[fm]][ix], Af[ix], Bf[ix], B_old[ix], K)
        data.table(subset = sn, family = fm, ell = l, mean_w = cbar,
                   error = v$error, gain = v$gain, gain_legacy = v$gain_legacy)
      }))
    }))
  }))
  cat(sprintf("  fold %d of %d done (%.1f min)\n", f, NFOLD, as.numeric(Sys.time() - t0, units = "mins")))
  out[, fold := f]
  fwrite(out,file.path(res,sprintf("Task51_Fold%02d.csv",f)))
  out
}))
fwrite(fr_folds, file.path(res, "Task51_Frontier_Folds.csv"))
fr <- fr_folds[, .(mean_w = r4(mean_w[1]), error = round(mean(error), 5), gain = round(mean(gain), 5),
                   gain_legacy = round(mean(gain_legacy), 5),
                   error_lo = round(min(error), 5), error_hi = round(max(error), 5),
                   gain_lo = round(min(gain), 5), gain_hi = round(max(gain), 5)), by = .(subset, family, ell)]
fwrite(fr, file.path(res, "Task51_Frontier.csv"))
po <- data.table(chrn = pairs$chrn, idx = st$idx[pairs$i], pos = pairs$pos, gap = pairs$gap, half = pairs$half,
                 heldout_error = A_tot / ncol(X), heldout_gain = B_tot / ncol(X),
                 heldout_gain_legacy = B_legacy_tot / ncol(X),
                 spread = (rowSds(X) + rowSds(Y)) / 2,
                 low_spread = pmin(rowSums(abs(X - rowMedians(X)) > 0.1), rowSums(abs(Y - rowMedians(Y)) > 0.1)) == 0)
for (nm in c("flag_not_unique","flag_poor_mapping","flag_SINE")) po[, (nm) := pairs[[nm]]]
saveRDS(po, file.path(dat, "pair_outcomes.rds"))
fwrite(data.table(pairs=nrow(po),pairs_changed=sum(abs(po$heldout_gain-po$heldout_gain_legacy)>1e-12),
  max_abs_gain_change=max(abs(po$heldout_gain-po$heldout_gain_legacy)),
  mean_gain_change=mean(po$heldout_gain-po$heldout_gain_legacy)),
  file.path(res,"Task51_ConstantPredictorImpact.csv"))

# ----------------------------------------------------------------
# 4. hypotheses
# ----------------------------------------------------------------
QUESTIONS <- list(H1 = c("Codex F", "F no distance", "more"), H2 = c("New C", "Codex F", "less"),
                  H3 = c("New C", "New C with dCor", "less"), H5 = c("F guarded", "Codex F", "less"),
                  H7 = c("F+ depth-weighted", "F+", "more"))
at_error <- function(A, B, target) {
  if (target < min(A) || target > max(A)) return(NA_real_)
  o <- order(A); approx(A[o], B[o], xout = target, ties = mean)$y
}
matched <- rbindlist(lapply(names(QUESTIONS), function(h) {
  q <- QUESTIONS[[h]]
  rbindlist(lapply(names(SUBSETS), function(sn) {
    cur <- fr[subset == sn]
    ra <- range(cur[family == q[1], error]); rb <- range(cur[family == q[2], error])
    lo <- max(ra[1], rb[1]); hi <- min(ra[2], rb[2])
    if (!(hi > lo)) return(data.table(level = NA_real_, hypothesis = h, a = q[1], b = q[2], subset = sn,
                                      folds = 0L, a_wins = 0L, gain_a = NA_real_, gain_b = NA_real_))
    lv <- seq(lo + 0.1 * (hi - lo), hi - 0.1 * (hi - lo), length.out = 6)
    ff <- fr_folds[subset == sn & family %in% q[1:2]]
    g <- ff[, .(level = lv, gain = vapply(lv, function(t) at_error(error, gain, t), 0)), by = .(fold, family)]
    w <- dcast(g, fold + level ~ family, value.var = "gain"); setnames(w, q[1:2], c("ga", "gb"))
    w[, .(hypothesis = h, a = q[1], b = q[2], subset = sn, folds = sum(!is.na(ga) & !is.na(gb)),
          a_wins = sum(ga > gb, na.rm = TRUE), gain_a = round(mean(ga, na.rm = TRUE), 5),
          gain_b = round(mean(gb, na.rm = TRUE), 5)), by = level]
  }))
}))
fwrite(matched, file.path(res, "Task51_Matched.csv"))
hyp <- matched[, .(cells = sum(folds), a_wins = sum(a_wins)), by = .(hypothesis, a, b, subset)]
hyp[, share := r4(a_wins / cells)]
hyp[, holds := mapply(function(h, s) if (QUESTIONS[[h]][3] == "more") s > 0.5 else s < 0.5, hypothesis, share)]
h4 <- dcast(fr[family %in% c("Codex F", "F+")], subset + ell ~ family, value.var = c("error", "gain"))
setnames(h4, c("subset", "ell", "eF", "eP", "gF", "gP"))
h4s <- h4[, .(hypothesis = "H4", a = "F+", b = "Codex F", cells = .N, a_wins = NA_integer_,
              share = round(max(abs(eP - eF), abs(gP - gF)), 6)), by = subset][, holds := share < 0.0005]

FLAGS <- c(not_unique = "flag_not_unique", poor_mapping = "flag_poor_mapping", SINE = "flag_SINE")
rel <- rbindlist(lapply(names(SUBSETS), function(sn) {
  ix <- SUBSETS[[sn]]
  sb <- 1L + floor(20 * (frank(po$spread[ix], ties.method = "first") - 1) / sum(ix))
  rbindlist(lapply(names(FLAGS), function(flg) {
    f <- pairs[[FLAGS[[flg]]]][ix]
    t <- data.table(sb, f, A = po$heldout_error[ix], B = po$heldout_gain[ix])[!is.na(f)][
      , .(A = mean(A), B = mean(B), n = .N), by = .(sb, f)]
    w <- merge(t[f == TRUE], t[f == FALSE], by = "sb", suffixes = c("_f", "_o"))
    data.table(subset = sn, flag = flg, flagged_pairs = sum(f, na.rm = TRUE),
               error_diff = round(sum(w$n_f * (w$A_f - w$A_o)) / sum(w$n_f), 5),
               gain_diff = round(sum(w$n_f * (w$B_f - w$B_o)) / sum(w$n_f), 5),
               bins_worse_error = sum(w$A_f > w$A_o), bins_lower_gain = sum(w$B_f < w$B_o), bins = nrow(w))
  }))
}))
rel[, holds := error_diff > 0 & gain_diff < 0 & bins_worse_error >= 12 & bins_lower_gain >= 12]
fwrite(rel, file.path(res, "Task51_Reliability_SpreadMatched.csv"))
h6 <- rel[, .(hypothesis = "H6", a = "flagged", b = "unflagged", cells = .N, a_wins = sum(holds),
              share = r4(mean(holds)), holds = all(holds)), by = subset]
hyps <- rbind(hyp, h4s, h6, use.names = TRUE); setorder(hyps, hypothesis, subset)
fwrite(hyps, file.path(res, "Task51_Hypotheses.csv"))
cat("\nFrozen hypotheses on the GTEx colon RRBS cohort:\n"); print(hyps)
cat("\nReplicates (holds genome-wide and in both halves):\n"); print(hyps[, .(replicates = all(holds)), by = hypothesis])
cat("\nH6 in detail:\n"); print(rel)

decay <- data.table(bin = cut(gap, c(0, 10, 20, 50, 100, 200, 500, Inf)), gap, sp = spearman_rows(X, Y))[
  !is.na(sp), .(pairs = .N, median_gap = as.numeric(median(gap)), median_spearman = r4(median(sp))), by = bin][order(bin)]
fwrite(decay, file.path(res, "Task51_DistanceDecay.csv"))
cat("\nExploratory: neighbour Spearman by gap:\n"); print(decay)

# ----------------------------------------------------------------
# 5. figure
# ----------------------------------------------------------------
SURF <- "#fcfcfb"; INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"; GRID <- "#e1e0d9"
FCOL <- c(`Codex F` = "#52514e", `F no distance` = "#52514e", `F guarded` = "#eb6834", `New C` = "#1baf7a",
          `New C with dCor` = "#4a3aa7", `F+ depth-weighted` = "#2a78d6", `F+` = "#2a78d6")
FLT <- c(`Codex F` = "solid", `F no distance` = "22", `F guarded` = "solid", `New C` = "solid",
         `New C with dCor` = "solid", `F+ depth-weighted` = "solid", `F+` = "13")
f37 <- fr[subset == "all"][, family := factor(family, levels = names(FCOL))][order(family, -ell)]
p37 <- ggplot(f37, aes(error, gain, colour = family, linetype = family)) +
  geom_path(linewidth = 0.8, lineend = "round") +
  geom_point(data = f37[ell == 200], shape = 21, fill = SURF, stroke = 1, size = 2.4, show.legend = FALSE) +
  scale_colour_manual(values = FCOL) + scale_linetype_manual(values = FLT) +
  labs(title = "The frozen test on a public RRBS cohort: 29 GTEx colon samples (GSE233417)",
       subtitle = paste0("Each curve moves from agreement alone toward agreement times co-methylation as the length ",
                         "scale shrinks. Open dot: 200 bp (the frozen F+).\nTop 10% of ", format(nrow(pairs), big.mark = ","),
                         " truly consecutive CpG pairs (10+ reads in every sample), mean over 10 patient folds."),
       x = "Held-out agreement error |x - y| (lower is better)", y = "Held-out prediction gain (higher is better)",
       caption = "Source: Results/Task51_Frontier.csv (Task 51, frozen plan of Task 47)") +
  theme_minimal(base_size = 11) +
  theme(plot.background = element_rect(fill = SURF, colour = NA),
        panel.grid.major = element_line(colour = GRID, linewidth = 0.3), panel.grid.minor = element_blank(),
        plot.title = element_text(colour = INK, face = "bold", size = 13), plot.title.position = "plot",
        plot.subtitle = element_text(colour = INK2, size = 9.5), plot.caption = element_text(colour = MUTED, hjust = 0),
        plot.caption.position = "plot", legend.position = "top", legend.justification = "left",
        legend.title = element_blank(), legend.key.width = unit(28, "pt"),
        axis.text = element_text(colour = INK2), axis.title = element_text(colour = INK2))
ggsave(file.path(res, "Fig37_GTExColon_Frontier.png"), p37, width = 10, height = 6.6, dpi = 200, bg = SURF)
cat("\ndone\n")
