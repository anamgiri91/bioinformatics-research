# ================================================================
# Task 56 - External cohort for the shared-read noise test (plan
#           phase F): are the 57 white-blood-cell samples of GSE233417
#           independent of the 29 colon donors, and of each other?
# Date: Oct 1, 2026
#
# WHY THIS COHORT (from metadata only; Results/Task56_CohortChoice.md)
#   It is the only public set found with read-level files and at
#   least 30 people of one tissue who are not development donors. The
#   other GTEx tissues share donors with the colon set, and no other
#   read-level series has 30 eligible donors of one tissue. The samples
#   carry no donor IDs, so independence is checked from the data.
#
# THE SCREEN (rules fixed here, before any data were read)
#   Genotype-like CpGs: autosomal CpGs with at least 15 reads in at
#   least 40 of the 57 blood samples, where at least 10% of samples
#   sit in each of three clusters (beta < 0.15; 0.3 to 0.7; > 0.85)
#   and at most 10% sit between them. These behave like genotypes
#   (CpG-SNPs, genotype-driven methylation), which do not change with
#   tissue. If fewer than 200 sites pass, the middle cluster widens to
#   0.2 to 0.8 and the gap limit to 15%.
#   Calls: with at least 15 reads, 0 (beta < 0.15), 1 (0.3 to 0.7) or
#   2 (> 0.85); otherwise no call.
#   Concordance of two samples: the share of sites called in both with
#   the same call (at least 100 such sites).
#   Same person: concordance >= 0.8.
#   Positive controls: 8 esophagus samples from 5 donors who are also
#   in the colon set. All same-donor control pairs must reach 0.8 and
#   every pair of different colon donors must stay below it, or the
#   screen fails and external validation is declared unmet.
#   Exclusions: a blood sample matching a colon donor is dropped; of
#   two matching blood samples, the later GEO accession is dropped. If
#   fewer than 30 blood samples remain, external validation is unmet.
#
# OUTPUTS
#   Results/Task56_IdentityScreen.csv   concordance by pair type
#   Results/Task56_ExternalSamples.csv  blood samples kept, with reasons
#   Results/Fig42_IdentityScreen.png
# ================================================================

suppressPackageStartupMessages({ library(data.table); library(jsonlite); library(ggplot2) })

pick_dir <- function(...) { for (d in c(...)) if (dir.exists(d)) return(d)
  stop("no candidate directory exists") }
repo_dir <- pick_dir("/mmfs1/home/wln26/Experiments.Outlier.July31.2026",
                     path.expand("~/Desktop/bioinformatics-research"))
res <- file.path(repo_dir, "Results"); D0 <- file.path(repo_dir, "Data")
AUTO <- paste0("chr", 1:22); MINR <- 15L
cpg <- readRDS(file.path(D0, "hg19_seq", "cpg_positions_hg19.rds")); NCPG <- sum(lengths(cpg))

wbc <- fromJSON(file.path(D0, "gtex_wbc_rrbs", "wbc_samples.json"))
col <- fromJSON(file.path(D0, "gtex_colon_rrbs", "colon_samples.json"))
eso <- fromJSON(file.path(D0, "gtex_eso_controls", "eso_controls.json"))
smp <- rbind(
  data.table(group = "blood", gsm = wbc$gsm, id = sub(".*\\[(WBC-RRBS-[0-9]+)\\].*", "\\1", wbc$title),
             donor = NA_character_, file = file.path(D0, "gtex_wbc_rrbs", basename(sapply(wbc$files, `[`, 2)))),
  data.table(group = "colon", gsm = col$gsm, id = sub("^GSM[0-9]+_(GTEX-[A-Z0-9]+)-.*", "\\1", basename(sapply(col$files, `[`, 1))),
             donor = sub("^GSM[0-9]+_(GTEX-[A-Z0-9]+)-.*", "\\1", basename(sapply(col$files, `[`, 1))),
             file = file.path(D0, "gtex_colon_rrbs", basename(sapply(col$files, `[`, 2)))),
  data.table(group = "esophagus control", gsm = eso$gsm, id = paste0(eso$donor, " (", eso$gsm, ")"), donor = eso$donor,
             file = file.path(D0, "gtex_eso_controls", basename(sapply(eso$files, `[`, 2)))))
stopifnot(nrow(smp) == 57 + 29 + 8, all(file.exists(smp$file)))

site_counts <- function(f) {                # per-CpG reads and methylated reads, autosomes, global index
  p <- fread(cmd = sprintf("gzip -dc '%s'", f), header = FALSE, sep = "\t", col.names = c("chr", "idx", "pat", "n"),
             colClasses = c("character", "integer", "character", "integer"))[chr %in% AUTO]
  L <- nchar(p$pat); codes <- utf8ToInt(paste0(p$pat, collapse = ""))
  row <- rep.int(seq_len(nrow(p)), L); off <- sequence(L) - 1L; k <- which(codes != 46L)
  data.table(idx = p$idx[row[k]] + off[k], n = p$n[row[k]], m = codes[k] == 67L)[, .(N = sum(n), M = sum(n[m])), keyby = idx]
}

# ----------------------------------------------------------------
# 1. find genotype-like CpGs from the blood samples
# ----------------------------------------------------------------
cnt <- list(cov = integer(NCPG), low = integer(NCPG), mid = integer(NCPG), high = integer(NCPG),
            mid2 = integer(NCPG), high2 = integer(NCPG), low2 = integer(NCPG))
blood <- smp[group == "blood"]
t0 <- Sys.time()
for (s in seq_len(nrow(blood))) {
  sc <- site_counts(blood$file[s])[N >= MINR]
  b <- sc$M / sc$N; i <- sc$idx
  cnt$cov[i] <- cnt$cov[i] + 1L
  cnt$low[i] <- cnt$low[i] + (b < 0.15); cnt$mid[i] <- cnt$mid[i] + (b >= 0.3 & b <= 0.7)
  cnt$high[i] <- cnt$high[i] + (b > 0.85)
  cnt$low2[i] <- cnt$low2[i] + (b < 0.15); cnt$mid2[i] <- cnt$mid2[i] + (b >= 0.2 & b <= 0.8)
  cnt$high2[i] <- cnt$high2[i] + (b > 0.85)
  if (s %% 10 == 0) cat(sprintf("blood %d/%d (%.1f min)\n", s, nrow(blood), as.numeric(Sys.time() - t0, units = "mins")))
}
pick <- function(lo, mi, hi, gap_max) {
  n <- cnt$cov; gap <- n - lo - mi - hi
  which(n >= 40 & lo >= 0.1 * n & mi >= 0.1 * n & hi >= 0.1 * n & gap <= gap_max * n)
}
sites <- pick(cnt$low, cnt$mid, cnt$high, 0.10); rule <- "primary (middle 0.3-0.7, gap <= 10%)"
if (length(sites) < 200) { sites <- pick(cnt$low2, cnt$mid2, cnt$high2, 0.15); rule <- "fallback (middle 0.2-0.8, gap <= 15%)" }
cat("genotype-like CpGs:", length(sites), "using the", rule, "rule\n")
rm(cnt); invisible(gc())

# ----------------------------------------------------------------
# 2. calls for every sample at those sites, and concordance
# ----------------------------------------------------------------
G <- matrix(NA_integer_, length(sites), nrow(smp))
for (s in seq_len(nrow(smp))) {
  sc <- site_counts(smp$file[s])[N >= MINR]
  at <- match(sites, sc$idx); h <- which(!is.na(at)); b <- sc$M[at[h]] / sc$N[at[h]]
  G[h, s] <- fifelse(b < 0.15, 0L, fifelse(b >= 0.3 & b <= 0.7, 1L, fifelse(b > 0.85, 2L, NA_integer_)))
}
pairs <- CJ(a = seq_len(nrow(smp)), b = seq_len(nrow(smp)))[a < b]
pairs[, c("compared", "concordance") := {
  v <- mapply(function(x, y) { k <- !is.na(G[, x]) & !is.na(G[, y]); c(sum(k), mean(G[k, x] == G[k, y])) }, a, b)
  list(v[1, ], v[2, ])
}]
pairs[, `:=`(ga = smp$group[a], gb = smp$group[b], da = smp$donor[a], db = smp$donor[b])]
pairs[, type := fcase(
  ga == "blood" & gb == "blood", "blood vs blood",
  (ga == "blood") != (gb == "blood") & (ga == "colon" | gb == "colon"), "blood vs colon",
  (ga == "blood") != (gb == "blood"), "blood vs esophagus control",
  !is.na(da) & !is.na(db) & da == db, "same GTEx donor (positive control)",
  ga == "colon" & gb == "colon", "colon vs colon (different donors)",
  default = "other different GTEx donors")]
pairs <- pairs[compared >= 100]
SAME <- 0.8
summ <- pairs[, .(pairs = .N, min = round(min(concordance), 3), median = round(median(concordance), 3),
                  max = round(max(concordance), 3), at_or_above_0.8 = sum(concordance >= SAME),
                  median_sites = median(compared)), by = type][order(type)]
fwrite(summ, file.path(res, "Task56_IdentityScreen.csv"))
cat("\nconcordance by pair type:\n"); print(summ)

ctrl_ok <- all(pairs[type == "same GTEx donor (positive control)", concordance >= SAME]) &&
           all(pairs[type == "colon vs colon (different donors)", concordance < SAME]) &&
           pairs[type == "same GTEx donor (positive control)", .N] >= 8
cat("positive controls separate same from different people:", ctrl_ok, "\n")

# ----------------------------------------------------------------
# 3. which blood samples are kept
# ----------------------------------------------------------------
keep <- smp[group == "blood", .(gsm, id, kept = TRUE, reason = "independent")]
bc <- pairs[type == "blood vs colon" & concordance >= SAME]
for (k in seq_len(nrow(bc))) { s <- if (smp$group[bc$a[k]] == "blood") bc$a[k] else bc$b[k]
  keep[gsm == smp$gsm[s], `:=`(kept = FALSE, reason = "matches a colon donor")] }
bb <- pairs[type == "blood vs blood" & concordance >= SAME]
for (k in seq_len(nrow(bb))) { later <- max(smp$gsm[bb$a[k]], smp$gsm[bb$b[k]])
  keep[gsm == later & kept == TRUE, `:=`(kept = FALSE, reason = "duplicate of an earlier blood sample")] }
if (!ctrl_ok) keep[, `:=`(kept = FALSE, reason = "screen failed its positive controls")]
fwrite(keep, file.path(res, "Task56_ExternalSamples.csv"))
cat("\nblood samples kept:", sum(keep$kept), "of", nrow(keep),
    if (sum(keep$kept) >= 30 && ctrl_ok) "-> external cohort accepted\n" else "-> EXTERNAL VALIDATION UNMET\n")

# ----------------------------------------------------------------
# 4. figure
# ----------------------------------------------------------------
SURF <- "#fcfcfb"; INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"; GRID <- "#e1e0d9"
lv <- c("same GTEx donor (positive control)", "colon vs colon (different donors)", "other different GTEx donors",
        "blood vs blood", "blood vs colon", "blood vs esophagus control")
pairs[, type := factor(type, levels = rev(lv[lv %in% unique(type)]))]
p42 <- ggplot(pairs, aes(concordance, type)) +
  geom_vline(xintercept = SAME, colour = "#d03b3b", linewidth = 0.5, linetype = "22") +
  geom_jitter(height = 0.18, width = 0, size = 0.9, alpha = 0.45, colour = "#2a78d6") +
  labs(title = "Are the blood samples different people from the colon donors, and from each other?",
       subtitle = paste0("Share of genotype-like CpGs (", format(length(sites), big.mark = ","),
                         " sites) with the same call. Dashed line: 0.8, the same-person rule."),
       x = "Concordance of genotype-like calls", y = NULL,
       caption = "Source: Results/Task56_IdentityScreen.csv (Task 56)") +
  theme_minimal(base_size = 10.5) +
  theme(plot.background = element_rect(fill = SURF, colour = NA), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = GRID, linewidth = 0.3),
        plot.title = element_text(colour = INK, face = "bold", size = 12), plot.title.position = "plot",
        plot.subtitle = element_text(colour = INK2, size = 9.5), plot.caption = element_text(colour = MUTED, hjust = 0),
        plot.caption.position = "plot", axis.text = element_text(colour = INK2))
ggsave(file.path(res, "Fig42_IdentityScreen.png"), p42, width = 10, height = 4.8, dpi = 200, bg = SURF)
cat("\ndone\n")
