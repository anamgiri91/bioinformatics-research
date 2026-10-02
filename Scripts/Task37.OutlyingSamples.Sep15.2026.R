# ================================================================
# Task 37 - Do the outlier flags identify outlying samples?
# Date: Sep 15, 2026
#
# THE QUESTION (supervisor, 2026-09-15, item 2)
#
#   "Do the outlier flags help to identify some outlying samples
#    (e.g., N14-N17)? Outlying sample checking based on the current
#    analysis."
#
# WHY IT IS ANSWERABLE NOW
#
#   Task 33 showed N15, N17, N14 and N48 lead genome-wide ext flag
#   burden, and that the ranking survives a 0.10 magnitude floor.
#   plan.md item 4 left open what drives it. Task 36 then found that
#   N14, N15 and N17 read about +0.08 higher on type II probes only.
#   So there are two things to separate:
#
#     * does flag burden pick out samples that are ALSO outlying by
#       measures that never use a flag or a threshold?
#     * when it does, is the outlyingness biology or the array?
#
# WHAT IS MEASURED, PER SAMPLE, GENOME-WIDE (Alive Normal and Tumour)
#
#   A. Flag burden. ext and self flags from Task 13, total and split
#      by Infinium type I / type II, plus the type II share of each
#      sample's ext flags against the type II share of all probes.
#      The type-I-only ext burden is the key number: a sample whose
#      burden comes from a type II artefact drops back into the pack
#      once type II probes are left out.
#
#   B. Flag-free outlyingness, four ways, each turned into a robust z
#      across the 53 samples (median / MAD):
#        1. type II minus type I offset over L-state probes (Task 36),
#           and the same over H-state probes;
#        2. global distance: median |beta - site median| over all CpGs;
#        3. 1 - Pearson r with the cohort's median profile;
#        4. PCA on the 20,000 most variable CpGs: the largest robust z
#           over PC1-PC5.
#      A sample is called outlying on a measure at |robust z| > 3.
#
#   C. Agreement: Spearman between ext burden and each measure; which
#      samples are outlying on each; where each ranks on burden.
#
#   D. The same four flag-free measures recomputed on type I probes
#      only, to see which samples stay outlying once the type II
#      artefact cannot contribute.
#
# INPUTS
#   source matrices, Alive Normal / Tumour
#   Results/Task13_{Normal,Tumor}_{extref,selfref}_flags.csv
#   Data/zhou_HM450/HM450.hg38.manifest.tsv.gz (Infinium type)
#
# OUTPUT
#   Results/Task37_SampleOutlyingness.csv  - one row per sample: counts
#     and summary scores only, no value at an identified CpG (public,
#     same class as Task 23 / Task 33 sample summaries)
#   Results/Fig15_FlagBurdenVsOutlyingness.png
# ================================================================

options(stringsAsFactors = FALSE)
set.seed(20260915)
suppressPackageStartupMessages(library(data.table))

pick_dir <- function(...) { for (d in c(...)) if (dir.exists(d)) return(d)
  stop("no candidate directory exists") }
src_dir <- pick_dir(
  "/home/s_s355/research3.data/TCGA/filter.BRCA.UCEC.Aug17.2022/sorted.TCGA.380355cg.files",
  path.expand("~/Downloads/tcga_brca_380355cg"))
repo_dir <- pick_dir("/mmfs1/home/wln26/Experiments.Outlier.July31.2026",
                     path.expand("~/Desktop/bioinformatics-research"))
out_dir <- file.path(repo_dir, "Results")

inputs <- list(
  Normal = file.path(src_dir, "Sorted.BRCA.53Alive.Normal.380355cg.75col.May28.2026.txt"),
  Tumor  = file.path(src_dir, "Sorted.BRCA.53Alive.Tumor.380355cg.75col.May28.2026.txt"))
flag_file <- function(tissue, ref) file.path(out_dir, sprintf("Task13_%s_%s_flags.csv", tissue, ref))

Z_OUT  <- 3        # robust z beyond which a sample is called outlying
N_PCA  <- 20000    # most variable CpGs used for PCA
N_PC   <- 5
r4 <- function(x) round(x, 4)
rz <- function(x) (x - median(x, na.rm = TRUE)) / mad(x, na.rm = TRUE)
rowmed <- function(M) apply(M, 1, median, na.rm = TRUE)

man <- fread(file.path(repo_dir, "Data/zhou_HM450/HM450.hg38.manifest.tsv.gz"),
             select = c("Probe_ID", "type"), showProgress = FALSE)

cat("=========================================================\n")
cat("TASK 37 - DO THE FLAGS IDENTIFY OUTLYING SAMPLES?\n")
cat("=========================================================\n")

# flag-free outlyingness on a subset of rows
outlyingness <- function(B, rows, tag) {
  Bs <- B[rows, , drop = FALSE]
  med <- rowmed(Bs)
  dist <- apply(abs(Bs - med), 2, median, na.rm = TRUE)
  ok <- is.finite(med)
  r_med <- apply(Bs[ok, ], 2, function(b) cor(b, med[ok], use = "complete.obs"))
  v <- apply(Bs, 1, var, na.rm = TRUE)
  top <- order(v, decreasing = TRUE)[seq_len(min(N_PCA, length(v)))]
  X <- t(Bs[top, ]); X[is.na(X)] <- 0
  pc <- prcomp(X, center = TRUE, scale. = FALSE)
  S <- pc$x[, seq_len(N_PC), drop = FALSE]
  pcz <- apply(S, 2, rz)
  out <- data.frame(dist = dist, one.minus.r = 1 - r_med,
    pca.max.abs.z = apply(abs(pcz), 1, max),
    pca.worst.pc = apply(abs(pcz), 1, which.max))
  names(out) <- paste0(names(out), tag)
  attr(out, "pc.var") <- round(100 * pc$sdev[seq_len(N_PC)]^2 / sum(pc$sdev^2), 1)
  out
}

all_rows <- list(); fig_data <- list()

for (tissue in c("Normal", "Tumor")) {
  cat("\n#########################################################\n")
  cat("### ", tissue, "\n", sep = "")
  cat("#########################################################\n")

  d <- fread(inputs[[tissue]], sep = "\t", showProgress = FALSE)
  cols <- grep(if (tissue == "Normal") "^N[0-9]+$" else "^T[0-9]+$", names(d), value = TRUE)
  B <- as.matrix(d[, ..cols]); rownames(B) <- d$Composite.Element.REF
  st <- d$methy.state; ids <- d$Composite.Element.REF; rm(d); invisible(gc())
  typ <- man$type[match(ids, man$Probe_ID)]
  cat("CpGs:", nrow(B), " type I:", sum(typ %in% "I"), " type II:", sum(typ %in% "II"),
      " no manifest type:", sum(is.na(typ)), "\n")
  share_II_probes <- mean(typ[!is.na(typ)] == "II")

  # ---- B1. type II - type I offsets in L and H states ----
  offset <- function(state) {
    k <- which(st == state); dev <- B[k, ] - rowmed(B[k, ])
    list(I  = apply(dev[typ[k] %in% "I", ], 2, median, na.rm = TRUE),
         II = apply(dev[typ[k] %in% "II", ], 2, median, na.rm = TRUE))
  }
  oL <- offset("L"); oH <- offset("H")

  # ---- B2-4. flag-free outlyingness, all probes and type I only ----
  o_all <- outlyingness(B, seq_len(nrow(B)), "")
  o_I   <- outlyingness(B, which(typ %in% "I"), ".typeI")
  cat("PCA variance explained, PC1-PC5 (all probes):", attr(o_all, "pc.var"), "%\n")
  rm(B); invisible(gc())

  # ---- A. flag burden, total and by type ----
  burden <- list()
  for (ref in c("extref", "selfref")) {
    fx <- fread(flag_file(tissue, ref), showProgress = FALSE)
    stopifnot(identical(fx$cgID, ids))
    F <- as.matrix(fx[, ..cols]); rm(fx); invisible(gc())
    nz <- F != 0 & !is.na(F)
    burden[[ref]] <- list(all = colSums(nz),
      I = colSums(nz[typ %in% "I", ]), II = colSums(nz[typ %in% "II", ]),
      hyper = colSums(F == 1, na.rm = TRUE), hypo = colSums(F == -1, na.rm = TRUE))
    rm(F, nz); invisible(gc())
  }
  e <- burden$extref; s <- burden$selfref

  res <- data.frame(tissue = tissue, sample = cols,
    ext.flags = e$all, ext.rank = rank(-e$all, ties.method = "min"),
    ext.hyper = e$hyper, ext.hypo = e$hypo,
    ext.flags.typeI = e$I, ext.rank.typeI = rank(-e$I, ties.method = "min"),
    ext.flags.typeII = e$II, ext.rank.typeII = rank(-e$II, ties.method = "min"),
    ext.typeII.share = r4(e$II / (e$I + e$II)),
    self.flags = s$all, self.rank = rank(-s$all, ties.method = "min"),
    self.flags.typeI = s$I, self.rank.typeI = rank(-s$I, ties.method = "min"),
    offset.L.typeI = r4(oL$I), offset.L.typeII = r4(oL$II),
    L.typeII.minus.typeI = r4(oL$II - oL$I),
    offset.H.typeI = r4(oH$I), offset.H.typeII = r4(oH$II),
    H.typeII.minus.typeI = r4(oH$II - oH$I),
    global.dist = r4(o_all$dist), one.minus.r = r4(o_all$one.minus.r),
    pca.max.abs.z = round(o_all$pca.max.abs.z, 2), pca.worst.pc = o_all$pca.worst.pc,
    global.dist.typeI = r4(o_I$dist.typeI), one.minus.r.typeI = r4(o_I$one.minus.r.typeI),
    pca.max.abs.z.typeI = round(o_I$pca.max.abs.z.typeI, 2))

  # robust z of each flag-free measure; PCA is already a robust z
  res$z.L.shift    <- round(rz(res$L.typeII.minus.typeI), 1)
  res$z.H.shift    <- round(rz(res$H.typeII.minus.typeI), 1)
  res$z.dist       <- round(rz(res$global.dist), 1)
  res$z.one.minus.r <- round(rz(res$one.minus.r), 1)
  res$z.dist.typeI <- round(rz(res$global.dist.typeI), 1)
  res$z.one.minus.r.typeI <- round(rz(res$one.minus.r.typeI), 1)
  res$outlying.on <- apply(cbind(
      L.shift = abs(res$z.L.shift) > Z_OUT, H.shift = abs(res$z.H.shift) > Z_OUT,
      dist = res$z.dist > Z_OUT, corr = res$z.one.minus.r > Z_OUT,
      pca = res$pca.max.abs.z > Z_OUT), 1,
    function(v) paste(c("L.shift","H.shift","dist","corr","pca")[v], collapse = ";"))
  res$outlying.on.typeI <- apply(cbind(
      dist = res$z.dist.typeI > Z_OUT, corr = res$z.one.minus.r.typeI > Z_OUT,
      pca = res$pca.max.abs.z.typeI > Z_OUT), 1,
    function(v) paste(c("dist","corr","pca")[v], collapse = ";"))
  res$typeII.shifted <- abs(res$L.typeII.minus.typeI) > 0.02

  # ---- C. print ----
  cat("\ntype II share of probes:", round(100 * share_II_probes, 1), "%\n")
  cat("\nsamples ranked by ext flag burden (top 10) with every flag-free measure:\n")
  show <- res[order(res$ext.rank), ][1:10, c("sample","ext.flags","ext.rank","ext.rank.typeI",
    "ext.typeII.share","self.rank","L.typeII.minus.typeI","H.typeII.minus.typeI",
    "z.dist","z.one.minus.r","pca.max.abs.z","outlying.on","outlying.on.typeI")]
  print(show, row.names = FALSE)

  cat("\nevery sample outlying on at least one flag-free measure:\n")
  print(res[res$outlying.on != "" | res$outlying.on.typeI != "",
    c("sample","ext.rank","ext.rank.typeI","self.rank","z.L.shift","z.H.shift","z.dist",
      "z.one.minus.r","pca.max.abs.z","pca.worst.pc","outlying.on","outlying.on.typeI")],
    row.names = FALSE)

  cat("\nnamed samples N/T 14-17:\n")
  print(res[res$sample %in% paste0(substr(tissue, 1, 1), 14:17),
    c("sample","ext.flags","ext.rank","ext.rank.typeI","ext.typeII.share","self.rank",
      "L.typeII.minus.typeI","z.dist","z.one.minus.r","pca.max.abs.z","outlying.on")],
    row.names = FALSE)

  cat("\nSpearman correlation of ext flag burden with each flag-free measure:\n")
  sp <- function(a, b) round(cor(a, b, method = "spearman"), 3)
  agree <- data.frame(measure = c("L type II - type I", "global distance", "1 - r with median",
                                  "PCA max |z|", "global distance (type I)",
                                  "1 - r with median (type I)", "self flag burden"),
    all.probe.ext = c(sp(res$ext.flags, res$L.typeII.minus.typeI), sp(res$ext.flags, res$global.dist),
      sp(res$ext.flags, res$one.minus.r), sp(res$ext.flags, res$pca.max.abs.z),
      sp(res$ext.flags, res$global.dist.typeI), sp(res$ext.flags, res$one.minus.r.typeI),
      sp(res$ext.flags, res$self.flags)),
    typeI.only.ext = c(sp(res$ext.flags.typeI, res$L.typeII.minus.typeI),
      sp(res$ext.flags.typeI, res$global.dist), sp(res$ext.flags.typeI, res$one.minus.r),
      sp(res$ext.flags.typeI, res$pca.max.abs.z), sp(res$ext.flags.typeI, res$global.dist.typeI),
      sp(res$ext.flags.typeI, res$one.minus.r.typeI), sp(res$ext.flags.typeI, res$self.flags.typeI)))
  print(agree, row.names = FALSE)

  shifted <- res$sample[res$typeII.shifted]
  cat("\ntype-II-shifted samples:", paste(shifted, collapse = ", "), "\n")
  cat("  their ext ranks (all probes):   ", paste(res$ext.rank[res$typeII.shifted], collapse = ", "), "\n")
  cat("  their ext ranks (type I only):  ", paste(res$ext.rank.typeI[res$typeII.shifted], collapse = ", "), "\n")
  cat("  type II share of their ext flags:", paste(res$ext.typeII.share[res$typeII.shifted], collapse = ", "),
      " (cohort median", median(res$ext.typeII.share), ")\n")
  top5 <- res$sample[order(res$ext.rank)][1:5]; top5I <- res$sample[order(res$ext.rank.typeI)][1:5]
  cat("  top 5 by ext burden, all probes:", paste(top5, collapse = ", "), "\n")
  cat("  top 5 by ext burden, type I only:", paste(top5I, collapse = ", "), "\n")

  all_rows[[tissue]] <- res
}

out <- rbindlist(all_rows)
fwrite(out, file.path(out_dir, "Task37_SampleOutlyingness.csv"))
cat("\n  written: ", file.path(out_dir, "Task37_SampleOutlyingness.csv"), " (", nrow(out), " rows)\n", sep = "")

# ------------------------------------------------
# Figure 15: burden against the type II shift, all probes vs type I only
# ------------------------------------------------
png(file.path(out_dir, "Fig15_FlagBurdenVsOutlyingness.png"), width = 2700, height = 2100, res = 190)
par(mfrow = c(2, 3), mar = c(4.4, 4.8, 3.2, 1), mgp = c(2.8, 0.7, 0), las = 1,
    cex.axis = 0.8, cex.main = 0.9, font.main = 1)
COL_SHIFT <- "#F28E2B"
for (tissue in c("Normal", "Tumor")) {
  r <- all_rows[[tissue]]
  col <- ifelse(r$typeII.shifted, COL_SHIFT, "grey40")
  lab <- r$typeII.shifted | r$ext.rank <= 6 | r$ext.rank.typeI <= 4 | r$outlying.on != ""

  xr <- range(r$L.typeII.minus.typeI); xr[2] <- xr[2] + 0.15 * diff(xr)
  plot(r$L.typeII.minus.typeI, r$ext.flags / 1000, pch = 16, col = col, xlim = xr,
       xlab = "type II minus type I offset, L-state probes", ylab = "ext flags (thousands)",
       main = sprintf("%s: ext burden vs type II shift", tissue))
  text(r$L.typeII.minus.typeI[lab], r$ext.flags[lab] / 1000, r$sample[lab], pos = 4, cex = 0.6)

  plot(r$ext.rank, r$ext.rank.typeI, pch = 16, col = col, xlim = c(1, 53), ylim = c(53, 1),
       xlab = "rank by ext burden, all probes (1 = most)",
       ylab = "rank by ext burden, type I probes only",
       main = sprintf("%s: does the ranking survive dropping type II?", tissue))
  abline(0, 1, lty = 3, col = "grey60")
  text(r$ext.rank[lab], r$ext.rank.typeI[lab], r$sample[lab], pos = 4, cex = 0.6)

  xr <- range(c(r$pca.max.abs.z, Z_OUT)); xr[2] <- xr[2] + 0.12 * diff(xr)
  plot(r$pca.max.abs.z, r$ext.flags / 1000, pch = 16, col = col, xlim = xr,
       xlab = "largest robust z on PC1-PC5 (flag-free)", ylab = "ext flags (thousands)",
       main = sprintf("%s: ext burden vs PCA outlyingness", tissue))
  abline(v = Z_OUT, lty = 2, col = "grey50")
  text(r$pca.max.abs.z[lab], r$ext.flags[lab] / 1000, r$sample[lab], pos = 4, cex = 0.6)
}
legend("bottomright", bty = "n", cex = 0.75, pch = 16, col = c(COL_SHIFT, "grey40"),
       legend = c("type-II-shifted array", "other samples"))
dev.off()
cat("  written: ", file.path(out_dir, "Fig15_FlagBurdenVsOutlyingness.png"), "\n")

cat("\n=========================================================\n")
cat("TASK 37 COMPLETE\n")
cat("=========================================================\n")
