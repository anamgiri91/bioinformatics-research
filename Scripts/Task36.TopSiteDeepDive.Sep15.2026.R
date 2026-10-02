# ================================================================
# Task 36 - Are the "top contrary-flag sites" really outliers?
# Date: Sep 15, 2026
#
# THE QUESTION (supervisor, on pages 5-6 of Sep_3_REPORT.docx)
#
#   The report listed the CpG sites carrying the most contrary `ext`
#   flags in the 100-CpG window (Task 29). Check them more carefully:
#
#   (1) look at the raw beta values - boxplot, the outliers.coef2 /
#       outliers.coef3 columns, the variance, the six-number summary
#       (Min, Q1, Median, Mean, Q3, Max) - and decide whether the
#       flagged values should be called outliers at all;
#   (2) look at CG sites close by (within 200 / 500 bp) and their
#       beta values. cg01731685 (17,084,565) and cg03304299
#       (17,084,852) are 287 bp apart: what do the 53 normals look
#       like at both?
#   (3) any other way of digging into why these sites are flagged.
#
# WHAT IS CHECKED, AND WHY EACH CHECK
#
#   A. Distribution at the site. The source files already carry
#      Min/Q1/Median/Mean/Q3/Max, st.dev, diff.Q3Q1 and the Tukey
#      counts outliers.coef2/3. These are re-derived from the 53
#      betas to confirm what they mean, then split by side (low vs
#      high) - a contrary flag at an L site is a LOW value, so a
#      Tukey count driven by HIGH values says nothing about it.
#
#   B. Where each flagged value sits: rank among the 53, robust z
#      (MAD units), Tukey fence status, and how far past the
#      reconstructed external threshold it lies, in units of the
#      technical noise Task 32 measured for that state.
#
#   C. Replication in an independent cohort: the percentile of each
#      flagged value inside the 32-sample "Dead" cohort of the same
#      tissue. Nobody in that cohort shaped the threshold or the flag.
#
#   D. Neighbours. Every 450K probe within 500 bp (Zhou hg38 manifest,
#      so probes filtered out of the 380k matrix are still listed),
#      its state, median, correlation with the focal site across the
#      53 samples, and whether the flagged samples are also low there.
#      Plus the number of CG dinucleotides in the hg38 reference
#      sequence within 200 / 500 bp, since most CGs are not on the
#      array at all.
#
#   E. Probe design and mapping: Infinium type, colour channel,
#      mapping quality, CpGs inside the 50-mer probe body, and the
#      Zhou et al. (2017) hg38 masks (non-unique mapping, SNPs).
#
#   F. Sample-level explanations. Is a flagged sample simply low
#      everywhere? Each sample's median deviation from the site median
#      over every L-state probe genome-wide, split by Infinium type.
#      And, if N_k and T_k are the same patient (tested, not assumed),
#      what the same patient's tumour shows at the site.
#
# GENOME BUILD
#
#   The matrices' Start column is hg38: the reference base at Start
#   and Start+1 is "CG" for these probes in hg38 and not in hg19
#   (checked against the UCSC API). The Zhou manifest is hg38 too.
#   The older Data/ tracks (UCSC.CpGI.chr22, SNV.all, IMR90) are not
#   used here because their build is not documented.
#
# INPUTS
#   source matrices (Alive + Dead, Normal + Tumour)
#   Results/Task13_{Normal,Tumor}_extref_flags.csv
#   Results/Task32_MeasuredNoise.csv
#   Data/zhou_HM450/HM450.hg38.manifest.tsv.gz
#   Data/zhou_HM450/HM450.hg38.mask.tsv.gz
#   Data/zhou_HM450/HM450.hg38.manifest.gencode.v36.tsv.gz
#     (github.com/zhou-lab/InfiniumAnnotationV1, Anno/HM450/)
#   hg38 sequence from api.genome.ucsc.edu, cached in Data/hg38_seq/
#
# PRIVACY
#   *Detail* outputs and the site summary (it carries Min and Max,
#   i.e. single individuals' values at a named CpG) are private under
#   the same rule as Task 29's detail tables - see .gitignore.
#   Task36_SampleOffset, Task36_PairingCheck and Task36_SequenceCG
#   carry no individual value at an identified site.
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
out_dir  <- file.path(repo_dir, "Results")
data_dir <- file.path(repo_dir, "Data")
zhou_dir <- file.path(data_dir, "zhou_HM450")
seq_dir  <- file.path(data_dir, "hg38_seq"); dir.create(seq_dir, showWarnings = FALSE)

inputs <- list(
  AliveNormal = "Sorted.BRCA.53Alive.Normal.380355cg.75col.May28.2026.txt",
  AliveTumor  = "Sorted.BRCA.53Alive.Tumor.380355cg.75col.May28.2026.txt",
  DeadNormal  = "Sorted.BRCA.32Dead.Normal.380355cg.54col.May28.2026.txt",
  DeadTumor   = "Sorted.BRCA.32Dead.Tumor.380355cg.54col.May28.2026.txt")
inputs <- lapply(inputs, function(f) file.path(src_dir, f))
sample_re <- c(AliveNormal = "^N[0-9]+$", AliveTumor = "^T[0-9]+$",
               DeadNormal  = "^ND[0-9]+$", DeadTumor = "^TD[0-9]+$")
flag_ext <- list(Normal = file.path(out_dir, "Task13_Normal_extref_flags.csv"),
                 Tumor  = file.path(out_dir, "Task13_Tumor_extref_flags.csv"))

# The sites exactly as the Sep 3 report tabulated them (pages 5-6).
TOP <- data.frame(
  tissue = c(rep("Normal", 8), rep("Tumor", 5)),
  state  = c("L","L","LM","LM","H","H","HM","HM",  "L","LM","LM","HM","HM"),
  contrary.dir = c(-1,-1,-1,-1, 1, 1, 1, 1,  -1,-1,-1, 1, 1),
  cgID = c("cg01731685","cg03304299","cg25563983","cg20009101",
           "cg17204569","cg03966099","cg13472645","cg01474260",
           "cg03304299","cg15668074","cg16234298","cg20353484","cg07885834"),
  report.flags = c(4,2,8,4,4,3,5,5,  2,19,6,10,7))
FOCAL <- unique(TOP$cgID)
PAIR  <- c("cg01731685", "cg03304299")
NB_BP <- 500
PAIR_SD <- 0.25
SHIFT_THR <- 0.02

write_out <- function(tbl, f) {
  p <- file.path(out_dir, f); fwrite(tbl, p)
  cat("  written: ", p, "  (", nrow(tbl), " rows)\n", sep = "")
}
r4 <- function(x) round(x, 4)
rowmed <- function(M) apply(M, 1, median, na.rm = TRUE)

# ------------------------------------------------
# 1. Annotation: Zhou hg38 manifest, masks, gene context
# ------------------------------------------------
cat("=========================================================\n")
cat("TASK 36 - TOP CONTRARY-FLAG SITES: ARE THEY OUTLIERS?\n")
cat("=========================================================\n\n")

man <- fread(file.path(zhou_dir, "HM450.hg38.manifest.tsv.gz"),
             select = c("CpG_chrm","CpG_beg","Probe_ID","type","channel",
                        "mapQ_A","AlleleA_ProbeSeq","AlleleB_ProbeSeq"),
             showProgress = FALSE)
msk <- fread(file.path(zhou_dir, "HM450.hg38.mask.tsv.gz"), showProgress = FALSE)
gen <- fread(file.path(zhou_dir, "HM450.hg38.manifest.gencode.v36.tsv.gz"),
             select = c("probeID","genesUniq","distToTSS","CGI","CGIposition"),
             showProgress = FALSE)
man[msk, on = c(Probe_ID = "Probe_ID"), `:=`(maskUniq = i.maskUniq, M_general = i.M_general)]
man[, pos := CpG_beg + 1L]    # 1-based C position, matches the matrices' Start

# CpGs inside the probe body other than the target. Type II probes write
# an underlying CpG as R; type I probes carry CG in the methylated allele.
man[, body.cpg := ifelse(type == "II",
        pmax(nchar(gsub("[^R]", "", AlleleA_ProbeSeq)), 0L),
        pmax(lengths(regmatches(AlleleB_ProbeSeq, gregexpr("CG", AlleleB_ProbeSeq))) - 1L, 0L))]

# ------------------------------------------------
# 2. Noise reference (Task 32): sd bound from adjacent probe pairs
# ------------------------------------------------
noise <- fread(file.path(out_dir, "Task32_MeasuredNoise.csv"))
noise <- noise[pair.set == "<= 100 bp"]
# argument names must not be column names, or data.table scopes them to columns
sigma_of <- function(tis, st) {
  v <- noise$sigma.median[noise$tissue == tis & noise$methy.state == st]
  stopifnot(length(v) == 1); v
}

# ------------------------------------------------
# 3. Load the four cohorts: keep chr22, compute sample-level offsets,
#    keep the SNP-driven probes needed for the N_k <-> T_k pairing test
# ------------------------------------------------
coh <- list(); offsets <- list(); var_mat <- list()

for (nm in names(inputs)) {
  d <- fread(inputs[[nm]], sep = "\t", showProgress = FALSE)
  cols <- grep(sample_re[[nm]], names(d), value = TRUE)
  B <- as.matrix(d[, ..cols]); rownames(B) <- d$Composite.Element.REF
  cat(nm, ": ", nrow(B), " CpGs x ", ncol(B), " samples\n", sep = "")

  # F. sample-wide offset over L-state probes, all and by Infinium type
  typ <- man$type[match(rownames(B), man$Probe_ID)]
  L <- which(d$methy.state == "L")
  devL <- B[L, ] - rowmed(B[L, ])
  off <- data.frame(cohort = nm, sample = cols,
    offset.L.all  = r4(apply(devL, 2, median, na.rm = TRUE)),
    offset.L.typeI  = r4(apply(devL[typ[L] %in% "I", ], 2, median, na.rm = TRUE)),
    offset.L.typeII = r4(apply(devL[typ[L] %in% "II", ], 2, median, na.rm = TRUE)),
    n.L.probes = length(L))
  off$rank.low.L.all <- rank(off$offset.L.all, ties.method = "min")
  offsets[[nm]] <- off; rm(devL)

  if (nm %in% c("AliveNormal", "AliveTumor"))
    var_mat[[nm]] <- B[apply(B, 1, sd, na.rm = TRUE) > PAIR_SD, ]

  k22 <- which(d$Chromosome == "chr22")
  info <- as.data.frame(d[k22, c("Composite.Element.REF","Start","Min","Q1","Median",
            "Mean","Q3","Max","st.dev","diff.Q3Q1","outliers.coef2",
            "outliers.coef3","methy.state"), with = FALSE])
  coh[[nm]] <- list(B = B[k22, , drop = FALSE], info = info)
  rm(d, B); invisible(gc())
}
offsets <- rbindlist(offsets)
offsets[, typeII.minus.typeI := r4(offset.L.typeII - offset.L.typeI)]
offsets[, typeII.shifted := abs(typeII.minus.typeI) > SHIFT_THR]
write_out(offsets, "Task36_SampleOffset.csv")

# A sample whose type II L-probes read high while its type I L-probes do
# not has a chip-level problem, not a biological one: an unmethylated
# CpG is unmethylated whichever chemistry measures it.
cat("\n--- F. type II minus type I offset over L-state probes, genome-wide ---\n")
print(as.data.frame(offsets[, .(median = median(typeII.minus.typeI), MAD = r4(mad(typeII.minus.typeI)),
  shifted = sum(typeII.shifted)), by = cohort]), row.names = FALSE)
cat("samples with |type II - type I| >", SHIFT_THR, ":\n")
print(as.data.frame(offsets[typeII.shifted == TRUE, .(cohort, sample, offset.L.typeI,
  offset.L.typeII, typeII.minus.typeI)]), row.names = FALSE)

# ------------------------------------------------
# 4. Is N_k the same patient as T_k? The README says the 53 are matched
#    pairs, but not that the column index is the pairing. Tested on the
#    probes that vary most in BOTH tissues (sd > 0.25): at that level the
#    variation is genotype-like and survives tumourigenesis. Common-SNP
#    masked probes were tried first, but only 34 survive the 380k filter.
#    At sd > 0.20 (402 probes) tumour-acquired variation swamps it.
# ------------------------------------------------
cat("\n--- pairing check: is N_k the same person as T_k? ---\n")
ids <- intersect(rownames(var_mat$AliveNormal), rownames(var_mat$AliveTumor))
Nm <- var_mat$AliveNormal[ids, ]; Tm <- var_mat$AliveTumor[ids, ]
keep <- rep(TRUE, length(ids))
C <- cor(Nm, Tm, use = "pairwise.complete.obs")
best <- apply(C, 1, which.max)
pairing <- data.frame(normal = rownames(C), best.tumour = colnames(C)[best],
  r.same.index = r4(diag(C)),
  r.best.other = r4(sapply(seq_len(nrow(C)), function(i) max(C[i, -i]))))
pairing$paired <- pairing$best.tumour == sub("^N", "T", pairing$normal)
pairing$rank.of.same.index <- sapply(seq_len(nrow(C)), function(i) rank(-C[i, ])[i])
cat("probes with sd >", PAIR_SD, "in both tissues:", sum(keep), "\n")
cat("N_k whose best-matching tumour is T_k:", sum(pairing$paired), "of", nrow(pairing), "\n")
cat("median r (same index) =", median(pairing$r.same.index),
    " median r (best other) =", median(pairing$r.best.other), "\n")
write_out(pairing, "Task36_PairingCheck.csv")
cat("median rank of T_k among the 53 tumours for N_k:", median(pairing$rank.of.same.index), "\n")
PAIRED <- mean(pairing$paired) > 0.5    # chance is 1 in 53
cat("treating N_k / T_k as the same patient:", PAIRED, "\n")
rm(var_mat, Nm, Tm, C); invisible(gc())

# ------------------------------------------------
# 5. ext flags at chr22, both tissues
# ------------------------------------------------
F <- list()
for (tissue in c("Normal", "Tumor")) {
  fx <- fread(flag_ext[[tissue]], showProgress = FALSE)
  nmc <- paste0("Alive", tissue)
  fx <- fx[match(rownames(coh[[nmc]]$B), fx$cgID)]
  M <- as.matrix(fx[, colnames(coh[[nmc]]$B), with = FALSE])
  rownames(M) <- rownames(coh[[nmc]]$B); F[[tissue]] <- M
  rm(fx); invisible(gc())
}

# threshold bracket from the flags: -1 iff beta < N, +1 iff beta > P
bracket <- function(b, fl) {
  up <- which(fl == 1); dn <- which(fl == -1); ot <- which(fl != -1)
  nlo <- if (length(dn)) max(b[dn], na.rm = TRUE) else NA
  nhi <- if (length(dn) && length(ot)) min(b[ot], na.rm = TRUE) else NA
  ot2 <- which(fl != 1)
  plo <- if (length(up) && length(ot2)) max(b[ot2], na.rm = TRUE) else NA
  phi <- if (length(up)) min(b[up], na.rm = TRUE) else NA
  c(N.lo = nlo, N.hi = nhi, N.mid = mean(c(nlo, nhi)),
    P.lo = plo, P.hi = phi, P.mid = mean(c(plo, phi)))
}

# ------------------------------------------------
# 6. A + E. per-site distribution, Tukey by side, annotation
# ------------------------------------------------
cat("\n--- A. six-number summary, variance and Tukey counts by side ---\n")
tukey_side <- function(b, q1, q3, k) {
  iqr <- q3 - q1
  c(low = sum(b < q1 - k * iqr, na.rm = TRUE), high = sum(b > q3 + k * iqr, na.rm = TRUE))
}
site_rows <- list()
for (nm in names(coh)) {
  tissue <- sub("^(Alive|Dead)", "", nm)
  B <- coh[[nm]]$B; inf <- coh[[nm]]$info
  for (cg in FOCAL) {
    i <- match(cg, rownames(B)); b <- B[i, ]
    q <- quantile(b, c(0, .25, .5, .75, 1), na.rm = TRUE, type = 7)
    t15 <- tukey_side(b, inf$Q1[i], inf$Q3[i], 1.5)
    t2  <- tukey_side(b, inf$Q1[i], inf$Q3[i], 2)
    t3  <- tukey_side(b, inf$Q1[i], inf$Q3[i], 3)
    a <- man[Probe_ID == cg]; g <- gen[probeID == cg]
    row <- data.frame(cohort = nm, cgID = cg, pos = inf$Start[i], state = inf$methy.state[i],
      n = sum(!is.na(b)), Min = r4(inf$Min[i]), Q1 = r4(inf$Q1[i]), Median = r4(inf$Median[i]),
      Mean = r4(inf$Mean[i]), Q3 = r4(inf$Q3[i]), Max = r4(inf$Max[i]),
      sd = r4(inf$st.dev[i]), variance = signif(inf$st.dev[i]^2, 3),
      IQR = r4(inf$diff.Q3Q1[i]), MAD = r4(mad(b, na.rm = TRUE)),
      summary.max.absdiff.vs.file = signif(max(abs(unname(q) - c(inf$Min[i], inf$Q1[i],
                               inf$Median[i], inf$Q3[i], inf$Max[i]))), 2),
      coef2.file = inf$outliers.coef2[i], coef3.file = inf$outliers.coef3[i],
      tukey1.5.low = t15[["low"]], tukey1.5.high = t15[["high"]],
      tukey2.low = t2[["low"]], tukey2.high = t2[["high"]],
      tukey3.low = t3[["low"]], tukey3.high = t3[["high"]],
      ext.neg1 = NA_integer_, ext.pos1 = NA_integer_,
      extN.lo = NA_real_, extN.hi = NA_real_, extP.lo = NA_real_, extP.hi = NA_real_,
      infinium = a$type, channel = a$channel, mapQ = a$mapQ_A, cpg.in.probe.body = a$body.cpg,
      zhou.mask = ifelse(is.na(a$maskUniq) | a$maskUniq == "", "none", a$maskUniq),
      gene = g$genesUniq, dist.to.TSS = g$distToTSS, CGI.position = g$CGIposition)
    if (grepl("^Alive", nm)) {
      fl <- F[[tissue]][cg, ]; br <- bracket(b, fl)
      row$ext.neg1 <- sum(fl == -1, na.rm = TRUE); row$ext.pos1 <- sum(fl == 1, na.rm = TRUE)
      row$extN.lo <- r4(br[["N.lo"]]); row$extN.hi <- r4(br[["N.hi"]])
      row$extP.lo <- r4(br[["P.lo"]]); row$extP.hi <- r4(br[["P.hi"]])
    }
    site_rows[[length(site_rows) + 1]] <- row
  }
}
sites <- rbindlist(site_rows)
# Genome-wide on Alive Normal, coef2 reproduces as the count beyond
# Q1/Q3 -/+ 2 x IQR (both sides) at 99.9995% of CpGs and coef3 at 3 x IQR
# at 100%. Any focal-site disagreement is printed rather than hidden.
sites[, coef2.reproduced := coef2.file == tukey2.low + tukey2.high]
sites[, coef3.reproduced := coef3.file == tukey3.low + tukey3.high]
cat("outliers.coef2 / coef3 = count beyond Q1/Q3 -/+ 2 / 3 x IQR, both sides;",
    "reproduced at", sum(sites$coef2.reproduced), "/", sum(sites$coef3.reproduced),
    "of", nrow(sites), "site x cohort rows\n")
if (!all(sites$coef2.reproduced & sites$coef3.reproduced))
  print(as.data.frame(sites[!(coef2.reproduced & coef3.reproduced),
    .(cohort, cgID, Q1, Q3, IQR, coef2.file, tukey2.low, tukey2.high,
      coef3.file, tukey3.low, tukey3.high)]), row.names = FALSE)
cat("six-number summary: largest |file - R type-7 quantile| per cohort:\n")
print(as.data.frame(sites[, .(max.absdiff = max(summary.max.absdiff.vs.file)), by = cohort]),
      row.names = FALSE)
cat("\n")
print(as.data.frame(sites[cohort %in% c("AliveNormal", "AliveTumor"),
  .(cohort, cgID, state, Min, Q1, Median, Mean, Q3, Max, variance, IQR,
    coef2.file, coef3.file, tukey2.low, tukey2.high, ext.neg1, ext.pos1)]), row.names = FALSE)
cat("\n--- E. probe design and mapping ---\n")
print(unique(as.data.frame(sites[, .(cgID, infinium, channel, mapQ, cpg.in.probe.body,
  zhou.mask, gene, dist.to.TSS, CGI.position)])), row.names = FALSE)
write_out(sites, "Task36_TopSites_Summary.csv")

# ------------------------------------------------
# 7. D. neighbours within 500 bp, and hg38 CG density
# ------------------------------------------------
cat("\n--- D. neighbours within", NB_BP, "bp ---\n")
fetch_seq <- function(pos, flank = 500) {
  f <- file.path(seq_dir, sprintf("chr22_%d_pm%d.txt", pos, flank))
  if (file.exists(f)) return(readLines(f))
  u <- sprintf("https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr22;start=%d;end=%d",
               pos - 1L - flank, pos + flank)
  s <- tryCatch(toupper(jsonlite::fromJSON(u)$dna), error = function(e) NA_character_)
  if (!is.na(s)) writeLines(s, f)
  s
}
man22 <- man[CpG_chrm == "chr22"]
nb_rows <- list(); seq_rows <- list()
for (cg in FOCAL) {
  p0 <- coh$AliveNormal$info$Start[match(cg, rownames(coh$AliveNormal$B))]

  s <- fetch_seq(p0)
  if (!is.na(s)) {
    stopifnot(substr(s, 501, 502) == "CG")          # build check, every site
    off <- gregexpr("CG", s)[[1]] - 501L; off <- off[off != 0]
    seq_rows[[cg]] <- data.frame(cgID = cg, pos = p0,
      CG.within.200bp = sum(abs(off) <= 200), CG.within.500bp = sum(abs(off) <= 500),
      array.probes.within.200bp = man22[abs(pos - p0) <= 200 & Probe_ID != cg, .N],
      array.probes.within.500bp = man22[abs(pos - p0) <= 500 & Probe_ID != cg, .N])
  }

  nbs <- man22[abs(pos - p0) <= NB_BP & Probe_ID != cg][order(pos)]
  for (k in seq_len(nrow(nbs))) {
    nb <- nbs$Probe_ID[k]
    row <- data.frame(focal = cg, focal.pos = p0, neighbour = nb, neighbour.pos = nbs$pos[k],
      distance.bp = nbs$pos[k] - p0, infinium = nbs$type[k],
      zhou.mask = ifelse(is.na(nbs$maskUniq[k]) | nbs$maskUniq[k] == "", "none", nbs$maskUniq[k]),
      in.380k.matrix = nb %in% rownames(coh$AliveNormal$B))
    for (nm in c("AliveNormal", "AliveTumor", "DeadNormal")) {
      B <- coh[[nm]]$B; tag <- c(AliveNormal = "AN", AliveTumor = "AT", DeadNormal = "DN")[[nm]]
      if (row$in.380k.matrix) {
        row[[paste0("state.", tag)]]  <- coh[[nm]]$info$methy.state[match(nb, rownames(B))]
        row[[paste0("median.", tag)]] <- r4(median(B[nb, ], na.rm = TRUE))
        row[[paste0("pearson.", tag)]]  <- r4(cor(B[cg, ], B[nb, ], use = "complete.obs"))
        row[[paste0("spearman.", tag)]] <- r4(cor(B[cg, ], B[nb, ], use = "complete.obs",
                                               method = "spearman"))
      }
    }
    if (row$in.380k.matrix) {
      for (tissue in c("Normal", "Tumor")) {
        fl <- F[[tissue]][nb, ]
        row[[paste0("ext.neg1.", substr(tissue, 1, 1))]] <- sum(fl == -1, na.rm = TRUE)
        row[[paste0("ext.pos1.", substr(tissue, 1, 1))]] <- sum(fl ==  1, na.rm = TRUE)
      }
    }
    nb_rows[[length(nb_rows) + 1]] <- row
  }
}
nbt <- rbindlist(nb_rows, fill = TRUE)
seqt <- rbindlist(seq_rows)
cat("\nCG dinucleotides in hg38 vs probes on the 450K array:\n")
print(as.data.frame(seqt), row.names = FALSE)
cat("\nneighbours of the pair (all 450K probes within 500 bp):\n")
print(as.data.frame(nbt[focal %in% PAIR, .(focal, neighbour, distance.bp, infinium, zhou.mask,
  in.380k.matrix, state.AN, median.AN, pearson.AN, spearman.AN, median.DN, pearson.DN,
  ext.neg1.N, ext.pos1.N)]), row.names = FALSE)
write_out(nbt, "Task36_TopSites_Neighbours_Detail.csv")
write_out(seqt, "Task36_SequenceCG.csv")

# ------------------------------------------------
# 8. B + C + F. every contrary flag at the reported sites
# ------------------------------------------------
cat("\n--- B/C/F. every contrary flag at the reported sites ---\n")
fl_rows <- list()
for (r in seq_len(nrow(TOP))) {
  tissue <- TOP$tissue[r]; cg <- TOP$cgID[r]; dir <- TOP$contrary.dir[r]
  nm <- paste0("Alive", tissue); other <- paste0("Alive", setdiff(c("Normal","Tumor"), tissue))
  B <- coh[[nm]]$B; b <- B[cg, ]; fl <- F[[tissue]][cg, ]
  inf <- coh[[nm]]$info[match(cg, rownames(B)), ]
  med <- median(b, na.rm = TRUE); md <- mad(b, na.rm = TRUE); iqr <- inf$Q3 - inf$Q1
  br <- bracket(b, fl); thr <- if (dir == -1) br[["N.mid"]] else br[["P.mid"]]
  sig <- sigma_of(tissue, TOP$state[r])
  dead <- coh[[paste0("Dead", tissue)]]$B[cg, ]
  nb_ids <- nbt[focal == cg & in.380k.matrix == TRUE, neighbour]
  nbdev <- if (length(nb_ids)) B[nb_ids, , drop = FALSE] - rowmed(B[nb_ids, , drop = FALSE]) else NULL
  off <- offsets[cohort == nm]

  for (j in which(fl == dir)) {
    smp <- colnames(B)[j]
    fl_rows[[length(fl_rows) + 1]] <- data.frame(
      tissue = tissue, state = TOP$state[r], direction = ifelse(dir == 1, "hyper", "hypo"),
      cgID = cg, sample = smp, beta = r4(b[j]), cohort.median = r4(med),
      abs.dbeta = r4(abs(b[j] - med)),
      rank.from.flag.side = if (dir == -1) rank(b)[j] else rank(-b)[j],
      robust.z = r4((b[j] - med) / md),
      beyond.tukey1.5 = if (dir == -1) b[j] < inf$Q1 - 1.5 * iqr else b[j] > inf$Q3 + 1.5 * iqr,
      beyond.tukey3   = if (dir == -1) b[j] < inf$Q1 - 3 * iqr else b[j] > inf$Q3 + 3 * iqr,
      ext.threshold.mid = r4(thr), margin.past.threshold = r4(abs(b[j] - thr)),
      noise.sigma = sig, margin.in.sigma = round(abs(b[j] - thr) / sig, 2),
      dbeta.in.sigma = round(abs(b[j] - med) / sig, 2),
      dead.cohort.pct.at.or.beyond = round(100 * mean(if (dir == -1) dead <= b[j] else dead >= b[j],
                                                      na.rm = TRUE), 1),
      dead.cohort.min = r4(min(dead, na.rm = TRUE)), dead.cohort.max = r4(max(dead, na.rm = TRUE)),
      neighbours.in.matrix = length(nb_ids),
      neighbour.median.dev = if (is.null(nbdev)) NA else r4(median(nbdev[, j], na.rm = TRUE)),
      neighbour.dev.rank = if (is.null(nbdev)) NA else
        rank(apply(nbdev, 2, median, na.rm = TRUE) * (if (dir == -1) 1 else -1),
             ties.method = "min")[j],
      sample.offset.L = off$offset.L.all[j], sample.offset.L.rank.low = off$rank.low.L.all[j],
      sample.typeII.shifted = off$typeII.shifted[j],
      paired.other.tissue.beta = if (PAIRED) r4(coh[[other]]$B[cg, sub("^[NT]", substr(other, 6, 6), smp)]) else NA,
      paired.other.tissue.median = if (PAIRED) r4(median(coh[[other]]$B[cg, ], na.rm = TRUE)) else NA)
  }
}
flags <- rbindlist(fl_rows)
print(as.data.frame(flags[, .(tissue, state, cgID, sample, beta, cohort.median, abs.dbeta,
  rank.from.flag.side, robust.z, beyond.tukey1.5, margin.in.sigma, dbeta.in.sigma,
  dead.cohort.pct.at.or.beyond, neighbour.dev.rank, sample.offset.L.rank.low,
  paired.other.tissue.beta)]), row.names = FALSE)
write_out(flags, "Task36_TopSites_ContraryFlags_Detail.csv")

cat("\nroll-up over all", nrow(flags), "contrary flags at the reported sites:\n")
cat("  beyond Tukey 1.5 x IQR on the flagged side:", sum(flags$beyond.tukey1.5), "\n")
cat("  beyond Tukey 3 x IQR on the flagged side:  ", sum(flags$beyond.tukey3), "\n")
cat("  |robust z| >= 3:                           ", sum(abs(flags$robust.z) >= 3), "\n")
cat("  past the threshold by < 1 noise sigma:     ", sum(flags$margin.in.sigma < 1), "\n")
cat("  |dbeta| < 0.05 / < 0.10:                   ", sum(flags$abs.dbeta < 0.05), "/",
    sum(flags$abs.dbeta < 0.10), "\n")
cat("  flagged value matched or exceeded by >=10% of the Dead cohort:",
    sum(flags$dead.cohort.pct.at.or.beyond >= 10), "\n")
cat("  flagged sample is type-II-shifted:         ", sum(flags$sample.typeII.shifted), "\n")
cat("\nper site: rank range of flagged samples and robust z range\n")
print(as.data.frame(flags[, .(flags = .N, best.rank = min(rank.from.flag.side),
  worst.rank = max(rank.from.flag.side), max.abs.z = max(abs(robust.z)),
  median.abs.dbeta = median(abs.dbeta), median.margin.sigma = median(margin.in.sigma),
  median.dead.pct = median(dead.cohort.pct.at.or.beyond)), by = .(tissue, state, cgID)]),
  row.names = FALSE)

# ------------------------------------------------
# 9. The pair, sample by sample
# ------------------------------------------------
cat("\n--- the pair cg01731685 / cg03304299 (287 bp apart) ---\n")
pair_rows <- list()
for (nm in c("AliveNormal", "AliveTumor", "DeadNormal", "DeadTumor")) {
  B <- coh[[nm]]$B; tissue <- sub("Alive|Dead", "", nm)
  x <- B[PAIR[1], ]; y <- B[PAIR[2], ]
  cat(sprintf("  %-12s Pearson r = %6.3f  Spearman rho = %6.3f  (n = %d)\n", nm,
      cor(x, y, use = "complete.obs"), cor(x, y, use = "complete.obs", method = "spearman"),
      sum(complete.cases(x, y))))
  ok <- !offsets[cohort == nm, typeII.shifted]
  cat(sprintf("  %-12s without type-II-shifted samples: Pearson r = %6.3f  Spearman rho = %6.3f  (n = %d)\n",
      "", cor(x[ok], y[ok], use = "complete.obs"),
      cor(x[ok], y[ok], use = "complete.obs", method = "spearman"), sum(ok)))
  pr <- data.frame(cohort = nm, sample = colnames(B),
    beta.cg01731685 = r4(x), rank.cg01731685 = rank(x, ties.method = "min"),
    beta.cg03304299 = r4(y), rank.cg03304299 = rank(y, ties.method = "min"),
    typeII.shifted = offsets[cohort == nm, typeII.shifted])
  if (grepl("^Alive", nm)) {
    pr$ext.cg01731685 <- F[[tissue]][PAIR[1], ]; pr$ext.cg03304299 <- F[[tissue]][PAIR[2], ]
  }
  pair_rows[[nm]] <- pr
}
pair <- rbindlist(pair_rows, fill = TRUE)
write_out(pair, "Task36_Pair_cg01731685_cg03304299_Detail.csv")
cat("\nAlive Normal samples flagged at either site of the pair:\n")
print(as.data.frame(pair[cohort == "AliveNormal" & (ext.cg01731685 != 0 | ext.cg03304299 != 0)]),
      row.names = FALSE)

# ------------------------------------------------
# 10. Figures (base graphics, as Tasks 25 / 34)
# ------------------------------------------------
png_open <- function(f, w, h) {
  png(file.path(out_dir, f), width = w, height = h, res = 190)
}
COL_CON <- "#E15759"; COL_AGR <- "#4E79A7"; COL_SHIFT <- "#F28E2B"

# Fig13: one panel per reported site; Alive / Dead boxes per tissue
png_open("Fig13_TopSites_RawBeta.png", 3400, 2500)
par(mfrow = c(4, 4), mar = c(3.2, 3.8, 3.4, 0.8), mgp = c(2.2, 0.6, 0), las = 1,
    cex.axis = 0.75, cex.main = 0.85, font.main = 1)
for (cg in FOCAL) {
  vals <- list(`N alive` = coh$AliveNormal$B[cg, ], `N dead` = coh$DeadNormal$B[cg, ],
               `T alive` = coh$AliveTumor$B[cg, ],  `T dead` = coh$DeadTumor$B[cg, ])
  boxplot(vals, range = 1.5, outline = FALSE, col = c("grey92","grey97","grey92","grey97"),
          border = "grey40", ylim = range(unlist(vals), na.rm = TRUE), ylab = "beta")
  for (k in 1:4) {
    v <- vals[[k]]; xj <- k + runif(length(v), -0.22, 0.22)
    col <- rep("grey55", length(v)); pch <- rep(16, length(v)); cexs <- rep(0.55, length(v))
    if (k %in% c(1, 3)) {
      tissue <- if (k == 1) "Normal" else "Tumor"
      st <- sites[cohort == paste0("Alive", tissue) & cgID == cg, state]
      cdir <- if (st %in% c("L", "LM")) -1 else if (st %in% c("H", "HM")) 1 else 0
      fl <- F[[tissue]][cg, ]
      # M and R carry no expected direction, so none of their flags is contrary
      if (cdir != 0) { col[fl %in% cdir] <- COL_CON; col[fl %in% -cdir] <- COL_AGR }
      else col[fl %in% c(-1, 1)] <- "black"
      cexs[fl %in% c(-1, 1)] <- 1
      br <- bracket(v, fl)
      if (is.finite(br[["N.mid"]])) segments(k - .4, br[["N.mid"]], k + .4, lty = 2, col = "grey20")
      if (is.finite(br[["P.mid"]])) segments(k - .4, br[["P.mid"]], k + .4, lty = 3, col = "grey20")
    }
    points(xj, v, pch = pch, col = col, cex = cexs)
  }
  s <- sites[cohort == "AliveNormal" & cgID == cg]
  title(sprintf("%s  N:%s T:%s  %s\n%s", cg, s$state, sites[cohort == "AliveTumor" & cgID == cg, state],
                s$infinium, ifelse(s$zhou.mask == "none", "no Zhou mask", s$zhou.mask)))
}
plot.new()
legend("center", bty = "n", cex = 0.85,
  legend = c("contrary ext flag", "state-consistent ext flag", "ext flag at an M / R site", "not flagged",
             "reconstructed N (hypo) threshold", "reconstructed P (hyper) threshold",
             "box = Q1-median-Q3, whiskers 1.5 x IQR"),
  col = c(COL_CON, COL_AGR, "black", "grey55", "grey20", "grey20", "grey40"),
  pch = c(16, 16, 16, 16, NA, NA, NA), lty = c(NA, NA, NA, NA, 2, 3, 1))
dev.off(); cat("  written: ", file.path(out_dir, "Fig13_TopSites_RawBeta.png"), "\n")

# Fig14: the pair and its neighbourhood
png_open("Fig14_Pair_Neighbourhood.png", 2700, 1150)
par(mfrow = c(1, 2), mar = c(4.2, 4.2, 3.2, 1), mgp = c(2.6, 0.7, 0), las = 1,
    cex.axis = 0.8, cex.main = 0.9, font.main = 1)
p1 <- coh$AliveNormal$info$Start[match(PAIR[1], rownames(coh$AliveNormal$B))]
win <- nbt[focal == PAIR[1] & in.380k.matrix == TRUE, .(neighbour, neighbour.pos)]
win <- rbind(win, data.frame(neighbour = PAIR[1], neighbour.pos = p1))
win <- unique(rbind(win, nbt[focal == PAIR[2] & in.380k.matrix == TRUE, .(neighbour, neighbour.pos)]))
win <- win[order(neighbour.pos)]
Bw <- coh$AliveNormal$B[win$neighbour, , drop = FALSE]
shifted <- offsets[cohort == "AliveNormal" & typeII.shifted == TRUE, sample]
flagged_any <- setdiff(colnames(Bw)[F$Normal[PAIR[1], ] != 0 | F$Normal[PAIR[2], ] != 0], shifted)
plot(NA, xlim = range(win$neighbour.pos), ylim = c(0, 1.3 * max(Bw, na.rm = TRUE)),
     xlab = "chr22 position (hg38)", ylab = "beta (53 Alive Normal)",
     main = "Every array probe within 500 bp of the pair, one line per sample")
cgi <- gen[probeID == PAIR[1], CGI]
if (length(cgi) && !is.na(cgi)) {
  rng <- as.numeric(strsplit(sub(".*:", "", cgi), "-")[[1]])
  rect(rng[1], -1, rng[2], 2, col = "#F2F2F2", border = NA)
  text(mean(pmin(pmax(rng, par("usr")[1]), par("usr")[2])), max(Bw, na.rm = TRUE) * 0.97,
       "hg38 CpG island", cex = 0.7, col = "grey40")
}
for (j in setdiff(colnames(Bw), c(flagged_any, shifted))) lines(win$neighbour.pos, Bw[, j], col = "#00000022")
for (j in shifted) lines(win$neighbour.pos, Bw[, j], col = COL_SHIFT, lwd = 1.4)
for (j in flagged_any) lines(win$neighbour.pos, Bw[, j], col = COL_CON, lwd = 1.4)
typ <- man$type[match(win$neighbour, man$Probe_ID)]
lab <- paste0(win$neighbour, " (", typ, ")")
lab_y <- max(Bw, na.rm = TRUE) * (0.62 - 0.07 * (seq_along(lab) %% 3))
text(win$neighbour.pos, lab_y, lab, cex = 0.5, srt = 90, adj = c(0, 0.5),
     col = ifelse(win$neighbour %in% PAIR, "black", "grey45"),
     font = ifelse(win$neighbour %in% PAIR, 2, 1))
abline(v = win$neighbour.pos[win$neighbour %in% PAIR], col = "grey30", lty = 3)
legend("top", bty = "n", cex = 0.7, lwd = 1.4, col = c(COL_CON, COL_SHIFT, "grey70"),
       legend = c(paste("ext-flagged at the pair:", paste(flagged_any, collapse = ", ")),
                  paste("type-II-shifted arrays:", paste(shifted, collapse = ", ")),
                  "other Alive Normal samples"))
x <- coh$AliveNormal$B[PAIR[1], ]; y <- coh$AliveNormal$B[PAIR[2], ]
xd <- coh$DeadNormal$B[PAIR[1], ]; yd <- coh$DeadNormal$B[PAIR[2], ]
plot(x, y, pch = 16, col = "grey45", xlim = range(c(x, xd)) * c(1, 1.12), ylim = range(c(y, yd)),
     xlab = "cg01731685 beta", ylab = "cg03304299 beta",
     main = sprintf("The pair across samples (Alive Normal r = %.2f)", cor(x, y)))
points(xd, yd, pch = 1, col = "grey60")
f1 <- F$Normal[PAIR[1], ] != 0; f2 <- F$Normal[PAIR[2], ] != 0
sh <- names(x) %in% shifted; shd <- colnames(coh$DeadNormal$B) %in%
  offsets[cohort == "DeadNormal" & typeII.shifted == TRUE, sample]
points(xd[shd], yd[shd], pch = 1, col = COL_SHIFT, cex = 1.1)
pc <- ifelse(sh, COL_SHIFT, COL_CON)
points(x[f1 | f2 | sh], y[f1 | f2 | sh], pch = 16, col = pc[f1 | f2 | sh], cex = 1.2)
text(x[f1 | f2 | sh], y[f1 | f2 | sh], names(x)[f1 | f2 | sh], pos = 4, cex = 0.6, col = pc[f1 | f2 | sh])
b1 <- bracket(x, F$Normal[PAIR[1], ]); b2 <- bracket(y, F$Normal[PAIR[2], ])
abline(v = b1[c("N.mid", "P.mid")], h = b2[c("N.mid", "P.mid")], lty = 2, col = "grey30")
legend("topright", inset = c(0.02, 0), bty = "n", cex = 0.7, pch = c(16, 16, 16, 1, 1, NA),
       lty = c(NA, NA, NA, NA, NA, 2),
       col = c(COL_CON, COL_SHIFT, "grey45", "grey60", COL_SHIFT, "grey30"),
       legend = c("ext-flagged at either site", "type-II-shifted array (Alive)", "Alive Normal",
                  "Dead Normal", "type-II-shifted array (Dead)", "reconstructed ext thresholds"))
dev.off(); cat("  written: ", file.path(out_dir, "Fig14_Pair_Neighbourhood.png"), "\n")

# ------------------------------------------------
# 11. Neighbour summary for every listed site (added 2026-09-15, advisor
#     follow-up: "comparable neighbour summaries for the other sites")
# ------------------------------------------------
cat("\n--- neighbour summary, every listed site ---\n")
nb_sum <- rbindlist(lapply(FOCAL, function(cg) {
  r <- nbt[focal == cg]; m <- r[in.380k.matrix == TRUE]
  s <- seqt[cgID == cg]
  data.table(cgID = cg, zhou.mask = sites[cohort == "AliveNormal" & cgID == cg, zhou.mask],
    CG.within.200bp = s$CG.within.200bp, CG.within.500bp = s$CG.within.500bp,
    array.probes.500bp = nrow(r), in.380k = nrow(m), masked.neighbours = sum(r$zhou.mask != "none"),
    neighbour.states.N = paste(sort(unique(m$state.AN)), collapse = "/"),
    median.r.N = if (nrow(m)) r4(median(m$pearson.AN)) else NA_real_,
    max.abs.r.N = if (nrow(m)) r4(max(abs(m$pearson.AN))) else NA_real_,
    median.r.T = if (nrow(m)) r4(median(m$pearson.AT)) else NA_real_,
    max.abs.r.T = if (nrow(m)) r4(max(abs(m$pearson.AT))) else NA_real_,
    neighbour.ext.flags.N = sum(m$ext.neg1.N + m$ext.pos1.N),
    neighbour.ext.flags.T = sum(m$ext.neg1.T + m$ext.pos1.T))
}))
print(as.data.frame(nb_sum), row.names = FALSE)
write_out(nb_sum, "Task36_NeighbourSummary.csv")

# ------------------------------------------------
# 12. Second close pair: cg20009101 / cg16234298, 74 bp apart. Listed in
#     different tissues (Normal LM and Tumour LM) and their states swap
#     between tissues, so every comparison is made within one tissue.
# ------------------------------------------------
PAIR2 <- c("cg20009101", "cg16234298")
cat("\n--- second close pair", PAIR2[1], "/", PAIR2[2], "(74 bp), per tissue ---\n")
p2_sum <- list(); p2_det <- list()
for (tissue in c("Normal", "Tumor")) {
  nm <- paste0("Alive", tissue); dn <- paste0("Dead", tissue)
  B <- coh[[nm]]$B; x <- B[PAIR2[1], ]; y <- B[PAIR2[2], ]
  fx1 <- F[[tissue]][PAIR2[1], ]; fy1 <- F[[tissue]][PAIR2[2], ]
  sh <- offsets[cohort == nm, typeII.shifted]
  st <- coh[[nm]]$info$methy.state[match(PAIR2, rownames(B))]
  xd <- coh[[dn]]$B[PAIR2[1], ]; yd <- coh[[dn]]$B[PAIR2[2], ]
  rx <- rank(x, ties.method = "min"); ry <- rank(y, ties.method = "min")
  both <- which(fx1 != 0 & fy1 != 0)
  same_dir <- both[fx1[both] == fy1[both]]
  flag_txt <- function(f, v, r) { k <- which(f != 0)
    if (!length(k)) return("none")
    paste(sprintf("%s(%s,b=%.3f,rank %d)", names(v)[k], ifelse(f[k] > 0, "+", "-"), v[k], r[k]), collapse = "; ") }
  other_rank <- function(f, rother, n = length(rother)) { k <- which(f != 0)
    if (!length(k)) return("none")
    paste(sprintf("%s:%d", names(rother)[k], rother[k]), collapse = "; ") }
  s1 <- sites[cohort == nm & cgID == PAIR2[1]]; s2 <- sites[cohort == nm & cgID == PAIR2[2]]
  p2_sum[[tissue]] <- data.table(tissue = tissue,
    state.1 = st[1], state.2 = st[2],
    median.1 = s1$Median, median.2 = s2$Median, IQR.1 = s1$IQR, IQR.2 = s2$IQR,
    variance.1 = s1$variance, variance.2 = s2$variance,
    pearson = r4(cor(x, y)), spearman = r4(cor(x, y, method = "spearman")),
    pearson.no.shifted = r4(cor(x[!sh], y[!sh])), spearman.no.shifted = r4(cor(x[!sh], y[!sh], method = "spearman")),
    pearson.dead = r4(cor(xd, yd)), spearman.dead = r4(cor(xd, yd, method = "spearman")),
    median.abs.diff = r4(median(abs(x - y))),
    ext.flags.1 = sum(fx1 != 0), ext.flags.2 = sum(fy1 != 0),
    flagged.both = length(both), flagged.both.same.direction = length(same_dir),
    flagged.1 = flag_txt(fx1, x, rx), flagged.2 = flag_txt(fy1, y, ry),
    rank.at.2.of.flagged.at.1 = other_rank(fx1, ry), rank.at.1.of.flagged.at.2 = other_rank(fy1, rx),
    shifted.arrays = paste(sprintf("%s(rank %d / %d)", names(x)[sh], rx[sh], ry[sh]), collapse = "; "))
  p2_det[[tissue]] <- data.table(tissue = tissue, sample = names(x),
    beta.cg20009101 = r4(x), rank.cg20009101 = rx, ext.cg20009101 = fx1,
    beta.cg16234298 = r4(y), rank.cg16234298 = ry, ext.cg16234298 = fy1, typeII.shifted = sh)
  cat(sprintf("\n  %s: states %s / %s; median %.3f / %.3f; Pearson %.3f (no shifted %.3f; Dead %.3f); Spearman %.3f\n",
      tissue, st[1], st[2], s1$Median, s2$Median, cor(x, y), cor(x[!sh], y[!sh]), cor(xd, yd),
      cor(x, y, method = "spearman")))
}
p2_sum <- rbindlist(p2_sum)
print(t(as.data.frame(p2_sum)))
write_out(p2_sum, "Task36_Pair2_Summary_Detail.csv")
write_out(rbindlist(p2_det), "Task36_Pair2_cg20009101_cg16234298_Detail.csv")

# Fig19: neighbourhood and scatter for the second pair, one row per tissue
png_open("Fig19_Pair2_Neighbourhood.png", 2700, 2200)
par(mfrow = c(2, 2), mar = c(4.2, 4.2, 3.2, 1), mgp = c(2.6, 0.7, 0), las = 1,
    cex.axis = 0.8, cex.main = 0.9, font.main = 1)
p1 <- coh$AliveNormal$info$Start[match(PAIR2, rownames(coh$AliveNormal$B))]
win2 <- unique(rbind(nbt[focal %in% PAIR2 & in.380k.matrix == TRUE, .(neighbour, neighbour.pos)],
                     data.table(neighbour = PAIR2, neighbour.pos = p1)))[order(neighbour.pos)]
for (tissue in c("Normal", "Tumor")) {
  nm <- paste0("Alive", tissue)
  Bw <- coh[[nm]]$B[win2$neighbour, , drop = FALSE]
  shs <- offsets[cohort == nm & typeII.shifted == TRUE, sample]
  fl <- colnames(Bw)[F[[tissue]][PAIR2[1], ] != 0 | F[[tissue]][PAIR2[2], ] != 0]
  fl <- setdiff(fl, shs)
  plot(NA, xlim = range(win2$neighbour.pos), ylim = c(0, 1.25 * max(Bw)),
       xlab = "chr22 position (hg38)", ylab = sprintf("beta (53 Alive %s)", tolower(tissue)),
       main = sprintf("%s: array probes within 500 bp of the pair", tissue))
  for (j in setdiff(colnames(Bw), c(fl, shs))) lines(win2$neighbour.pos, Bw[, j], col = "#00000022")
  for (j in shs) lines(win2$neighbour.pos, Bw[, j], col = COL_SHIFT, lwd = 1.3)
  for (j in fl) lines(win2$neighbour.pos, Bw[, j], col = COL_CON, lwd = 1.3)
  typ2 <- man$type[match(win2$neighbour, man$Probe_ID)]
  msk2 <- ifelse(win2$neighbour %in% nbt[zhou.mask != "none", neighbour] |
                 win2$neighbour %in% sites[zhou.mask != "none", cgID], "*", "")
  text(win2$neighbour.pos, 1.02 * max(Bw) * (1 + 0.07 * (seq_len(nrow(win2)) %% 2)),
       paste0(win2$neighbour, " (", typ2, ")", msk2), cex = 0.5, srt = 90, adj = c(0, 0.5),
       font = ifelse(win2$neighbour %in% PAIR2, 2, 1))
  abline(v = p1, lty = 3, col = "grey30")
  legend("topleft", bty = "n", cex = 0.65, lwd = 1.3, col = c(COL_CON, COL_SHIFT),
         legend = c(paste("ext-flagged at either site:", paste(fl, collapse = ", ")),
                    paste("type-II-shifted arrays:", paste(shs, collapse = ", "))))
  mtext("* = Zhou-masked probe", side = 1, line = 3, adj = 1, cex = 0.55)

  x <- coh[[nm]]$B[PAIR2[1], ]; y <- coh[[nm]]$B[PAIR2[2], ]
  xd <- coh[[paste0("Dead", tissue)]]$B[PAIR2[1], ]; yd <- coh[[paste0("Dead", tissue)]]$B[PAIR2[2], ]
  f1 <- F[[tissue]][PAIR2[1], ] != 0; f2 <- F[[tissue]][PAIR2[2], ] != 0
  shm <- names(x) %in% shs
  plot(x, y, pch = 16, col = "grey45", xlim = range(c(x, xd)) * c(1, 1.1), ylim = range(c(y, yd)),
       xlab = "cg20009101 beta", ylab = "cg16234298 beta (M_nonuniq, mapQ 4)",
       main = sprintf("%s: the pair across samples (r = %.2f)", tissue, cor(x, y)))
  points(xd, yd, pch = 1, col = "grey65")
  cl <- ifelse(shm, COL_SHIFT, ifelse(f1 & f2, "#7B3294", ifelse(f1, COL_CON, COL_AGR)))
  k <- f1 | f2 | shm
  points(x[k], y[k], pch = 16, col = cl[k], cex = 1.2)
  text(x[k], y[k], names(x)[k], pos = 4, cex = 0.55, col = cl[k])
  b1 <- bracket(x, F[[tissue]][PAIR2[1], ]); b2 <- bracket(y, F[[tissue]][PAIR2[2], ])
  abline(v = b1[c("N.mid", "P.mid")], h = b2[c("N.mid", "P.mid")], lty = 2, col = "grey30")
  legend("bottomright", bty = "n", cex = 0.65, pch = c(16, 16, 16, 16, 16, 1, NA), lty = c(rep(NA, 6), 2),
         col = c(COL_CON, COL_AGR, "#7B3294", COL_SHIFT, "grey45", "grey65", "grey30"),
         legend = c("ext flag at cg20009101 only", "ext flag at cg16234298 only", "flagged at both",
                    "type-II-shifted array", "Alive, not flagged", "Dead cohort", "reconstructed ext cutoffs"))
}
dev.off(); cat("  written: ", file.path(out_dir, "Fig19_Pair2_Neighbourhood.png"), "\n")

cat("\n=========================================================\n")
cat("TASK 36 COMPLETE\n")
cat("=========================================================\n")
