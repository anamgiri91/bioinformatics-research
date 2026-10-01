# ================================================================
# Task 50 - Build the public RRBS cohort for sections 5 and 6 of the
#           frozen plan: 29 GTEx colon samples (GSE233417)
# Date: Sep 29, 2026
#
# The cohort was chosen by the frozen rules from GEO metadata only
# (Results/Task47_DataChoice_Section5.md, commit 4481952).
#
# WHAT THIS DOES
#   The GEO per-CpG files (.beta.bed.gz) hold beta values but no read
#   counts. The read-level files (.pat.gz, wgbstools format: chr, first
#   CpG index, one letter per CpG, read count) record every fragment.
#   So, per sample:
#   1. per-CpG methylated and total read counts are rebuilt from the
#      .pat file (C = methylated, T = unmethylated, . = no call);
#   2. each CpG index is given its hg19 position from the genome: the
#      wgbstools index counts CpGs in genome order (chr1, chr2, ...), so
#      on chromosome k, index = (CpGs on chr1..k-1) + rank of the CpG on
#      chr k. CpG positions come from the UCSC hg19 sequence. The GEO
#      .beta.bed file is a separate product (about 2% more CpGs, some
#      betas differ, so it cannot be lined up row by row); it is used
#      only to check the map: most CpGs in both must share the beta;
#   3. for every pair of adjacent CpGs (index i and i + 1) the reads that
#      call both sites give the read-level linkage of section 6: the
#      share with the same state at both, minus p_i p_j + (1 - p_i)(1 - p_j)
#      from the same reads. Kept when at least 10 reads call both.
#   The .beta.bed "start" column is the 1-based C position (checked
#   against the hg19 sequence: chr1 10489, 10493, 10497 are CpG Cs).
#
# OUTPUT (Data/gtex_colon_rrbs/, gitignored with Data/)
#   cache_cohort.rds   sites with >= 10 reads in all 29 samples
#                      (autosomes): positions, methylated and total count
#                      matrices; the index-to-position map; per-pair
#                      read-level linkage averaged over samples
#   Results/Task50_BuildSummary.csv   per-sample checks (no values)
# ================================================================

suppressPackageStartupMessages({ library(data.table); library(jsonlite) })

pick_dir <- function(...) { for (d in c(...)) if (dir.exists(d)) return(d)
  stop("no candidate directory exists") }
repo_dir <- pick_dir("/mmfs1/home/wln26/Experiments.Outlier.July31.2026",
                     path.expand("~/Desktop/bioinformatics-research"))
dat <- file.path(repo_dir, "Data", "gtex_colon_rrbs"); res <- file.path(repo_dir, "Results")
AUTO <- paste0("chr", 1:22); MIN_COV <- 10L; MIN_SPAN <- 10L

seq_dir <- file.path(repo_dir, "Data", "hg19_seq"); cpg_file <- file.path(seq_dir, "cpg_positions_hg19.rds")
if (!file.exists(cpg_file)) {
  cpg <- lapply(1:22, function(k) {
    x <- Biostrings::readDNAStringSet(file.path(seq_dir, sprintf("chr%d.fa.gz", k)))[[1]]
    as.integer(BiocGenerics::start(Biostrings::matchPattern("CG", x)))  # 1-based C positions
  })
  saveRDS(cpg, cpg_file)
}
cpg <- readRDS(cpg_file)
I0 <- c(0L, cumsum(lengths(cpg))[-22])       # index offset of each autosome
allpos <- unlist(cpg, use.names = FALSE)      # allpos[global index] = position (autosomes in order)
cat("hg19 autosomal CpGs:", sum(lengths(cpg)), "\n")
samples <- fromJSON(file.path(dat, "colon_samples.json"))
stopifnot(nrow(samples) == 29)
sid <- sub("^GSM[0-9]+_(GTEX-[A-Z0-9]+)-.*", "\\1", basename(sapply(samples$files, `[`, 1)))
stopifnot(!anyDuplicated(sid))

one_sample <- function(s) {
  fp <- file.path(dat, basename(samples$files[[s]][2])); fb <- file.path(dat, basename(samples$files[[s]][1]))
  p <- fread(cmd = sprintf("gzip -dc '%s'", fp), header = FALSE, sep = "\t",
             col.names = c("chr", "idx", "pat", "n"), colClasses = c("character", "integer", "character", "integer"))
  p <- p[chr %in% AUTO]
  chrn <- match(p$chr, AUTO); L <- nchar(p$pat)
  codes <- utf8ToInt(paste0(p$pat, collapse = ""))
  stopifnot(all(codes %in% c(46L, 67L, 84L)))           # only ".", "C", "T"
  line <- rep.int(seq_len(nrow(p)), L); off <- sequence(L) - 1L
  call <- codes != 46L; meth <- codes == 67L
  k <- which(call)
  cg <- data.table(chrn = chrn[line[k]], idx = p$idx[line[k]] + off[k], n = p$n[line[k]], m = meth[k])
  cnt <- cg[, .(meth = sum(n[m]), total = sum(n)), keyby = .(chrn, idx)]
  rm(cg); invisible(gc())
  # adjacent CpGs called on the same read
  K <- length(codes)
  k <- which(c(line[-1] == line[-K], FALSE) & call & c(call[-1], FALSE))
  pr <- data.table(chrn = chrn[line[k]], idx = p$idx[line[k]] + off[k], n = p$n[line[k]], mi = meth[k], mj = meth[k + 1L])
  span <- pr[, .(N = sum(n), both = sum(n[mi & mj]), ni = sum(n[mi]), nj = sum(n[mj])), keyby = .(chrn, idx)][N >= MIN_SPAN]
  rm(pr, codes, line, off, call, meth, p); invisible(gc())
  span[, `:=`(pi = ni / N, pj = nj / N)]
  span[, linkage := (both + (N - ni - nj + both)) / N - (pi * pj + (1 - pi) * (1 - pj))]
  # positions from the genome, checked against the GEO beta file
  rk <- cnt$idx - I0[cnt$chrn]
  stopifnot(all(rk >= 1L), all(rk <= lengths(cpg)[cnt$chrn]))
  cnt[, pos := allpos[idx]]
  b <- fread(cmd = sprintf("gzip -dc '%s'", fb), header = FALSE, sep = "\t", col.names = c("chr", "pos", "end", "beta"))
  b <- b[chr %in% AUTO][, chrn := match(chr, AUTO)]
  chk <- merge(cnt[, .(chrn, pos, beta_pat = meth / total)], b[, .(chrn, pos, beta)], by = c("chrn", "pos"))
  agree <- mean(abs(chk$beta_pat - chk$beta) < 1e-6)
  stopifnot(agree > 0.8)                                               # a wrong map would agree rarely
  dbeta <- 1 - agree
  list(cnt = cnt, span = span[, .(chrn, idx, N, linkage)],
       qc = data.table(sample = sid[s], gsm = samples$gsm[s],
                       covered_cpgs = nrow(cnt), cpgs_10x = sum(cnt$total >= MIN_COV),
                       cpgs_in_geo_bed_too = nrow(chk), share_same_beta_as_geo = round(agree, 4),
                       adjacent_pairs_10_reads = nrow(span)))
}

t0 <- Sys.time()
acc <- NULL; map <- NULL; link <- NULL; qc <- list()
for (s in seq_len(nrow(samples))) {
  o <- one_sample(s)
  keep <- o$cnt[total >= MIN_COV, .(chrn, idx, meth, total)]
  setnames(keep, c("meth", "total"), paste0(c("m_", "n_"), s))
  acc <- if (is.null(acc)) keep else merge(acc, keep, by = c("chrn", "idx"))     # stays: >= 10 reads in every sample so far
  mp <- o$cnt[, .(chrn, idx, pos)]
  if (is.null(map)) map <- mp else {
    both <- merge(map, mp, by = c("chrn", "idx"))
    stopifnot(identical(both$pos.x, both$pos.y))                                 # an index has one position in every sample
    map <- unique(rbind(map, mp))
  }
  sp <- o$span[, .(chrn, idx, s1 = linkage, c1 = 1L, n1 = N)]
  link <- if (is.null(link)) sp else {
    z <- merge(link, sp, by = c("chrn", "idx"), all = TRUE, suffixes = c("", ".new"))
    for (v in c("s1", "c1", "n1")) set(z, j = v, value = rowSums(cbind(z[[v]], z[[paste0(v, ".new")]]), na.rm = TRUE))
    z[, c("s1.new", "c1.new", "n1.new") := NULL]; z
  }
  qc[[s]] <- o$qc
  cat(sprintf("sample %2d/%d %s: %d CpGs, %d at >= 10 reads; cohort sites left %d (%.1f min)\n", s, nrow(samples), sid[s],
              o$qc$covered_cpgs, o$qc$cpgs_10x, nrow(acc), as.numeric(Sys.time() - t0, units = "mins")))
  rm(o); invisible(gc())
}
setkey(map, chrn, idx); setkey(acc, chrn, idx)
acc <- map[acc]
M <- as.matrix(acc[, paste0("m_", seq_along(sid)), with = FALSE]); N <- as.matrix(acc[, paste0("n_", seq_along(sid)), with = FALSE])
colnames(M) <- colnames(N) <- sid
link[, `:=`(linkage = s1 / c1, samples = c1, mean_reads = n1 / c1)][, c("s1", "c1", "n1") := NULL]
cache <- list(sites = acc[, .(chrn, idx, pos)], M = M, N = N, map = map, linkage = link, samples = sid,
              built = format(Sys.time()), note = "hg19; pos = 1-based C; sites have >= 10 reads in all 29 samples")
saveRDS(cache, file.path(dat, "cache_cohort.rds"))
qc <- rbindlist(qc)
fwrite(qc, file.path(res, "Task50_BuildSummary.csv"))
cat("\ncohort sites with >= 10 reads in all 29 samples:", nrow(acc), "\nadjacent pairs with linkage in >= 1 sample:",
    nrow(link), "\nshare of shared CpGs with the same beta as GEO, lowest sample:", min(qc$share_same_beta_as_geo), "\ndone\n")
