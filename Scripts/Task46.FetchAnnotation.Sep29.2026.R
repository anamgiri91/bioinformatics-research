# ================================================================
# Task 46 (part 0) - fetch hg38 chr22 annotation for the reliability
#                    factor of the similarity score
# Date: Sep 29, 2026
#
# Downloads, once, into Data/annotation/ (gitignored with Data/):
#   rmsk_chr22.tsv.gz          UCSC RepeatMasker, hg38, chr22 (REST API)
#   bismap_k50_multi_chr22.tsv.gz   Bismap multi-read mappability, 50-mers,
#                              after bisulfite conversion (bigWig, 0..1)
#   bismap_k100_multi_chr22.tsv.gz  the same for 100-mers (sensitivity)
#   bismap_k50_C2T_unique_chr22.tsv.gz  Bismap single-read mappability,
#                              50-mers, forward strand (C->T converted):
#                              regions where a 50-mer maps uniquely
#   annotation_manifest.json   URLs and download times
#
# Bismap: Karimzadeh, Ernst, Kundaje & Hoffman 2018, Nucleic Acids Res,
# "Umap and Bismap: quantifying genome and methylome mappability".
# Our RRBS file is the forward strand, so the C->T track is the one
# that applies to single-read uniqueness.
# ================================================================

suppressPackageStartupMessages({ library(data.table); library(rtracklayer); library(jsonlite) })

pick_dir <- function(...) { for (d in c(...)) if (dir.exists(d)) return(d)
  stop("no candidate directory exists") }
repo_dir <- pick_dir("/mmfs1/home/wln26/Experiments.Outlier.July31.2026",
                     path.expand("~/Desktop/bioinformatics-research"))
ann_dir <- file.path(repo_dir, "Data", "annotation")
dir.create(ann_dir, showWarnings = FALSE, recursive = TRUE)
chr22 <- GRanges("chr22", IRanges(1, 50818468))   # hg38 chr22 length
log <- list()

# RepeatMasker via the UCSC REST API
f_rmsk <- file.path(ann_dir, "rmsk_chr22.tsv.gz")
if (!file.exists(f_rmsk)) {
  url <- "https://api.genome.ucsc.edu/getData/track?genome=hg38;track=rmsk;chrom=chr22;maxItemsOutput=-1"
  j <- fromJSON(url)
  rm <- as.data.table(j$rmsk)[, .(chrom = genoName, start = genoStart, end = genoEnd, strand,
                                  repName, repClass, repFamily, milliDiv)]
  stopifnot(nrow(rm) > 1000, all(rm$chrom == "chr22"))
  fwrite(rm, f_rmsk)
  log$rmsk <- list(url = url, items = nrow(rm), downloaded = format(Sys.time()),
                   ucsc_data_time = j$dataTime)
}

# Bismap tracks, read straight from UCSC's big files for chr22 only
base <- "https://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/"
get_track <- function(file, out, type) {
  f <- file.path(ann_dir, out)
  if (file.exists(f)) return(invisible(NULL))
  url <- paste0(base, file)
  g <- if (type == "bw") import(BigWigFile(url), which = chr22) else import(BigBedFile(url), which = chr22)
  d <- data.table(chrom = as.character(seqnames(g)), start = start(g), end = end(g))
  if (type == "bw") d[, score := g$score]
  stopifnot(nrow(d) > 1000)
  fwrite(d, f)
  log[[out]] <<- list(url = url, rows = nrow(d), downloaded = format(Sys.time()))
}
get_track("k50.Bismap.MultiTrackMappability.bw",  "bismap_k50_multi_chr22.tsv.gz",  "bw")
get_track("k100.Bismap.MultiTrackMappability.bw", "bismap_k100_multi_chr22.tsv.gz", "bw")
get_track("k50.C2T-Converted.bb", "bismap_k50_C2T_unique_chr22.tsv.gz", "bb")

if (length(log)) write_json(log, file.path(ann_dir, "annotation_manifest.json"), pretty = TRUE, auto_unbox = TRUE)
for (f in list.files(ann_dir, full.names = TRUE)) cat(sprintf("%-40s %8.1f KB\n", basename(f), file.size(f) / 1024))
cat("done\n")
