# ================================================================
# Task 40 - OutlierMeth checks for the 2026-09-28 email
#           item 1: probe-sequence problems and the beta offset
#           item 3: does the output depend on the CpG ID?
# Date: Sep 28, 2026
#
# THE REQUEST
#
#   1. Check OutlierMeth (Downs, Thursby & Cope 2023, Epigenetics
#      18(1):2213874): (1) did the authors consider the reference
#      sequence (repetitive regions, other sequence features)? (2) what
#      does "outlier" mean to them? (3) did they consider the offset in
#      beta = M / (M + U + 100)?
#   3. Run OutlierMeth (1) with the CG ID replaced by the position and
#      (2) with the CG ID removed, and see if the output is exactly the
#      same.
#
#   Items 1(1)-(3) are answered from the paper. This script adds the
#   data behind two of those answers, and runs item 3.
#
# WHAT IS RUN
#
#   Part A (item 1(1)). The paper never masks probes. How many of the
#     Task 13 ext flags sit on probes that Zhou et al. 2017 mask? The
#     mask file has no repeat column; M_nonuniq (the probe sequence also
#     matches another place in the genome) is the closest to
#     "repetitive region".
#   Part B (item 1(3)). Arithmetic only: how far the +100 offset moves
#     beta for a given true methylation fraction and total signal
#     M + U. The files hold beta only, so the offset cannot be tested
#     on this data.
#   Part C (item 3). The package, sourced from Packages/OutlierMeth
#     (it is not installed on this machine), run on all 380,355 CpGs
#     with four kinds of row name:
#       cgID           as in Task 13; must reproduce Task 13 exactly
#       chr:position   e.g. "chr1:15865"; unique
#       position only  e.g. "15865"; 860 CpGs share a position with a
#                      CpG on another chromosome
#       none           row names removed
#     Each is run two ways: Task 13's call,
#       flagMeth(beta[rownames(ref), ], reference = ref, p = 0.01),
#     and the plain call flagMeth(beta, reference = ref, p = 0.01).
#     Self reference (referenceMeth on the same 53 samples), Normal and
#     Tumour. The external panel (tcga.rda) is not on this machine, so
#     the external case uses a stand-in with the same structure: a
#     reference built from the 32 Dead Normal samples, keyed by cgID,
#     applied to the 53 Alive Normal samples.
# ================================================================

suppressPackageStartupMessages(library(data.table))

pick_dir <- function(...) { for (d in c(...)) if (dir.exists(d)) return(d)
  stop("no candidate directory exists") }
src_dir <- pick_dir(
  "/home/s_s355/research3.data/TCGA/filter.BRCA.UCEC.Aug17.2022/sorted.TCGA.380355cg.files",
  path.expand("~/Downloads/tcga_brca_380355cg"))
repo_dir <- pick_dir("/mmfs1/home/wln26/Experiments.Outlier.July31.2026",
                     path.expand("~/Desktop/bioinformatics-research"))
out_dir <- file.path(repo_dir, "Results")
pkg_dir <- file.path(repo_dir, "Packages", "OutlierMeth")
r4 <- function(x) round(x, 4)

# the package is not installed here; its four functions are plain R,
# so sourcing them runs the same code
for (f in list.files(file.path(pkg_dir, "R"), pattern = "\\.R$", full.names = TRUE)) source(f)
cat("OutlierMeth version", read.dcf(file.path(pkg_dir, "DESCRIPTION"), fields = "Version")[1], "\n")

src_file <- function(cohort, tissue) file.path(src_dir, if (cohort == "Alive")
  sprintf("Sorted.BRCA.53Alive.%s.380355cg.75col.May28.2026.txt", tissue) else
  sprintf("Sorted.BRCA.32Dead.%s.380355cg.54col.May28.2026.txt", tissue))
read_beta <- function(cohort, tissue) {
  d <- fread(src_file(cohort, tissue))
  sc <- 5:(which(names(d) == "Min") - 1)
  b <- as.matrix(d[, ..sc])
  rownames(b) <- d$Composite.Element.REF
  list(beta = b, meta = d[, .(cgID = Composite.Element.REF, chr = Chromosome, pos = Start)])
}
flag_file <- function(tissue, ref) file.path(out_dir, sprintf("Task13_%s_%sref_flags.csv", tissue, ref))
flag_matrix <- function(f) as.matrix(f[, grep("^[NT][0-9]+$", names(f), value = TRUE), with = FALSE])

# ----------------------------------------------------------------
# Part A. ext flags on probes with sequence problems (Zhou et al. 2017)
# ----------------------------------------------------------------
msk <- fread(file.path(repo_dir, "Data", "zhou_HM450", "HM450.hg38.mask.tsv.gz"))
groups <- list(
  "maps to more than one place (M_nonuniq)"   = msk[grepl("M_nonuniq", maskUniq), Probe_ID],
  "common SNP in the probe (SNPcommon_5pt)"   = msk[grepl("SNPcommon_5pt", maskUniq), Probe_ID],
  "any problem in the general mask (M_general)" = msk[M_general == TRUE, Probe_ID])
partA <- rbindlist(lapply(c("Normal", "Tumor"), function(tissue) {
  f <- fread(flag_file(tissue, "ext"))
  fm <- flag_matrix(f)
  n_eval <- rowSums(!is.na(fm)); n_flag <- rowSums(fm != 0, na.rm = TRUE)
  clean <- !(f$cgID %in% groups[[3]])
  rate_clean <- sum(n_flag[clean]) / sum(n_eval[clean])
  rbindlist(lapply(c(names(groups), "not in the general mask"), function(g) {
    in_g <- if (g == "not in the general mask") clean else f$cgID %in% groups[[g]]
    data.table(tissue, probes = g, CpGs = sum(in_g), share_of_CpGs = r4(mean(in_g)),
               ext_flags = sum(n_flag[in_g]),
               share_of_ext_flags = r4(sum(n_flag[in_g]) / sum(n_flag)),
               flag_rate = r4(sum(n_flag[in_g]) / sum(n_eval[in_g])),
               rate_vs_unmasked = round(sum(n_flag[in_g]) / sum(n_eval[in_g]) / rate_clean, 2))
  }))
}))
fwrite(partA, file.path(out_dir, "Task40_MaskedProbeFlags.csv"))
cat("\nPart A. ext flags on masked probes (Task 13 flag matrices):\n"); print(partA)

# ----------------------------------------------------------------
# Part B. how far the +100 offset moves beta (arithmetic)
# ----------------------------------------------------------------
# true_beta = M / (M + U); with the offset, beta = M / (M + U + 100)
#           = true_beta * total / (total + 100), where total = M + U
partB <- CJ(total_signal = c(500, 1000, 2000, 5000, 10000, 20000),
            true_beta = c(0.05, 0.5, 0.95))
partB[, beta_with_offset := r4(true_beta * total_signal / (total_signal + 100))]
partB[, shift_down := r4(true_beta - beta_with_offset)]
fwrite(partB, file.path(out_dir, "Task40_OffsetArithmetic.csv"))
cat("\nPart B. beta with and without the +100 offset:\n"); print(partB)

# ----------------------------------------------------------------
# Part C. does the output change when the CpG ID changes?
# ----------------------------------------------------------------
catch <- function(expr) tryCatch(expr, error = function(e) e)
pad_task13 <- function(beta, flag) {       # Task 13's padding step
  full <- matrix(NA, nrow = nrow(beta), ncol = ncol(beta), dimnames = dimnames(beta))
  full[rownames(flag), ] <- flag
  full
}
describe <- function(res, base) {
  if (inherits(res, "error"))
    return(list(result = paste("error:", conditionMessage(res))))
  same <- (is.na(res) & is.na(base)) | (!is.na(res) & !is.na(base) & res == base)
  list(result = "ran", cells = length(res), identical_cells = sum(same),
       different_cells = sum(!same), rows_with_a_difference = sum(rowSums(!same) > 0),
       hypo_flags = sum(res == -1, na.rm = TRUE), hyper_flags = sum(res == 1, na.rm = TRUE),
       NA_cells = sum(is.na(res)))
}
name_rows <- function(b, meta, ids) {
  rownames(b) <- switch(ids, "cgID" = meta$cgID,
                        "chr:position" = paste0(meta$chr, ":", meta$pos),
                        "position only" = as.character(meta$pos),
                        "none" = NULL)
  b
}
id_kinds <- c("cgID", "chr:position", "position only", "none")

partC <- list()
for (tissue in c("Normal", "Tumor")) {
  x <- read_beta("Alive", tissue)
  stopifnot(anyDuplicated(paste(x$meta$chr, x$meta$pos)) == 0, is.integer(x$meta$pos))
  base <- NULL
  for (ids in id_kinds) {
    t0 <- Sys.time()
    b <- name_rows(x$beta, x$meta, ids)
    ref <- catch(referenceMeth(b))
    ref_note <- if (inherits(ref, "error")) paste("error:", conditionMessage(ref)) else
      paste0(nrow(ref), " rows; first names ", paste(head(rownames(ref), 2), collapse = ", "))
    matched <- if (inherits(ref, "error")) NA_integer_ else length(intersect(rownames(b), rownames(ref)))
    t13  <- if (inherits(ref, "error")) ref else catch(pad_task13(b, flagMeth(b[rownames(ref), ], reference = ref, p = 0.01)))
    flat <- if (inherits(ref, "error")) ref else catch(flagMeth(b, reference = ref, p = 0.01))
    if (ids == "cgID") {
      base <- t13
      saved <- fread(flag_file(tissue, "self"))
      stopifnot(identical(saved$cgID, rownames(base)))
      partC[[length(partC) + 1]] <- c(list(tissue = tissue, reference = "self (53 samples)",
        row_names = "cgID, this run", call = "vs the saved Task 13 file", reference_rows = ref_note,
        rows_matched = matched), describe(base, flag_matrix(saved)))
    }
    for (call in c("Task 13 call", "plain call")) {
      res <- if (call == "Task 13 call") t13 else flat
      partC[[length(partC) + 1]] <- c(list(tissue = tissue, reference = "self (53 samples)",
        row_names = ids, call = call, reference_rows = ref_note, rows_matched = matched),
        describe(res, base))
    }
    cat(tissue, ids, "done in", round(as.numeric(Sys.time() - t0, units = "secs")), "s\n")
    rm(b, ref, t13, flat); invisible(gc())
  }
  if (tissue == "Normal") {
    # external case, stand-in panel: 32 Dead Normal samples keyed by cgID
    dead <- read_beta("Dead", "Normal")
    stopifnot(identical(rownames(dead$beta), rownames(x$beta)))
    ref_dead <- referenceMeth(dead$beta)
    ext_base <- flagMeth(x$beta, reference = ref_dead, p = 0.01)
    ref_dead_pos <- ref_dead
    rownames(ref_dead_pos) <- paste0(x$meta$chr, ":", x$meta$pos)[match(rownames(ref_dead), x$meta$cgID)]
    ext_runs <- list(
      list("cgID", "panel keyed by cgID", ref_dead),
      list("chr:position", "panel keyed by cgID", ref_dead),
      list("chr:position", "panel renamed to chr:position", ref_dead_pos),
      list("none", "panel keyed by cgID", ref_dead))
    for (er in ext_runs) {
      b <- name_rows(x$beta, x$meta, er[[1]])
      res <- catch(flagMeth(b, reference = er[[3]], p = 0.01))
      partC[[length(partC) + 1]] <- c(list(tissue = tissue, reference = "external stand-in (32 Dead Normal)",
        row_names = er[[1]], call = er[[2]], reference_rows = paste(nrow(er[[3]]), "rows"),
        rows_matched = length(intersect(rownames(b), rownames(er[[3]])))),
        describe(res, ext_base))
    }
    rm(dead, ref_dead, ref_dead_pos, ext_base, b, res)
  }
  rm(x, base); invisible(gc())
}
partC <- rbindlist(partC, fill = TRUE)
fwrite(partC, file.path(out_dir, "Task40_IDTest.csv"))
cat("\nPart C. OutlierMeth output with the CpG ID replaced or removed:\n")
print(partC[, !"reference_rows"])
print(partC[, .(tissue, row_names, call, reference_rows)])
cat("\ndone\n")
