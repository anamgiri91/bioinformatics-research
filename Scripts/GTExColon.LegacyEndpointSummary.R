# Replay the identical hypothesis-summary code using the retained legacy gain.
# The frozen score, thresholds, fold allocation and selected pairs are unchanged.
suppressPackageStartupMessages(library(data.table))
arg <- grep("^--file=",commandArgs(FALSE),value=TRUE)
root <- dirname(dirname(normalizePath(sub("^--file=","",arg[1L]))))
res <- file.path(root,"Results"); r4 <- function(x) round(x,4)
fr <- fread(file.path(res,"Task51_Frontier.csv"))
fr_folds <- fread(file.path(res,"Task51_Frontier_Folds.csv"))
po <- readRDS(file.path(root,"Data","gtex_colon_rrbs","pair_outcomes.rds"))
fr[,gain:=gain_legacy]; fr_folds[,gain:=gain_legacy]
po[,heldout_gain:=heldout_gain_legacy]
pairs <- po
SUBSETS <- list(all=rep(TRUE,nrow(po)),first=po$half=="first",second=po$half=="second")
src <- readLines(file.path(root,"Scripts","Task51.GTExColon.FrozenTests.Sep29.2026.R"))
start <- grep("^QUESTIONS <-",src); end <- grep("^decay <-",src)
stopifnot(length(start)==1L,length(end)==1L,start<end)
section <- src[start:(end-1L)]
section <- gsub('"Task51_', '"Task51_LegacyEndpoint_',section,fixed=TRUE)
eval(parse(text=section))
cat("Legacy endpoint hypothesis summaries saved separately\n")
