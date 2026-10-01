# Final consistency checks and provenance; does not rerun the cohort analysis.
suppressPackageStartupMessages({library(data.table);library(jsonlite)})
arg <- grep("^--file=",commandArgs(FALSE),value=TRUE)
root <- dirname(dirname(normalizePath(sub("^--file=","",arg[1L]))))
res <- file.path(root,"Results")
read <- function(n) fread(file.path(res,n))
proto <- fromJSON(file.path(res,"Task47_FrozenProtocol.json"))
for(f in names(proto$files_sha256))
  stopifnot(digest::digest(file=file.path(root,f),algo="sha256")==proto$files_sha256[[f]])
ff <- read("Task51_Frontier_Folds.csv"); fr <- read("Task51_Frontier.csv")
stopifnot(nrow(ff)==10*18*7*3,!anyDuplicated(ff[,.(subset,family,ell,fold)]),
          setequal(ff$fold,1:10),all(is.finite(ff$error)),all(is.finite(ff$gain)),
          all(ff$error>=0 & ff$error<=1),all(ff$gain==ff$gain_legacy))
expected <- ff[,.(error=round(mean(error),5),gain=round(mean(gain),5)),by=.(subset,family,ell)]
joined <- merge(fr,expected,by=c("subset","family","ell"),suffixes=c("","_check"))
stopifnot(nrow(joined)==nrow(fr),max(abs(joined$error-joined$error_check))<1e-12,
          max(abs(joined$gain-joined$gain_check))<1e-12)
h <- read("Task51_Hypotheses.csv"); old <- read("Task51_LegacyEndpoint_Hypotheses.csv")
stopifnot(isTRUE(all.equal(h,old,check.attributes=FALSE)))
po <- readRDS(file.path(root,"Data","gtex_colon_rrbs","pair_outcomes.rds"))
stopifnot(nrow(po)==1574581,!anyDuplicated(po[,.(chrn,idx)]),
          all(is.finite(po$heldout_gain)),all(po$heldout_gain==po$heldout_gain_legacy),
          sum(po$flag_not_unique,na.rm=TRUE)==146998,
          sum(is.na(po$flag_not_unique))==414)
h8 <- read("Task52_H8_Detail.csv"); h9 <- read("Task52_H9_Detail.csv")
stopifnot(all(h8$holds==(h8$mean_linkage>0 & h8$spearman_with_distance<0)),
          all(h9$holds==(h9$difference>0)),
          h8[subset=="all",pairs]==sum(h8[subset!="all",pairs]),
          h9[subset=="all",pairs]==sum(h9[subset!="all",pairs]),
          sum(read("Task52_H8_ByDistance.csv")$pairs)==h8[subset=="all",pairs],
          all(h8$lo95_mean<=h8$hi95_mean),all(h9$lo95<=h9$hi95))
source(file.path(root,"Scripts","GTExColon.BlockBootstrap.R"))
check_block_bootstrap()
scripts <- c("Task47.FrozenScore.Sep29.2026.R","Task50.GTExColon.Build.Sep29.2026.R",
  "Task51.GTExColon.FrozenTests.Sep29.2026.R","Task52.GTExColon.SharedReads.Sep29.2026.R",
  "GTExColon.BlockBootstrap.R","GTExColon.Continuation.Checks.R",
  "GTExColon.LegacyEndpointSummary.R","GTExColon.Validation.Checks.R",
  "Item5.SharedReadNoise.Pilot.R")
for(f in scripts) invisible(parse(file.path(root,"Scripts",f)))
tables <- list.files(res,pattern="^(Task5[012]_.*\\.csv|GTExColon_.*\\.csv|Item5_SharedReadNoise_.*\\.(csv|json))$",full.names=TRUE)
files <- c(file.path(root,"Scripts",scripts),tables,file.path(root,"report.md"),
  file.path(root,"Results",c("Task47_FrozenProtocol.md","Task47_FrozenProtocol.json",
    "GTExColon_Validation_Report.md","GTExColon_Continuation_Deviations.md",
    "Task51_H6_AnnotationFix.md","Item5_SharedReadNoise_Design.md",
    "Fig37_GTExColon_Frontier.png","Fig38_SharedReads.png")),
  file.path(root,"Data","gtex_colon_rrbs",c("colon_samples.json","cache_cohort.rds","pair_outcomes.rds")))
stopifnot(all(file.exists(files)))
files <- unique(files)
manifest <- data.table(path=substring(files,nchar(root)+2L),bytes=file.info(files)$size,
  sha256=vapply(files,function(f) digest::digest(file=f,algo="sha256"),""))
write_json(list(checked_at=format(Sys.time(),tz="UTC",usetz=TRUE),status="all checks passed",
  frozen_score_unchanged=TRUE,files=manifest),file.path(res,"GTExColon_Validation_Manifest.json"),pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),file.path(res,"GTExColon_Validation_SessionInfo.txt"))
cat("PASS: frozen hashes; 3,780 fold records; fold means; identical legacy results; corrected annotation counts; H8/H9 summaries; bootstrap checks; script parsing\n")
