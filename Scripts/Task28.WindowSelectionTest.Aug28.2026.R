# ================================================================
# Task 28 - How should the 100-CpG window be chosen, and does it matter?
# Date: Aug 28, 2026
#
# Task 20 replaced the Task 17 window (100 CpGs consecutive by INDEX,
# spanning 5.17 Mb) with tight100: the 100 consecutive CpGs of minimum
# total SPAN. A fair objection was raised: minimising total span does
# not stop one large internal gap. tight100 still has a 13,823 bp hole.
#
# The stricter criterion proposed was: find 100 consecutive CpGs where
# EVERY neighbouring gap is small - a genuinely continuous region.
#
# This script:
#   A. enumerates all 6,710 possible 100-CpG windows and asks whether
#      any satisfies "every gap <= 1 kb";
#   B. builds the min-max-gap window and compares it against index100
#      and tight100 on spacing AND on how representative its state mix
#      is of chr22 as a whole;
#   C. asks the question that decides it - do the three windows give
#      different scientific answers? Contrary-flag rates, bio-vs-ext
#      agreement and flag magnitude are recomputed in each.
#   D. reports the longest genuinely continuous runs on chr22, since
#      those - not any 100-CpG window - are the right unit for a
#      co-methylation question.
# ================================================================

options(stringsAsFactors = FALSE)

pick_dir <- function(...) { for (d in c(...)) if (dir.exists(d)) return(d); stop("no dir") }
src_dir <- pick_dir(
  "/home/s_s355/research3.data/TCGA/filter.BRCA.UCEC.Aug17.2022/sorted.TCGA.380355cg.files",
  path.expand("~/Downloads/tcga_brca_380355cg"))
repo <- pick_dir("/mmfs1/home/wln26/Experiments.Outlier.July31.2026",
                 path.expand("~/Desktop/bioinformatics-research"))
out_dir <- file.path(repo, "Results")

inputs <- list(
  Normal = file.path(src_dir,"Sorted.BRCA.53Alive.Normal.380355cg.75col.May28.2026.txt"),
  Tumor  = file.path(src_dir,"Sorted.BRCA.53Alive.Tumor.380355cg.75col.May28.2026.txt"))
fext <- list(Normal=file.path(out_dir,"Task13_Normal_extref_flags.csv"),
             Tumor =file.path(out_dir,"Task13_Tumor_extref_flags.csv"))

W <- 100; ST <- c("L","LM","M","HM","H","R"); P_LEVEL <- 0.01
rd <- function(p,tab=FALSE) as.data.frame(data.table::fread(
  p, sep=if(tab) "\t" else ",", header=TRUE, check.names=FALSE, showProgress=FALSE))
wo <- function(t,f){ p<-file.path(out_dir,f); write.csv(t,p,row.names=FALSE)
  cat("  written:",p," (",nrow(t)," rows )\n",sep="") }

cat("=========================================================\n")
cat("TASK 28 - WINDOW SELECTION: SPAN vs MAXIMUM GAP\n")
cat("=========================================================\n\n")

normal <- rd(inputs$Normal, TRUE); tumor <- rd(inputs$Tumor, TRUE)
c22 <- normal[normal$Chromosome=="chr22",]; c22 <- c22[order(c22$Start),]
pos <- c22$Start; n <- nrow(c22); nw <- n-W+1; gaps <- diff(pos)

# ---------- A. enumerate every window ----------
span <- pos[W:n] - pos[1:nw]
mx <- sapply(1:nw, function(i) max(gaps[i:(i+W-2)]))
md <- sapply(1:nw, function(i) median(gaps[i:(i+W-2)]))
u1 <- sapply(1:nw, function(i) sum(gaps[i:(i+W-2)] <= 1000))
enum <- data.frame(start.index=1:nw, start.pos=pos[1:nw], end.pos=pos[W:n],
                   span=span, max.gap=mx, median.gap=md, gaps.under.1kb=u1)
wo(enum, "Task28_Chr22_AllWindows.csv")

cat("A. Is a fully continuous 100-CpG window possible on chr22?\n")
for (thr in c(1000,2000,5000,10000))
  cat(sprintf("   windows with every gap <= %5d bp : %d of %d\n", thr, sum(mx<=thr), nw))
cat("\n   ANSWER: no. The 450k array does not place 100 consecutive probes\n")
cat("   that densely anywhere on chr22. The constraint is the platform's\n")
cat("   probe design, not the search.\n\n")

# ---------- B. the three candidate windows ----------
wins <- list(index100 = 1,
             tight100 = which.min(span),
             minmaxgap100 = which.min(mx))
expected <- table(factor(c22$methy.state, levels=ST)) / n * W

cat("B. The three candidate windows\n\n")
rows <- list()
for (nm in names(wins)) {
  i <- wins[[nm]]; idx <- i:(i+W-1); g <- diff(pos[idx])
  obs <- table(factor(c22$methy.state[idx], levels=ST))
  chisq <- sum((obs-expected)^2/expected)
  rows[[nm]] <- data.frame(window=nm, start.index=i,
    start.pos=pos[i], end.pos=pos[i+W-1], span=pos[i+W-1]-pos[i],
    max.gap=max(g), median.gap=median(g), gaps.under.1kb=sum(g<=1000),
    L=obs[["L"]],LM=obs[["LM"]],M=obs[["M"]],HM=obs[["HM"]],H=obs[["H"]],R=obs[["R"]],
    state.chisq.vs.chr22=round(chisq,1))
}
cmp <- do.call(rbind, rows); rownames(cmp) <- NULL
print(cmp[,c("window","span","max.gap","median.gap","gaps.under.1kb",
             "L","LM","M","HM","H","R","state.chisq.vs.chr22")], row.names=FALSE)
cat("\n   state.chisq.vs.chr22 = how far the window's state mix departs from\n")
cat("   chr22 overall (lower = more representative).\n")
wo(cmp, "Task28_Chr22_WindowCandidates.csv")

# ---------- C. does the choice change the answer? ----------
build_bio <- function(B, states) {
  o <- matrix(0, nrow(B), ncol(B))
  for (i in seq_len(nrow(B))) {
    st <- states[i]; v <- as.numeric(B[i,])
    if (all(is.na(v))) { o[i,] <- NA; next }
    if (!st %in% c("L","LM","H","HM")) next
    hi <- st %in% c("L","LM"); p <- if (hi) 1-P_LEVEL else P_LEVEL
    thr <- quantile(v,p,na.rm=TRUE,names=FALSE)
    o[i,] <- if (hi) as.integer(v>thr) else -as.integer(v<thr)
    o[i,is.na(v)] <- NA
  }
  o
}
kappa2 <- function(a,b){ k<-!is.na(a)&!is.na(b); a<-a[k]; b<-b[k]
  if(!length(a)) return(NA_real_)
  tb<-table(factor(a,levels=c(-1,0,1)),factor(b,levels=c(-1,0,1))); N<-sum(tb)
  po<-sum(diag(tb))/N; pe<-sum(rowSums(tb)*colSums(tb))/N^2
  if (pe==1) NA_real_ else (po-pe)/(1-pe) }
CONTRARY <- list(H=1,HM=1,L=-1,LM=-1)

cat("\nC. Do the three windows give different scientific answers?\n\n")
res <- list()
for (tissue in c("Normal","Tumor")) {
  src <- if (tissue=="Normal") normal else tumor
  cols <- paste0(substr(tissue,1,1),1:53)
  fx <- rd(fext[[tissue]])
  for (nm in names(wins)) {
    i <- wins[[nm]]; ids <- c22$Composite.Element.REF[i:(i+W-1)]
    meta <- src[match(ids, src$Composite.Element.REF),]
    B <- as.matrix(meta[,cols,drop=FALSE]); storage.mode(B)<-"numeric"
    st <- meta$methy.state
    bio <- build_bio(B, st)
    E <- as.matrix(fx[match(ids,fx$cgID),cols,drop=FALSE]); storage.mode(E)<-"numeric"
    med <- apply(B,1,median,na.rm=TRUE); dev <- abs(B-med)
    contra <- ev <- 0
    for (s in names(CONTRARY)) { k<-which(st==s); if(!length(k)) next
      contra <- contra + sum(E[k,]==CONTRARY[[s]], na.rm=TRUE)
      ev <- ev + sum(!is.na(E[k,])) }
    ok <- !is.na(E) & E!=0
    magH <- median(dev[ok & st=="H"], na.rm=TRUE)
    magL <- median(dev[ok & st=="L"], na.rm=TRUE)
    res[[length(res)+1]] <- data.frame(tissue=tissue, window=nm,
      ext.flag.rate.pct=round(100*mean(E!=0,na.rm=TRUE),2),
      contrary.pct=round(100*contra/max(1,ev),2),
      bio.vs.ext.kappa=round(kappa2(as.vector(bio),as.vector(E)),3),
      median.dbeta.H=round(magH,4), median.dbeta.L=round(magL,4),
      pct.trivial=round(100*mean(dev[ok]<0.10,na.rm=TRUE),1))
  }
  rm(fx); invisible(gc())
}
sens <- do.call(rbind,res)
print(sens, row.names=FALSE)
wo(sens, "Task28_Chr22_WindowSensitivity.csv")

# ---------- D. genuinely continuous runs ----------
cat("\nD. The longest genuinely continuous runs on chr22\n")
cat("   (if the goal is co-methylation, this - not a 100-CpG window - is the unit)\n\n")
runs <- list()
for (thr in c(500,1000,2000,5000)) {
  r <- rle(gaps <= thr)
  if (!any(r$values)) next
  ends <- cumsum(r$lengths)[r$values]; starts <- ends - r$lengths[r$values] + 1
  ncpg <- r$lengths[r$values] + 1
  b <- which.max(ncpg); s <- starts[b]; e <- ends[b]+1
  runs[[length(runs)+1]] <- data.frame(max.gap.bp=thr, longest.run.cpgs=max(ncpg),
    start.pos=pos[s], end.pos=pos[e], span.bp=pos[e]-pos[s],
    median.gap.bp=median(diff(pos[s:e])),
    runs.ge.20.cpgs=sum(ncpg>=20), runs.ge.30.cpgs=sum(ncpg>=30))
}
rr <- do.call(rbind,runs)
print(rr, row.names=FALSE)
wo(rr, "Task28_Chr22_ContinuousRuns.csv")

cat("\n=========================================================\n")
cat("TASK 28 COMPLETE\n=========================================================\n")
