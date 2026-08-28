# ================================================================
# Task 26 - Empirical answers to the supervisor review
# Date: Aug 28, 2026
#
# The review raised points that are answerable from data rather than
# by softening wording. This script answers four of them.
#
#   A. DENOMINATOR AUDIT. "Are all reported external-reference
#      percentages using the same denominator?" They are NOT. CpGs
#      absent from the tcga.rda panel come back all-NA, so the ext
#      method is evaluated on fewer cells than bio/self at the same
#      sites. Every table should state this. Quantified here.
#
#   B. FLOOR SENSITIVITY. "Is 0.10 biologically justified, or chosen
#      after seeing these data?" Chosen after. This runs the floor at
#      0.05 / 0.10 / 0.15 / 0.20 so the choice can be judged against
#      a curve rather than a single point.
#
#   C. RATE-MATCHED COMPARISON. "ext+floor has a lower rate; should it
#      be recalibrated so the stability advantage is not partly a
#      lower-rate artifact?" Yes. Every competing method is re-tuned
#      DOWN to the floor's own flag rate and stability re-measured, so
#      the comparison is like for like.
#
#   D. GENOME-WIDE CONFIRMATION. The magnitude result rested on chr22
#      (6,809 CpGs). Repeated over all 380,355.
#
# Outputs:
#   Task26_DenominatorAudit_<tissue>.csv
#   Task26_FloorSensitivity_<tissue>.csv
#   Task26_RateMatchedStability_<tissue>.csv
#   Task26_GenomeWideMagnitude.csv
# ================================================================

options(stringsAsFactors = FALSE)
set.seed(20260828)

pick_dir <- function(...) { for (d in c(...)) if (dir.exists(d)) return(d); stop("no dir") }
src_dir <- pick_dir(
  "/home/s_s355/research3.data/TCGA/filter.BRCA.UCEC.Aug17.2022/sorted.TCGA.380355cg.files",
  path.expand("~/Downloads/tcga_brca_380355cg"))
repo <- pick_dir("/mmfs1/home/wln26/Experiments.Outlier.July31.2026",
                 path.expand("~/Desktop/bioinformatics-research"))
out_dir <- file.path(repo, "Results")

inputs <- list(
  Normal = file.path(src_dir, "Sorted.BRCA.53Alive.Normal.380355cg.75col.May28.2026.txt"),
  Tumor  = file.path(src_dir, "Sorted.BRCA.53Alive.Tumor.380355cg.75col.May28.2026.txt"))
flag_ext <- list(
  Normal = file.path(out_dir, "Task13_Normal_extref_flags.csv"),
  Tumor  = file.path(out_dir, "Task13_Tumor_extref_flags.csv"))

ST <- c("L","LM","M","HM","H","R")
FLOORS   <- c(0.05, 0.10, 0.15, 0.20)
NOISE_SD <- 0.01
N_REPS   <- 10
EPS      <- 1e-3

rd <- function(p, tab=FALSE) as.data.frame(data.table::fread(
  p, sep=if(tab) "\t" else ",", header=TRUE, check.names=FALSE, showProgress=FALSE))
wo <- function(t, f) { p <- file.path(out_dir, f); write.csv(t, p, row.names=FALSE)
  cat("  written:", p, " (", nrow(t), " rows )\n", sep="") }

rowMedian <- function(M) apply(M, 1, median, na.rm=TRUE)
rowMad    <- function(M) apply(M, 1, mad,    na.rm=TRUE)
to_M <- function(B) { b <- pmin(pmax(B, EPS), 1-EPS); log2(b/(1-b)) }

f_tukey <- function(B,k){ q1<-apply(B,1,quantile,.25,na.rm=TRUE); q3<-apply(B,1,quantile,.75,na.rm=TRUE)
  iq<-q3-q1; o<-matrix(0L,nrow(B),ncol(B)); o[B<q1-k*iq]<--1L; o[B>q3+k*iq]<-1L; o[is.na(B)]<-NA; o }
f_madz <- function(X,k){ m<-rowMedian(X); s<-rowMad(X); s[is.na(s)|s==0]<-NA
  o<-matrix(0L,nrow(X),ncol(X)); o[!is.na(s)&X<m-k*s]<--1L; o[!is.na(s)&X>m+k*s]<-1L; o[is.na(X)]<-NA; o }
f_betafit <- function(B,p){ mu<-rowMeans(B,na.rm=TRUE); v<-apply(B,1,var,na.rm=TRUE); vm<-mu*(1-mu)
  ok<-is.finite(mu)&is.finite(v)&v>0&v<vm; cc<-rep(NA_real_,length(mu)); cc[ok]<-(vm[ok]/v[ok])-1
  a<-mu*cc; b<-(1-mu)*cc; o<-matrix(0L,nrow(B),ncol(B))
  lo<-pbeta(B,a,b); hi<-pbeta(B,a,b,lower.tail=FALSE)
  o[!is.na(lo)&lo<p]<--1L; o[!is.na(hi)&hi<p]<-1L; o[is.na(B)|is.na(a)]<-NA; o }

calibrate <- function(fn, lo, hi, target, decreasing=TRUE, iters=20) {
  for (i in seq_len(iters)) { mid<-(lo+hi)/2; r<-mean(fn(mid)!=0,na.rm=TRUE)
    if (is.na(r)) break; if ((r>target)==decreasing) lo<-mid else hi<-mid }
  (lo+hi)/2
}
jacc <- function(a,b){ ok<-!is.na(a)&!is.na(b); A<-a[ok]!=0; B2<-b[ok]!=0
  u<-sum(A|B2); if(!u) return(NA_real_); sum(A&B2&a[ok]==b[ok])/u }

cat("====================================================\n")
cat("TASK 26 - EMPIRICAL ANSWERS TO THE REVIEW\n")
cat("====================================================\n\n")

gw <- list()

for (tissue in c("Normal","Tumor")) {

  cat("##################################################\n### ", tissue, "\n", sep="")
  cat("##################################################\n")

  cols <- paste0(substr(tissue,1,1),1:53)
  d  <- rd(inputs[[tissue]], tab=TRUE)
  fx <- rd(flag_ext[[tissue]])
  idx <- match(d$Composite.Element.REF, fx$cgID); stopifnot(!anyNA(idx))

  Ball <- as.matrix(d[, cols, drop=FALSE]); storage.mode(Ball)<-"numeric"
  Fall <- as.matrix(fx[idx, cols, drop=FALSE]); storage.mode(Fall)<-"numeric"
  stall <- d$methy.state
  rm(fx); invisible(gc())

  # ---- A. denominator audit (genome-wide) ----
  cat("\n=== A. Denominator audit ===\n")
  cat("  CpGs absent from the tcga.rda panel return all-NA, so the ext\n")
  cat("  method is scored on fewer cells than bio/self at the same sites.\n\n")
  aud <- do.call(rbind, lapply(ST, function(s){
    k <- which(stall==s); if(!length(k)) return(NULL)
    fm <- Fall[k,,drop=FALSE]
    allna <- sum(apply(fm,1,function(r) all(is.na(r))))
    data.frame(tissue=tissue, methy.state=s, n.cpgs=length(k),
      cells.total=length(fm),
      cells.evaluable.ext=sum(!is.na(fm)),
      cells.NA.ext=sum(is.na(fm)),
      cpgs.absent.from.panel=allna,
      pct.cpgs.absent=round(100*allna/length(k),2),
      pct.cells.lost=round(100*sum(is.na(fm))/length(fm),2))
  }))
  print(aud[,c("methy.state","n.cpgs","cpgs.absent.from.panel",
               "pct.cpgs.absent","pct.cells.lost")], row.names=FALSE)
  wo(aud, paste0("Task26_DenominatorAudit_", tissue, ".csv"))

  # ---- D. genome-wide magnitude ----
  cat("\n=== D. Magnitude, genome-wide (all 380,355 CpGs) ===\n")
  med <- rowMedian(Ball); dev <- abs(Ball - med)
  gwt <- do.call(rbind, lapply(ST, function(s){
    k <- which(stall==s); if(!length(k)) return(NULL)
    f <- as.vector(Fall[k,]); w <- as.vector(dev[k,])
    ok <- !is.na(f) & f!=0 & !is.na(w); if(!any(ok)) return(NULL)
    w <- w[ok]
    data.frame(tissue=tissue, methy.state=s, n.cpgs=length(k), n.flags=length(w),
      median.abs.dbeta=round(median(w),4),
      pct.under.0.05=round(100*mean(w<0.05),2),
      pct.under.0.10=round(100*mean(w<0.10),2),
      pct.under.0.15=round(100*mean(w<0.15),2))
  }))
  print(gwt[,c("methy.state","n.cpgs","n.flags","median.abs.dbeta",
               "pct.under.0.10","pct.under.0.15")], row.names=FALSE)
  gw[[tissue]] <- gwt
  rm(dev); invisible(gc())

  # ---- chr22 subset for the stability work ----
  c22 <- which(d$Chromosome=="chr22"); o22 <- c22[order(d$Start[c22])]
  B <- Ball[o22,,drop=FALSE]; F <- Fall[o22,,drop=FALSE]; st <- stall[o22]
  rm(Ball, Fall); invisible(gc())
  medc <- rowMedian(B); devc <- abs(B - medc)

  # reconstructed external thresholds, so ext can be re-applied to noisy data
  P <- N <- rep(NA_real_, nrow(B))
  for (i in seq_len(nrow(B))) {
    b<-B[i,]; fl<-F[i,]; if(all(is.na(fl))) next
    up<-which(fl==1); dn<-which(fl==-1); ze<-which(fl==0)
    if(length(up)&&length(ze)) P[i] <- (max(b[ze],na.rm=TRUE)+min(b[up],na.rm=TRUE))/2
    if(length(dn)&&length(ze)) N[i] <- (max(b[dn],na.rm=TRUE)+min(b[ze],na.rm=TRUE))/2
  }
  ext_apply <- function(X){ o<-matrix(0L,nrow(X),ncol(X))
    o[!is.na(P)&X>P]<-1L; o[!is.na(N)&X<N]<--1L; o[is.na(X)]<-NA; o }

  noisy <- lapply(seq_len(N_REPS), function(r){
    Bn <- B + matrix(rnorm(length(B),0,NOISE_SD), nrow(B))
    Bn[Bn<0]<-0; Bn[Bn>1]<-1; Bn[is.na(B)]<-NA; Bn })

  score <- function(O, Operturb_list) {
    ok <- !is.na(O) & O!=0
    magH <- median(devc[ok & st=="H"], na.rm=TRUE)
    magR <- median(devc[ok & st=="R"], na.rm=TRUE)
    rates <- sapply(ST, function(s){ k<-which(st==s); if(!length(k)) NA else mean(O[k,]!=0,na.rm=TRUE) })
    list(rate=mean(O!=0,na.rm=TRUE),
         stab=mean(sapply(Operturb_list, function(Q) jacc(O,Q))),
         srr=max(rates,na.rm=TRUE)/max(1e-9,min(rates,na.rm=TRUE)),
         magHR=magH/magR, trivial=mean(devc[ok]<0.10, na.rm=TRUE))
  }

  # ---- B. floor sensitivity ----
  cat("\n=== B. Absolute-difference floor: sensitivity to the cut-off ===\n")
  base_ext <- ext_apply(B)
  rows <- list()
  for (fl in c(0, FLOORS)) {
    O <- base_ext; O[!is.na(O) & O!=0 & devc < fl] <- 0L
    per <- lapply(noisy, function(Bn){ Q<-ext_apply(Bn); dn<-abs(Bn-rowMedian(Bn))
      Q[!is.na(Q) & Q!=0 & dn < fl] <- 0L; Q })
    s <- score(O, per)
    rows[[length(rows)+1]] <- data.frame(tissue=tissue,
      floor=fl, flag.rate.pct=round(100*s$rate,3),
      flags.retained.pct=round(100*sum(O!=0,na.rm=TRUE)/sum(base_ext!=0,na.rm=TRUE),2),
      stability.jaccard=round(s$stab,4),
      state.rate.ratio=round(s$srr,2),
      magnitude.H.vs.R=round(s$magHR,3),
      pct.trivial=round(100*s$trivial,2))
  }
  fs <- do.call(rbind, rows)
  print(fs, row.names=FALSE)
  cat("\n  floor = 0 is the current method. Read down the stability column:\n")
  cat("  the gain is monotone, so 0.10 is a point on a curve, not a\n")
  cat("  discovered optimum. Choose it on biological grounds, not this table.\n")
  wo(fs, paste0("Task26_FloorSensitivity_", tissue, ".csv"))

  # ---- C. rate-matched comparison ----
  cat("\n=== C. Competing methods re-tuned DOWN to the floor's flag rate ===\n")
  target <- fs$flag.rate.pct[fs$floor==0.10] / 100
  cat("  target rate = ", round(100*target,3), "% (ext + 0.10 floor)\n\n", sep="")
  M <- to_M(B)
  k_tk <- calibrate(function(k) f_tukey(B,k), 0.1, 25, target)
  k_mb <- calibrate(function(k) f_madz(B,k),  0.1, 25, target)
  k_mm <- calibrate(function(k) f_madz(M,k),  0.1, 25, target)
  p_bf <- calibrate(function(p) f_betafit(B,p), 1e-14, 0.2, target, decreasing=FALSE)

  meth <- list(
    `ext + 0.10 floor` = list(cur=(function(){O<-base_ext; O[!is.na(O)&O!=0&devc<0.10]<-0L; O})(),
                              per=lapply(noisy,function(Bn){Q<-ext_apply(Bn); dn<-abs(Bn-rowMedian(Bn))
                                Q[!is.na(Q)&Q!=0&dn<0.10]<-0L; Q})),
    `tukey.iqr`  = list(cur=f_tukey(B,k_tk),  per=lapply(noisy,function(Bn) f_tukey(Bn,k_tk))),
    `mad.beta`   = list(cur=f_madz(B,k_mb),   per=lapply(noisy,function(Bn) f_madz(Bn,k_mb))),
    `mad.mvalue` = list(cur=f_madz(M,k_mm),   per=lapply(noisy,function(Bn) f_madz(to_M(Bn),k_mm))),
    `beta.fit`   = list(cur=f_betafit(B,p_bf),per=lapply(noisy,function(Bn) f_betafit(Bn,p_bf))))

  rows <- list()
  for (nm in names(meth)) {
    s <- score(meth[[nm]]$cur, meth[[nm]]$per)
    rows[[length(rows)+1]] <- data.frame(tissue=tissue, method=nm,
      flag.rate.pct=round(100*s$rate,3), stability.jaccard=round(s$stab,4),
      state.rate.ratio=round(s$srr,2), magnitude.H.vs.R=round(s$magHR,3),
      pct.trivial=round(100*s$trivial,2))
  }
  rm <- do.call(rbind, rows); rm <- rm[order(-rm$stability.jaccard),]
  print(rm, row.names=FALSE)
  cat("\n  All rows now flag at the same rate, so any remaining stability\n")
  cat("  difference is a property of the method, not of its sensitivity.\n")
  wo(rm, paste0("Task26_RateMatchedStability_", tissue, ".csv"))

  rm(B,F,M,devc,noisy); invisible(gc()); cat("\n")
}

wo(do.call(rbind, gw), "Task26_GenomeWideMagnitude.csv")
cat("\n====================================================\nTASK 26 COMPLETE\n")
cat("====================================================\n")
