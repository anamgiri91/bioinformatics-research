# Read-only integrity and coordinate checks before running frozen sections 5/6.
suppressPackageStartupMessages({library(data.table); library(jsonlite); library(matrixStats)})
arg <- grep("^--file=",commandArgs(FALSE),value=TRUE)
root <- dirname(dirname(normalizePath(sub("^--file=","",arg[1L]))))
proto <- fromJSON(file.path(root,"Results","Task47_FrozenProtocol.json"))
for (f in names(proto$files_sha256))
  stopifnot(digest::digest(file=file.path(root,f),algo="sha256")==proto$files_sha256[[f]])
cc <- readRDS(file.path(root,"Data","gtex_colon_rrbs","cache_cohort.rds"))
stopifnot(identical(dim(cc$M),dim(cc$N)),ncol(cc$M)==29,!anyNA(cc$M),!anyNA(cc$N),
          all(cc$M>=0),all(cc$M<=cc$N),all(cc$N>=10),!anyDuplicated(cc$samples))
cpg <- readRDS(file.path(root,"Data","hg19_seq","cpg_positions_hg19.rds"))
allpos <- unlist(cpg,use.names=FALSE); nchr <- lengths(cpg); ends <- cumsum(nchr)
starts <- c(0L,head(ends,-1))
st <- as.data.table(cc$sites)
stopifnot(all(st$pos==allpos[st$idx]),all(st$idx>starts[st$chrn]),all(st$idx<=ends[st$chrn]))
samples <- fromJSON(file.path(root,"Data","gtex_colon_rrbs","colon_samples.json"))
qc <- fread(file.path(root,"Results","Task50_BuildSummary.csv"))
selected <- unique(c(1L,which.min(qc$share_same_beta_as_geo)))
checks <- list()
for (s in selected) {
  bed <- fread(file.path(root,"Data","gtex_colon_rrbs",basename(samples$files[[s]][1])),
    header=FALSE,col.names=c("chr","pos","end","beta"))
  bed <- bed[chr %in% paste0("chr",1:22)][,chrn:=match(chr,paste0("chr",1:22))]
  stopifnot(!anyDuplicated(bed[,.(chrn,pos)]),all(bed$beta>=0 & bed$beta<=1))
  setkey(bed,chrn,pos)
  beta_pat <- cc$M[,s]/cc$N[,s]
  for (offset in c(-1L,0L,1L,100L)) {
    ix <- st$idx+offset
    ok <- ix>starts[st$chrn] & ix<=ends[st$chrn]
    ref <- rep(NA_real_,nrow(st))
    ref[ok] <- bed[data.table(chrn=st$chrn[ok],pos=allpos[ix[ok]]),on=.(chrn,pos),beta]
    v <- data.table(chrn=st$chrn,pat=beta_pat,reference=ref)[is.finite(reference)]
    checks[[length(checks)+1L]] <- v[,{
      mid <- reference>.05 & reference<.95
      list(sites=.N,exact_share=mean(abs(pat-reference)<1e-6),mean_absolute_difference=mean(abs(pat-reference)),
        correlation=cor(pat,reference),intermediate_sites=sum(mid),
        intermediate_exact_share=mean(abs(pat[mid]-reference[mid])<1e-6),
        intermediate_mae=mean(abs(pat[mid]-reference[mid])))
    },by=chrn][,`:=`(sample=cc$samples[s],cpg_index_offset=offset)]
    cat("Validated coordinate comparison:",cc$samples[s],"offset",offset,"\n")
  }
}
fwrite(rbindlist(checks),file.path(root,"Results","GTExColon_MappingAudit.csv"))

# Test the intended straight-line endpoint when a training predictor is constant.
Xt <- matrix(rep(1/3,26),1); Yt <- matrix(seq(.1,.9,length.out=26),1)
mx <- rowMeans(Xt); my <- rowMeans(Yt)
vx <- rowMeans((Xt-mx)^2); vy <- rowMeans((Yt-my)^2)
cxy <- rowMeans((Xt-mx)*(Yt-my))
legacy_slope <- ifelse(vx>0,cxy/vx,0)
constant_x <- rowMins(Xt)==rowMaxs(Xt)
vx[constant_x] <- 0; cxy[constant_x] <- 0
corrected_slope <- ifelse(vx>0,cxy/vx,0)
stopifnot(corrected_slope==0,abs(legacy_slope)>1)
fwrite(data.table(case="constant one-third predictor; 26 training samples",legacy_slope=legacy_slope,
                 corrected_slope=corrected_slope,expected_intercept_only_slope=0),
       file.path(root,"Results","GTExColon_ConstantPredictorCheck.csv"))
cat("PASS: frozen hashes, counts, genome-index coordinates and constant-predictor check\n")
