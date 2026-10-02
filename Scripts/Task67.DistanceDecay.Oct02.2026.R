# Post-result development on already examined colon and blood. No plasma redo.
# Run from repo root. Main estimates and intervals use exactly the same pairs.
suppressPackageStartupMessages({library(data.table);library(jsonlite);library(digest);library(ggplot2)})
source("Scripts/SharedReadNoise.Decay.R")
setDTthreads(4L)
cpg <- readRDS("Data/hg19_seq/cpg_positions_hg19.rds")
offset <- c(0L,cumsum(lengths(cpg))[-22L])
spec <- data.table(cohort=c("Colon (development)","Blood (post-result)"),
  file=c("Data/gtex_colon_rrbs/srb_pair_results.rds","Data/gtex_wbc_rrbs/srb_pair_results.rds"))
BREAKS <- c(0,10,20,40,60,100,150,200)
all_results <- list(); supports <- list(); curves <- list()
for (k in seq_len(nrow(spec))) {
  cat("Reading",spec$cohort[k],"\n")
  x <- readRDS(spec$file[k]); e <- copy(as.data.table(x$est))
  e[, `:=`(pos=0L,gap=0L,density=0L)]
  for (ch in unique(e$chrn)) {
    ii <- which(e$chrn==ch); local <- e$idx[ii]-offset[ch]; pp <- cpg[[ch]]
    stopifnot(all(local>=1 & local<length(pp)))
    pos <- pp[local]; gap <- pp[local+1L]-pos; mid <- floor(pos+gap/2)
    dens <- findInterval(mid+500L,pp)-findInterval(mid-501L,pp)-2L
    stopifnot(all(dens>=0))
    set(e,ii,"pos",pos);set(e,ii,"gap",gap);set(e,ii,"density",dens)
  }
  e[, bin := cut(gap,BREAKS)]
  e[, cell_curve := paste(cut(depth_A,c(2,5,10,20,Inf),right=FALSE),
                         cut(density,c(0,10,30,Inf),right=FALSE),sep="|")]
  cc <- e[, .N,by=.(cell_curve,bin)]
  good <- cc[,.(bins=.N,smallest=min(N)),by=cell_curve][bins==7 & smallest>=20,cell_curve]
  # All seven bins share this same common-support distribution of count/sequence cells.
  target <- cc[cell_curve %in% good,.(mass=min(N)),by=cell_curve][,target:=mass/sum(mass)]
  ef <- merge(e[cell_curve %in% good],target,by="cell_curve",sort=FALSE)
  ef[, weight_curve := target/.N,by=.(cell_curve,bin)]
  rawcurve <- e[,.(pairs=.N,raw=mean(raw_cov_A),corrected=mean(cor_cov_A),reference=mean(ref_cov)),by=bin][,curve:="Unadjusted"]
  matched <- ef[,.(pairs=.N,raw=sum(weight_curve*raw_cov_A),corrected=sum(weight_curve*cor_cov_A),reference=sum(weight_curve*ref_cov)),by=bin][,curve:="Depth + CpG density"]
  curves[[k]] <- rbind(rawcurve,matched)[,cohort:=spec$cohort[k]]
  # Use existing modulo-25 cache, without pairing a full-data point with a subset CI.
  at <- match(paste(x$sub$chrn,x$sub$idx),paste(e$chrn,e$idx)); ok <- which(!is.na(at))
  meta <- e[at[ok]]; design <- decay_design(meta,min_per_band=2L,within_block=TRUE)
  ii <- ok[design$row_id]
  pre <- decay_prepare(lapply(x$subm,function(M) M[ii,,drop=FALSE]))
  v <- decay_values(pre)
  expected <- as.matrix(meta[design$row_id,.(raw_cov_A,cor_cov_A,ref_cov)])
  stopifnot(max(abs(v-expected))<1e-10)
  est <- decay_metrics(v,design)
  cat("Matched",nrow(design),"pairs in",uniqueN(design$block),"blocks\n")
  boot <- decay_bootstrap(pre,design,draws=2000L,seed=2026100200L+k)
  bad <- rowSums(!is.finite(boot))>0
  # Never discard failed draws and silently return an interval.
  ci <- if(any(bad)) matrix(NA_real_,2,7) else apply(boot,2,quantile,c(.025,.975))
  centered <- sweep(ci,2,colMeans(boot)-est,"-")
  all_results[[k]] <- data.table(cohort=spec$cohort[k],metric=names(est),estimate=unname(est),
    exploratory_lo95=centered[1,],exploratory_hi95=centered[2,],
    plain_lo95=ci[1,],plain_hi95=ci[2,],bootstrap_draws=nrow(boot),failed_draws=sum(bad),
    status="post-result development; interval calibration pending")
  supports[[k]] <- data.table(cohort=spec$cohort[k],full_pairs=nrow(e),cached_eligible_pairs=nrow(meta),
    available_near=sum(meta$gap<=10),available_far=sum(meta$gap>150),
    matched_near=sum(design$band=="near"),matched_far=sum(design$band=="far"),
    cells=uniqueN(design$cell),blocks=uniqueN(design$block),chromosomes=uniqueN(design$chrn))
  saveRDS(list(design=design,point=est,bootstrap=boot),sprintf("Data/Task67_%d_development.rds",k))
  rm(x,e,ef,pre,v,expected);invisible(gc())
}
ans <- rbindlist(all_results); support <- rbindlist(supports); curve <- rbindlist(curves)
fwrite(ans,"Results/Task67_DecayContrasts.csv");fwrite(support,"Results/Task67_CommonSupport.csv")
fwrite(curve,"Results/Task67_DistanceCurves.csv")
long <- melt(curve,id.vars=c("cohort","curve","bin","pairs"),measure.vars=c("raw","corrected","reference"))
p <- ggplot(long,aes(bin,value,colour=variable,group=variable))+geom_line()+geom_point(size=1.5)+
  facet_grid(curve~cohort)+scale_colour_manual(values=c(raw="#c75c31",corrected="#286da8",reference="#288254"))+
  labs(title="Co-methylation distance profiles after depth and CpG-density adjustment",
       subtitle="Post-result development; absolute covariance, not a causal biological decay rate",
       x="Consecutive CpG distance (bp)",y="Mean covariance",colour=NULL,
       caption="Primary near–far contrast additionally matches within 1-Mb blocks; see Task67_DecayContrasts.csv.\nEach adjusted curve uses common cells across all seven bins; no outcome values determine weights.")+
  theme_minimal(base_size=11)+theme(legend.position="bottom",axis.text.x=element_text(angle=40,hjust=1))
ggsave("Results/Fig48_AdjustedDistanceDecay.png",p,width=11,height=7,dpi=160,bg="white")
inputs <- c(spec$file,"Data/hg19_seq/cpg_positions_hg19.rds","Scripts/SharedReadNoise.Decay.R", "Scripts/Task67.DistanceDecay.Oct02.2026.R")
write_json(list(status="post-result development, not a confirmatory lock",created=format(Sys.time(),tz="UTC"),
  input_sha256=as.list(setNames(vapply(inputs,function(f) digest(file=f,algo="sha256"),character(1)),inputs))),
  "Results/Task67_RunManifest.json",pretty=TRUE,auto_unbox=TRUE)
print(ans);print(support);cat("done\n")
