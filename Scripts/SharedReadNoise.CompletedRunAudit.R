# Post-result audit of Tasks 59 and 65. Run from the repository root.
# Reads frozen outputs without rebuilding counts, splitting reads, or changing
# the estimators. New outputs are descriptive, not additional locked tests.
suppressPackageStartupMessages({library(data.table); library(jsonlite); library(digest); library(ggplot2)})
stopifnot(file.exists("Results/Task64_LockedProtocol.json"))
out <- "Results/SharedReadNoise_CompletedRunAudit"
spec <- data.table(task=c(59L,65L), lock=c(58L,64L),
                   cohort=c("Blood RRBS", "Plasma capture"),
                   rds=c("Data/gtex_wbc_rrbs/srb_pair_results.rds",
                         "Data/gse149438_plasma_mhap/s6_pair_results.rds"))
locks <- list(); inputs <- character(); summaries <- list(); distances <- list()
close <- function(x,y) stopifnot(isTRUE(all.equal(x,y,tolerance=1e-12)))
gap_breaks <- c(0,10,20,40,60,100,150,200)
for (k in seq_len(nrow(spec))) {
  z <- spec[k]
  lf <- sprintf("Results/Task%d_LockedProtocol.json",z$lock)
  lock <- fromJSON(lf)
  actual <- vapply(names(lock$sha256),function(f) digest(file=f,algo="sha256"),character(1))
  stopifnot(all(actual == unlist(lock$sha256)))
  locks[[k]] <- list(task=z$lock, files=length(actual), all_match=TRUE)
  pf <- sprintf("Results/Task%d_Primary.csv",z$task)
  p <- fread(pf); x <- readRDS(z$rds); e <- as.data.table(x$est)
  stopifnot(nrow(e)==p$pairs, all(e$n>=20), all(e$gap>0 & e$gap<=200),
            !anyDuplicated(e[,.(chrn,idx)]))
  raw_loss <- (e$raw_cov_A-e$ref_cov)^2
  cor_loss <- (e$cor_cov_A-e$ref_cov)^2
  delta <- cor_loss-raw_loss
  close(raw_loss,e$loss_raw); close(cor_loss,e$loss_cor); close(delta,e$dC)
  close(e$raw_cov_A-e$mean_noise_A,e$cor_cov_A)
  close(signif(mean(delta),4),p$delta_mse_cov)
  close(round(mean(delta)/mean(raw_loss),4),p$relative_change)
  close(round(mean(cor_loss<raw_loss),4),p$share_cor_closer)
  close(signif(mean(e$dD),4),p$delta_mse_d2)
  boot <- e$idx %% 25L == 0L
  saved_n <- if ("bootstrap_pairs" %in% names(p)) p$bootstrap_pairs else p$subsample_pairs
  stopifnot(sum(boot)==saved_n)
  close(signif(mean(delta[boot]),4),p$subsample_estimate)
  # Saved intervals are recorded, not recomputed by this audit.
  summaries[[k]] <- data.table(cohort=z$cohort, task=z$task, donors=p$donors,
    pairs=nrow(e), bootstrap_pairs=sum(boot), mean_raw=mean(e$raw_cov_A),
    mean_corrected=mean(e$cor_cov_A), mean_reference=mean(e$ref_cov),
    mean_noise=mean(e$mean_noise_A), delta_mse=mean(delta),
    relative_change=mean(delta)/mean(raw_loss),
    share_better=mean(cor_loss<raw_loss), share_equal=mean(cor_loss==raw_loss),
    saved_lo95=p$lo95,saved_hi95=p$hi95)
  e[, bin := cut(gap,gap_breaks)]
  distances[[k]] <- e[,.(pairs=.N, raw=mean(raw_cov_A),corrected=mean(cor_cov_A),
    reference=mean(ref_cov),noise=mean(mean_noise_A),
    pairs_with_represented_overlap=sum(overlap_A>0),
    pairs_with_nonzero_correction=sum(mean_noise_A!=0)),by=bin][,cohort:=z$cohort]
  if (z$task==65L) {
    close(signif(mean(delta[e$gap<=40]),4),p$near_delta_mse_cov)
    stopifnot(sum(e$gap<=40)==p$near_pairs,
              all(e$overlap_A[e$gap>150]==0),all(e$mean_noise_A[e$gap>150]==0))
  }
  inputs <- c(inputs,lf,pf,z$rds)
  cat("PASS Task",z$task,":",nrow(e),"pairs;",sum(boot),"bootstrap pairs; lock intact\n")
  rm(x,e); invisible(gc())
}
s <- rbindlist(summaries); d <- rbindlist(distances)
fwrite(s,paste0(out,"_Summary.csv")); fwrite(d,paste0(out,"_Distance.csv"))
plotdata <- melt(d,id.vars=c("cohort","bin","pairs"),
                 measure.vars=c("raw","corrected","reference"),
                 variable.name="estimate",value.name="covariance")
plotdata[,estimate:=factor(estimate,levels=c("raw","corrected","reference"),
                          labels=c("Observed A","Corrected A","Cross-partition B/C"))]
fig <- ggplot(plotdata,aes(bin,covariance,colour=estimate,group=estimate,shape=estimate))+
  geom_line(linewidth=.6)+geom_point(size=2)+facet_wrap(~cohort,ncol=2)+
  scale_colour_manual(values=c("#cc6633","#2369a2","#20805a"))+
  labs(title="Mean CpG covariance by genomic distance",
       subtitle="Descriptive summaries of completed tests; estimates use the same eligible pairs within each bin",
       x="Distance between consecutive reference CpGs (bp)",y="Mean covariance",colour=NULL,shape=NULL,
       caption=paste("No uncertainty intervals shown; pair composition and depth differ between bins.",
         "Plasma B/C may retain shared-fragment noise because some mates remain separate records.",sep="\n"))+
  theme_minimal(base_size=11)+theme(legend.position="bottom",axis.text.x=element_text(angle=40,hjust=1),
                                  plot.caption=element_text(hjust=0))
ggsave(paste0(out,"_Distance.png"),fig,width=12,height=5.3,dpi=160,bg="white")
inputs <- c(inputs,"Scripts/SharedReadNoise.CompletedRunAudit.R")
write_json(list(audited_at=format(Sys.time(),tz="UTC",usetz=TRUE),locks=locks,
  scope="Output arithmetic and current lock integrity; not a raw-read reanalysis or bootstrap rerun",
  input_sha256=as.list(setNames(vapply(inputs,function(f) digest(file=f,algo="sha256"),character(1)),inputs))),
  paste0(out,"_Manifest.json"),pretty=TRUE,auto_unbox=TRUE)
cat("Audit complete; frozen inputs and outputs were not modified.\n")
