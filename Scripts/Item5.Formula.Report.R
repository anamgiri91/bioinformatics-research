# Descriptive summaries and a standalone figure for the tested candidate.
arg <- grep("^--file=",commandArgs(FALSE),value=TRUE)
script <- normalizePath(sub("^--file=","",arg[1L]))
root <- dirname(dirname(script)); out <- file.path(root,"Results")
suppressPackageStartupMessages(library(data.table))
source(file.path(root,"Scripts","Item5.ContextGatedSimilarity.R"))
p <- fread(file.path(out,"Item5_Formula_PairScores.csv.gz"))
prior <- fread(file.path(out,"Item5_Spatial_PairScores.csv.gz"))
stopifnot(identical(p$pos1,prior$pos1),identical(p$pos2,prior$pos2))
p[,`:=`(pearson_raw=prior$pearson_raw,spearman_raw=prior$spearman_raw)]
fwrite(p,file.path(out,"Item5_Formula_Scorecard.csv.gz"))

f <- fread(file.path(out,"Item5_Formula_HeldPatientFolds.csv"))
wide <- dcast(f,held_patient~method,value.var="heldout_mean_absolute_difference")
comparison <- data.table(comparator=c("manhattan_raw","old_spatial_blend_100bp","formula_100bp_uncorrected"))
comparison[,`:=`(folds=18L,
  candidate_better_folds=vapply(comparator,function(m)sum(wide$formula_100bp<wide[[m]]),0L),
  mean_difference=vapply(comparator,function(m)mean(wide$formula_100bp-wide[[m]]),0.0))]
fwrite(comparison,file.path(out,"Item5_Formula_PairedComparisons.csv"))

# A small explicit counterexample: high agreement + inverse dependence.
x <- seq(.49,.51,length.out=18); y <- 1-x
inverse <- dependence_credit(x,y)
sx <- c(.95,rep(.05,17)); sy <- sx
singleton <- dependence_credit(sx,sy)
issues <- data.table(example=c("near_constant_inverse","shared_singleton"),gap_bp=c(500,500),
  agreement=c(1-mean(abs(x-y)),1),pearson=c(cor(x,y),1),
  dependence_credit=c(inverse["credit"],singleton["credit"]))
issues[,formula:=gated_affinity(agreement,dependence_credit,gap_bp)]
issues[,exact_singleton_alignment_p:=c(NA_real_,1/18)]
fwrite(issues,file.path(out,"Item5_Formula_Counterexamples.csv"))

s <- fread(file.path(out,"Item5_Formula_HeldPatientSummary.csv"))
methods <- c("manhattan_raw","formula_100bp","old_spatial_blend_100bp")
values <- s$mean_heldout_absolute_difference[match(methods,s$method)]
png(file.path(out,"Item5_Formula.png"),width=1600,height=1500,res=180)
par(mfrow=c(2,1),mar=c(4.3,4.8,3.6,1.5),mgp=c(2.9,.7,0),family="sans")
gap <- seq(0,500,length.out=501)
credits <- c(0,.25,.5,.75,1)
colors <- c("#7F8C8D","#56B4E9","#0072B2","#D55E00","#009E73")
plot(NA,xlim=c(0,550),ylim=c(0,1),xlab="Genomic separation (bp)",
  ylab="Candidate affinity",main="Distance changes the contribution of dependence",las=1)
mtext("Illustration: Manhattan similarity = 0.9; length scale = 100 bp",side=3,line=.4,cex=.85)
for(i in seq_along(credits)) {
  z <- gated_affinity(rep(.9,length(gap)),rep(credits[i],length(gap)),gap)
  lines(gap,z,col=colors[i],lwd=2.5)
  text(505,tail(z,1),labels=paste0("E=",credits[i]),col=colors[i],adj=0,cex=.75)
}
abline(v=100,lty=3,col="#999999")
par(mar=c(5,4.8,3.6,1.5))
bars <- barplot(values,names.arg=c("Manhattan alone","Candidate formula","Earlier fixed blend"),
  col=c("#0072B2","#009E73","#7F8C8D"),ylim=c(0,.047),border=NA,
  ylab="Held-out mean absolute beta difference",las=1,
  main="Internal agreement-ranking comparison")
text(bars,values+.002,labels=sprintf("%.4f",values),cex=.95)
mtext("Top 10% of eligible pairs; mean over 18 held-patient folds; lower is better",side=3,line=.4,cex=.8)
mtext("Exploratory reuse of this cohort; not independent biological validation",side=1,line=3,cex=.8)
dev.off()
jsonlite::write_json(list(script_md5=unname(tools::md5sum(script)),
  source_script_md5=unname(tools::md5sum(file.path(root,"Scripts","Item5.ContextGatedSimilarity.R")))),
  file.path(out,"Item5_Formula_ReportManifest.json"),pretty=TRUE,auto_unbox=TRUE)
print(comparison)
print(issues)
