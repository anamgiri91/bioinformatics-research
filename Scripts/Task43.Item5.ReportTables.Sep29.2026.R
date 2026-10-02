# Reproducible report summaries; no method fitting or threshold selection.
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- dirname(dirname(normalizePath(sub("^--file=", "", file_arg[1L]))))
suppressPackageStartupMessages(library(data.table))
source(file.path(root,"Scripts","Task41.RRBS.Item5.Feasibility.Sep28.2026.R"))
out <- file.path(root,"Results")
p <- fread(file.path(out,"Task41_PairMetrics.csv.gz"))
p <- p[pair_set=="true_consecutive" & n_common==18]
st <- fread(file.path(out,"Task39_RRBS_CompleteSiteStates.csv.gz"),select=c("pos","state"))
p[, s1:=st$state[match(pos1,st$pos)]][, s2:=st$state[match(pos2,st$pos)]]
stopifnot(nrow(p)==32664L,!anyNA(p$s1),!anyNA(p$s2))
p[, state_pair:=fifelse(s1==s2,paste0(s1,"-",s2),"different")]
tab <- p[, .(pairs=.N,median_exact=median(agreement_exact),median_ordinal=median(ordinal_similarity),
  pearson_defined=sum(is.finite(pearson_raw)),median_pearson=median(pearson_raw,na.rm=TRUE),
  xi_defined=sum(is.finite(xi_mean)),median_xi=median(xi_mean,na.rm=TRUE)),by=.(state_pair,gap_bin)]
fwrite(tab,file.path(out,"Task43_TrueAdjacentStateByGap.csv"))
quantization <- p[, .(pairs=.N,raw_constant_pairs=sum(constant_raw_any),
  binned_constant_pairs=sum(constant_levels_any),
  variable_raw_but_constant_after_binning=sum(!constant_raw_any & constant_levels_any),
  median_abs_pearson_in_lost_pairs=median(abs(pearson_raw[!constant_raw_any & constant_levels_any]),na.rm=TRUE))]
fwrite(quantization,file.path(out,"Task43_QuantizationAudit.csv"))

# A directly usable score table: no combined score or invented probability.
# Support counts describe sampling limitations, not an exclusion/tuning rule.
d <- fread(file.path(root,"Data","NN.hg38.18P.forw.chr22.w.header.txt"))
levels <- levels_beta(as.matrix(d[,3:20]))
x <- levels[p$row1,,drop=FALSE]; y <- levels[p$row2,,drop=FALSE]
nonmodal <- function(z) ncol(z)-apply(z,1,function(v) max(tabulate(v,nbins=10)))
card <- p[, .(pos1,pos2,gap_bp,n_common,s1,s2,
  ordinal_similarity,raw_beta_similarity=1-beta_MAE,beta_MAE,beta_RMSE,
  agreement_exact,agreement_within_one,pearson_raw,spearman_raw,spearman_levels,
  xi_xy,xi_yx,xi_mean,constant_raw_any,constant_levels_any)]
card[, min_nonmodal_level_samples:=pmin(nonmodal(x),nonmodal(y))]
card[, support_status:=fcase(min_nonmodal_level_samples==0L,"constant_level_profile",
                             min_nonmodal_level_samples==1L,"singleton_level_profile",
                             default="both_have_multiple_nonmodal_samples")]
fwrite(card,file.path(out,"Task43_PrimaryScorecard.csv.gz"))
fwrite(card[,.(pairs=.N),by=support_status],file.path(out,"Task43_ScoreSupport.csv"))

s <- fread(file.path(out,"Task43_BenchmarkSummary.csv"))
stopifnot(nrow(s)==99L,all(s$generated==400L))
fwrite(dcast(s,scenario+truth~method,value.var="operational_rejection_rate"),
  file.path(out,"Task43_BenchmarkComparison.csv"))
# Paired comparison: differences estimated on exactly the same replicates.
r <- fread(file.path(out,"Task43_BenchmarkReplicates.csv"))
comparisons <- list(c("xi_levels_mean","dcor_levels"),c("pearson_beta","pearson_levels"),
                    c("dcor_beta","dcor_levels"),c("spearman_beta","spearman_levels"))
comparison <- rbindlist(lapply(comparisons,function(methods) {
  a <- r[method==methods[1]]; b <- r[method==methods[2]]
  z <- merge(a,b,by=c("scenario","truth","replicate"),suffixes=c("_a","_b"))
  z[, .(method_a=methods[1],method_b=methods[2],replicates=.N,
        a_only=sum(rejected_a & !rejected_b),b_only=sum(!rejected_a & rejected_b),
        difference_a_minus_b=mean(as.numeric(rejected_a)-as.numeric(rejected_b))),by=.(scenario,truth)]
}))
fwrite(comparison,file.path(out,"Task43_PairedMethodDifferences.csv"))
cat("PASS: primary real-data population and 99 benchmark summary rows verified\n")
print(quantization)
