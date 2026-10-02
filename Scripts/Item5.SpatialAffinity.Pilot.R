# Exploratory item-5 extension: Manhattan agreement + distance correlation + bp gap.
# This is a composite affinity, not a new correlation or calibrated probability.
# Run: Rscript Scripts/Item5.SpatialAffinity.Pilot.R
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script <- normalizePath(sub("^--file=", "", file_arg[1L]))
root <- dirname(dirname(script))
source(file.path(root, "Scripts", "Task41.RRBS.Item5.Feasibility.Sep28.2026.R"))
out <- file.path(root, "Results")
write_tab <- function(x, name) fwrite(x, file.path(out, paste0("Item5_Spatial_", name, ".csv")))

# Biased sample dCor with Euclidean distances for univariate beta profiles.
# Matrix form is checked below against an independently expanded scalar formula.
dcor_pair <- function(x, y) {
  if (length(unique(x)) < 2L || length(unique(y)) < 2L) return(NA_real_)
  center <- function(z) {
    z <- abs(outer(z, z, "-"))
    z - rowMeans(z) - rep(colMeans(z), each = nrow(z)) + mean(z)
  }
  a <- center(x); b <- center(y)
  pmin(1, sqrt(pmax(sum(a*b), 0)/sqrt(sum(a*a)*sum(b*b))))
}
dcor_rows <- function(x, y) vapply(seq_len(nrow(x)), function(i) dcor_pair(x[i,],y[i,]), 0.0)
affinity <- function(x, y, gap, ell = 100) {
  exp(-gap/ell) * ((1-mean(abs(x-y))) + dcor_pair(x,y))/2
}

set.seed(202609291)
for (n in c(6L,17L,18L)) {
  x <- runif(n); y <- runif(n)
  a <- abs(outer(x,x,"-")); b <- abs(outer(y,y,"-"))
  dcov <- mean(a*b) + mean(a)*mean(b) - 2*mean(rowMeans(a)*rowMeans(b))
  va <- mean(a*a) + mean(a)^2 - 2*mean(rowMeans(a)^2)
  vb <- mean(b*b) + mean(b)^2 - 2*mean(rowMeans(b)^2)
  stopifnot(abs(dcor_pair(x,y)-sqrt(pmax(dcov,0)/sqrt(va*vb))) < 1e-10,
            abs(dcor_pair(x,y)-dcor_pair(y,x)) < 1e-12,
            abs(dcor_pair(x,x)-1) < 1e-12,
            is.na(dcor_pair(rep(.1,n),y)))
}

# Freeze choices before reading the new pilot outcomes. No tuning or optimality claim.
jsonlite::write_json(list(date="2026-09-29", seed=202609291L,
  question="Does a fixed contextual affinity improve held-out absolute agreement ranking?",
  formula="A=exp(-gap/ell)*(S_raw+dCor_raw)/2; S_raw=1-mean(abs(beta_i-beta_j))",
  primary_ell_bp=100, sensitivity_ell_bp=c(50,200), quality_weight=1,
  alpha=.5, selection_fraction=.1, folds="Each of 18 patients held out once",
  target="Absolute difference of the two CpGs in the held-out patient; lower is better",
  population="Original-row consecutive CpG pairs complete in all 18 patients",
  common_set="Within each fold, both training profiles must vary so dCor is defined",
  baselines=c("gap_only","manhattan_raw","dcor_raw","unweighted_mix"),
  tie_rule="Ascending genomic position, same deterministic rule for all methods",
  inference="Descriptive patient-fold summaries only; shared sites and folds are dependent",
  limits="Internal exploratory pilot on an already examined cohort; not external validation, biological truth, or a novel principle",
  script_md5=unname(tools::md5sum(script))),
  file.path(out,"Item5_Spatial_Protocol.json"), pretty=TRUE,auto_unbox=TRUE)

p <- fread(file.path(out,"Task41_PairMetrics.csv.gz"))[
  pair_set=="true_consecutive" & n_common==18]
setorder(p,pos1,pos2)
d <- fread(file.path(root,"Data","NN.hg38.18P.forw.chr22.w.header.txt"))
beta <- as.matrix(d[,3:20])
x <- beta[p$row1,,drop=FALSE]; y <- beta[p$row2,,drop=FALSE]
stopifnot(nrow(p)==32664L,all(p$row2==p$row1+1L),!anyNA(x),!anyNA(y),
          all(p$gap_bp==p$pos2-p$pos1),all(p$gap_bp>0))
cat("Computing raw-beta dependence for",nrow(p),"original-row pairs\n")
dc <- dcor_rows(x,y)
xi_xy <- vapply(seq_len(nrow(p)),function(i) xi_tie_average(x[i,],y[i,]),0.0)
xi_yx <- vapply(seq_len(nrow(p)),function(i) xi_tie_average(y[i,],x[i,]),0.0)
card <- p[,.(row1,row2,pos1,pos2,gap_bp,n_common,
  ordinal_similarity,raw_beta_similarity=1-beta_MAE,
  pearson_raw,spearman_raw,xi_levels_mean=xi_mean,
  constant_raw_any,constant_levels_any)]
card[, `:=`(dcor_raw=dc,xi_raw_xy=xi_xy,xi_raw_yx=xi_yx,
             xi_raw_mean=(xi_xy+xi_yx)/2)]
card[, unweighted_mix:=(raw_beta_similarity+dcor_raw)/2]
for (ell in c(50,100,200)) card[,paste0("affinity_",ell,"bp"):=exp(-gap_bp/ell)*unweighted_mix]
stopifnot(all(is.na(dc)==card$constant_raw_any))
fwrite(card,file.path(out,"Item5_Spatial_PairScores.csv.gz"))

# A mathematical check with ties: fixed spatial weighting preserves each pair's
# permutation ordering. Enumerate all 18 alignments of a singleton profile.
sx <- c(.95,rep(.05,17)); sy <- sx
stats <- vapply(seq_len(18),function(i) {
  v <- rep(.05,18); v[i] <- .95
  ((1-mean(abs(sx-v)))+dcor_pair(sx,v))/2
},0.0)
inv <- rbindlist(lapply(c(2,8,50,100,205),function(gap) {
  w <- exp(-gap/100)
  data.table(gap_bp=gap,observed=stats[1],weighted_observed=w*stats[1],
    exact_p=sum(stats >= stats[1]-1e-12)/18,
    weighted_exact_p=sum(w*stats >= w*stats[1]-w*1e-12)/18)
}))
stopifnot(all(inv$exact_p==inv$weighted_exact_p),all(inv$exact_p==1/18))
write_tab(inv,"PermutationInvariance")

# Common-set rank sensitivity. The primary 100-bp scale is not selected using this.
keep <- which(is.finite(card$unweighted_mix)); k <- ceiling(length(keep)*.1)
choose <- function(score,ids,k) ids[order(-score[ids],card$pos1[ids],card$pos2[ids])[seq_len(k)]]
ref <- choose(card$unweighted_mix,keep,k)
rank_audit <- rbindlist(lapply(c(50,100,200),function(ell) {
  v <- card[[paste0("affinity_",ell,"bp")]]
  selected <- choose(v,keep,k)
  data.table(ell_bp=ell,eligible=length(keep),top_pairs=k,
    overlap_with_unweighted_top=sum(selected %in% ref)/k,
    rank_correlation=cor(v[keep],card$unweighted_mix[keep],method="spearman"),
    median_gap_selected=median(card$gap_bp[selected]),
    median_gap_unweighted=median(card$gap_bp[ref]))
}))
write_tab(rank_audit,"RankSensitivity")

# No held-out beta enters score computation or within-fold eligibility. The
# candidate pair universe is complete-case, which still limits generalization.
fold_results <- list()
for (held in seq_len(ncol(x))) {
  tx <- x[,-held,drop=FALSE]; ty <- y[,-held,drop=FALSE]
  agreement <- 1-rowMeans(abs(tx-ty)); dependence <- dcor_rows(tx,ty)
  mix <- (agreement+dependence)/2
  eligible <- which(is.finite(dependence)); k <- ceiling(length(eligible)*.1)
  held_error <- abs(x[,held]-y[,held])
  scores <- list(gap_only=-p$gap_bp,manhattan_raw=agreement,dcor_raw=dependence,
    unweighted_mix=mix,spatial_50bp=exp(-p$gap_bp/50)*mix,
    spatial_100bp=exp(-p$gap_bp/100)*mix,spatial_200bp=exp(-p$gap_bp/200)*mix)
  fold_results[[held]] <- rbindlist(lapply(names(scores),function(method) {
    selected <- choose(scores[[method]],eligible,k)
    data.table(held_patient=held,method=method,total_pairs=nrow(p),
      eligible_pairs=length(eligible),selected_pairs=k,
      heldout_mean_absolute_difference=mean(held_error[selected]),
      heldout_median_absolute_difference=median(held_error[selected]),
      median_selected_gap_bp=median(p$gap_bp[selected]),
      all_eligible_mean_absolute_difference=mean(held_error[eligible]))
  }))
  cat("Completed held-patient fold",held,"of",ncol(x),"\n")
}
folds <- rbindlist(fold_results)
write_tab(folds,"HeldPatientFolds")
summ <- folds[,.(patient_folds=.N,min_eligible=min(eligible_pairs),max_eligible=max(eligible_pairs),
  mean_heldout_absolute_difference=mean(heldout_mean_absolute_difference),
  min_patient_mean_difference=min(heldout_mean_absolute_difference),
  max_patient_mean_difference=max(heldout_mean_absolute_difference),
  median_selected_gap_bp=median(median_selected_gap_bp)),by=method]
write_tab(summ,"HeldPatientSummary")
write_tab(card[,.(pairs=.N,dcor_defined=sum(is.finite(dcor_raw)),
  raw_xi_mean_defined=sum(is.finite(xi_raw_mean)),
  binned_xi_mean_defined=sum(is.finite(xi_levels_mean)),
  min_gap=min(gap_bp),median_gap=median(gap_bp),max_gap=max(gap_bp))],"Coverage")
capture.output(sessionInfo(),file=file.path(out,"Item5_Spatial_SessionInfo.txt"))
print(summ)
print(rank_audit)
cat("PASS: formula identity, symmetry, constants, true adjacency, and spatial permutation invariance\n")
