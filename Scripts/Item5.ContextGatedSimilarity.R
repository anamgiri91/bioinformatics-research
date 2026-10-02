# Candidate: Manhattan agreement supported by a spatial prior and dependence.
# F = S * [w + (1-w) E], w = exp(-gap/ell).
# E is positive excess dCor^2 above its exact conditional permutation mean.
# This is a proposed affinity adaptation, not a new dCor or a p-value.

centered_distances <- function(x) {
  a <- abs(outer(x,x,"-"))
  a-rowMeans(a)-rep(colMeans(a),each=length(x))+mean(a)
}

dependence_credit <- function(x,y) {
  stopifnot(length(x)==length(y),length(x)>=2L,
            all(is.finite(x)),all(is.finite(y)))
  if (length(unique(x))<2L || length(unique(y))<2L) {
    return(c(q=NA_real_,q0=NA_real_,excess=NA_real_,credit=0))
  }
  a <- centered_distances(x); b <- centered_distances(y)
  den <- sqrt(sum(a*a)*sum(b*b))
  q <- pmin(1,pmax(0,sum(a*b)/den))
  q0 <- sum(diag(a))*sum(diag(b))/((length(x)-1)*den)
  # n=2 always has q0=1: there is no discriminating permutation information.
  if (!is.finite(q0) || 1-q0<=1e-12) {
    return(c(q=q,q0=q0,excess=NA_real_,credit=0))
  }
  excess <- (q-q0)/(1-q0)
  c(q=q,q0=q0,excess=excess,credit=pmin(1,pmax(0,excess)))
}

profile_components <- function(x,y) {
  stopifnot(identical(dim(x),dim(y)),all(is.finite(x)),all(is.finite(y)),
            all(x>=0 & x<=1),all(y>=0 & y<=1))
  e <- t(vapply(seq_len(nrow(x)),function(i) dependence_credit(x[i,],y[i,]),numeric(4)))
  cbind(similarity=1-rowMeans(abs(x-y)),e)
}

gated_affinity <- function(similarity,credit,gap,ell=100) {
  stopifnot(all(is.finite(similarity)),all(is.finite(credit)),
            all(similarity>=0 & similarity<=1),all(credit>=0 & credit<=1),
            all(is.finite(gap) & gap>=0),length(ell)==1L,is.finite(ell),ell>0)
  w <- exp(-gap/ell)
  similarity*(w+(1-w)*credit)
}

all_permutations <- function(v) {
  if (length(v)==1L) return(matrix(v,nrow=1L))
  do.call(rbind,lapply(seq_along(v),function(i) cbind(v[i],all_permutations(v[-i]))))
}

check_formula <- function() {
  # Enumerate small tied and untied configurations to verify the analytic mean.
  cases <- list(list(x=c(.05,.2,.3,.6,.8,.95),y=c(.9,.1,.7,.5,.2,.4)),
                list(x=c(0,0,.2,.2,1,1),y=c(0,.1,.1,.5,.5,1)))
  perms <- all_permutations(1:6)
  diagnostics <- lapply(seq_along(cases),function(i) {
    x <- cases[[i]]$x; y <- cases[[i]]$y
    e <- dependence_credit(x,y)
    vals <- t(apply(perms,1,function(perm) dependence_credit(x,y[perm])))
    a <- abs(outer(x,x,"-")); b <- abs(outer(y,y,"-"))
    cov2 <- mean(a*b)+mean(a)*mean(b)-2*mean(rowMeans(a)*rowMeans(b))
    varx <- mean(a*a)+mean(a)^2-2*mean(rowMeans(a)^2)
    vary <- mean(b*b)+mean(b)^2-2*mean(rowMeans(b)^2)
    stopifnot(abs(mean(vals[,"q"])-e["q0"])<1e-12,
              abs(mean(vals[,"excess"]))<1e-12,
              abs(e["q"]-cov2/sqrt(varx*vary))<1e-12,
              max(abs(e-dependence_credit(y,x)))<1e-12)
    data.frame(case=i,permutations=nrow(perms),analytic_q0=e["q0"],
      enumerated_mean_q=mean(vals[,"q"]),signed_excess_mean=mean(vals[,"excess"]),
      positive_credit_mean=mean(vals[,"credit"]))
  })
  s <- c(0,.25,.5,.9,1); e <- c(0,.1,.5,.9,1)
  f0 <- gated_affinity(s,e,0); f1 <- gated_affinity(s,e,50)
  f2 <- gated_affinity(s,e,200)
  stopifnot(all(f0==s),all(f2<=f1),all(f1<=s),all(f2>=0),
            dependence_credit(rep(.1,18),rep(.1,18))["credit"]==0,
            dependence_credit(c(0,1),c(1,0))["credit"]==0,
            abs(dependence_credit(1:18,1:18)["credit"]-1)<1e-12)
  cat("PASS: exact permutation mean, signed centering, symmetry, constants, bounds, distance monotonicity\n")
  do.call(rbind,diagnostics)
}

run_formula <- function(root,script) {
  suppressPackageStartupMessages(library(data.table))
  input <- file.path(root,"Data","NN.hg38.18P.forw.chr22.w.header.txt")
  out <- file.path(root,"Results")
  write_tab <- function(x,name) fwrite(x,file.path(out,paste0("Item5_Formula_",name,".csv")))
  write_tab(check_formula(),"ExactNullChecks")
  jsonlite::write_json(list(date="2026-09-29",formula="F=S*(exp(-g/ell)+(1-exp(-g/ell))*E)",
    S="1-mean(abs(x-y))",q="ordinary sample dCor squared",
    q0="tr(A)*tr(B)/((n-1)*sqrt(sum(A^2)*sum(B^2)))",
    E="max(0,(q-q0)/(1-q0)); 0 credit if constant or q0=1; undefined statistic flagged",
    primary_ell_bp=100,sensitivity_ell_bp=c(50,200),seed=NULL,
    scope="All complete original-row adjacent pairs; held-patient comparisons on a common training-variable subset",
    endpoint="Top 10% ranked pairs: absolute beta difference in each held-out patient",
    ties="Ascending genomic position for every method",
    status="New formula evaluated on an already examined cohort; exploratory development, not independent validation",
    tuning="No parameter selected from this run; length scale copied from prior primary protocol",
    provenance=list(input_md5=unname(tools::md5sum(input)),script_md5=unname(tools::md5sum(script))),
    prior_art="Conditional null mean is the established RV expectation; see accompanying derivation/references",
    inference="No p-value, population-correlation, biological-validation or novelty claim"),
    file.path(out,"Item5_Formula_Protocol.json"),pretty=TRUE,auto_unbox=TRUE)

  p <- fread(file.path(out,"Item5_Spatial_PairScores.csv.gz"))
  setorder(p,pos1,pos2)
  d <- fread(input); beta <- as.matrix(d[,3:20])
  x <- beta[p$row1,,drop=FALSE]; y <- beta[p$row2,,drop=FALSE]
  stopifnot(nrow(p)==32664L,all(p$row2==p$row1+1L),!anyNA(x),!anyNA(y))
  cat("Computing candidate for",nrow(p),"pairs\n")
  z <- profile_components(x,y)
  stopifnot(max(abs(sqrt(z[,"q"])-p$dcor_raw),na.rm=TRUE)<1e-10,
            identical(is.na(z[,"q"]),is.na(p$dcor_raw)),
            max(abs(z[,"similarity"]-p$raw_beta_similarity))<1e-12)
  card <- p[,.(row1,row2,pos1,pos2,gap_bp,n_common,raw_beta_similarity,ordinal_similarity,
               dcor_raw,xi_raw_mean,xi_levels_mean,constant_raw_any)]
  card[,`:=`(null_mean_dcor_squared=z[,"q0"],signed_excess=z[,"excess"],
    dependence_credit=z[,"credit"],dependence_credit_defined=is.finite(z[,"excess"]))]
  for (ell in c(50,100,200)) card[,paste0("formula_",ell,"bp"):=
    gated_affinity(raw_beta_similarity,dependence_credit,gap_bp,ell)]
  card[,formula_100bp_uncorrected:=gated_affinity(raw_beta_similarity,
    fifelse(is.finite(dcor_raw),dcor_raw^2,0),gap_bp,100)]
  card[,old_spatial_blend_100bp:=exp(-gap_bp/100)*(raw_beta_similarity+dcor_raw)/2]

  # Reuse all training folds to audit patient influence, while leaving scores
  # and eligibility in each fold independent of the held-out beta values.
  folds <- vector("list",18)
  loo_min <- rep(Inf,nrow(p)); loo_max <- rep(-Inf,nrow(p)); loo_undefined <- integer(nrow(p))
  for (held in 1:18) {
    tx <- x[,-held,drop=FALSE]; ty <- y[,-held,drop=FALSE]
    zt <- profile_components(tx,ty); s <- zt[,"similarity"]; credit <- zt[,"credit"]
    dc <- sqrt(zt[,"q"])
    f <- gated_affinity(s,credit,p$gap_bp,100)
    loo_min <- pmin(loo_min,f); loo_max <- pmax(loo_max,f)
    loo_undefined <- loo_undefined+as.integer(!is.finite(zt[,"excess"]))
    eligible <- which(is.finite(dc)); k <- ceiling(length(eligible)/10)
    scores <- list(manhattan_raw=s,dcor_raw=dc,gap_only=-p$gap_bp,
      old_spatial_blend_100bp=exp(-p$gap_bp/100)*(s+dc)/2,
      formula_50bp=gated_affinity(s,credit,p$gap_bp,50),formula_100bp=f,
      formula_200bp=gated_affinity(s,credit,p$gap_bp,200),
      formula_100bp_uncorrected=gated_affinity(s,ifelse(is.finite(dc),dc^2,0),p$gap_bp,100))
    err <- abs(x[,held]-y[,held])
    folds[[held]] <- rbindlist(lapply(names(scores),function(method) {
      score <- scores[[method]]
      ids <- eligible[order(-score[eligible],p$pos1[eligible],p$pos2[eligible])[seq_len(k)]]
      data.table(held_patient=held,method=method,eligible_pairs=length(eligible),selected_pairs=k,
        heldout_mean_absolute_difference=mean(err[ids]),
        median_selected_gap_bp=median(p$gap_bp[ids]))
    }))
    cat("Completed held-patient fold",held,"of 18\n")
  }
  f <- rbindlist(folds); write_tab(f,"HeldPatientFolds")
  summary <- f[,.(folds=.N,mean_heldout_absolute_difference=mean(heldout_mean_absolute_difference),
    min_patient_mean_difference=min(heldout_mean_absolute_difference),
    max_patient_mean_difference=max(heldout_mean_absolute_difference),
    median_selected_gap_bp=median(median_selected_gap_bp)),by=method]
  write_tab(summary,"HeldPatientSummary")
  card[,`:=`(loo_min_formula_100bp=loo_min,loo_max_formula_100bp=loo_max,
             loo_undefined_credit_count=loo_undefined)]
  card[,max_loo_change:=pmax(abs(formula_100bp-loo_min),abs(formula_100bp-loo_max))]
  fwrite(card,file.path(out,"Item5_Formula_PairScores.csv.gz"))
  write_tab(card[,.(pairs=.N,median_score=median(formula_100bp),
    median_loo_change=median(max_loo_change),q95_loo_change=quantile(max_loo_change,.95)),
    by=.(dependence_credit_defined,loses_credit_on_deletion=loo_undefined_credit_count>0)],"SupportSummary")

  # Restore agreement and depict rare support explicitly, without inventing read counts.
  cases <- list(identical_variable=list(x=seq(.05,.95,length.out=18),y=seq(.05,.95,length.out=18)),
    shifted_variable=list(x=seq(.05,.25,length.out=18),y=seq(.70,.90,length.out=18)),
    inverse_variable=list(x=seq(.05,.95,length.out=18),y=seq(.95,.05,length.out=18)),
    identical_constant=list(x=rep(.1,18),y=rep(.1,18)),
    shared_singleton=list(x=c(.95,rep(.05,17)),y=c(.95,rep(.05,17))))
  examples <- rbindlist(lapply(names(cases),function(name) {
    xx <- cases[[name]]$x; yy <- cases[[name]]$y
    e <- dependence_credit(xx,yy); s <- 1-mean(abs(xx-yy))
    rbindlist(lapply(c(8,100,500),function(g) {
      data.table(example=name,gap_bp=g,manhattan_similarity=s,dcor=sqrt(e["q"]),
        q0=e["q0"],credit=e["credit"],formula=gated_affinity(s,e["credit"],g))
    }))
  }))
  write_tab(examples,"Examples")

  # Independent verification that unchanged baselines reproduce the previous pilot.
  prior <- fread(file.path(out,"Item5_Spatial_HeldPatientFolds.csv"))
  current <- f[method %in% c("manhattan_raw","dcor_raw","gap_only","old_spatial_blend_100bp")]
  current[method=="old_spatial_blend_100bp",method:="spatial_100bp"]
  compare <- merge(current,prior,by=c("held_patient","method"),suffixes=c("_new","_old"))
  stopifnot(nrow(compare)==72L,
    max(abs(compare$heldout_mean_absolute_difference_new-compare$heldout_mean_absolute_difference_old))<1e-10,
    all(card$formula_100bp>=0 & card$formula_100bp<=card$raw_beta_similarity+1e-12))
  capture.output(sessionInfo(),file=file.path(out,"Item5_Formula_SessionInfo.txt"))
  print(summary)
  cat("PASS: all baseline folds reproduced; scores bounded by Manhattan agreement\n")
}

if (sys.nframe()==0L) {
  arg <- grep("^--file=",commandArgs(FALSE),value=TRUE)
  script <- normalizePath(sub("^--file=","",arg[1L]))
  run_formula(dirname(dirname(script)),script)
}
