# Absolute covariance-decay contrasts; no covariance or correlation ratios.
suppressPackageStartupMessages(library(data.table))

decay_design <- function(meta, min_per_band=2L, within_block=TRUE) {
  z <- copy(as.data.table(meta))
  stopifnot(all(c("chrn","pos","gap","depth_A","density") %in% names(z)))
  z[, row_id := .I]
  z[, block := paste(chrn,(pos-1L)%/%1000000L,sep=":")]
  z[, depth_bin := cut(depth_A,c(2,5,10,20,Inf),right=FALSE)]
  z[, density_bin := cut(density,c(0,10,30,Inf),right=FALSE)]
  z[, band := fifelse(gap<=10,"near",fifelse(gap>150 & gap<=200,"far",NA_character_))]
  z[, cell := if (within_block) paste(block,depth_bin,density_bin,sep="|") else
    paste(chrn,depth_bin,density_bin,sep="|")]
  counts <- z[!is.na(band) & !is.na(depth_bin) & !is.na(density_bin),.N,by=.(cell,band)]
  cells <- dcast(counts,cell~band,value.var="N",fill=0)
  if (!all(c("near","far") %in% names(cells))) stop("no near/far support")
  cells <- cells[near>=min_per_band & far>=min_per_band]
  if (!nrow(cells)) stop("no common support")
  cells[, mass := pmin(near,far)]
  cells[, target := mass/sum(mass)]
  z <- merge(z[!is.na(band)],cells,by="cell",sort=FALSE)
  z[, weight := target/fifelse(band=="near",near,far)]
  stopifnot(abs(sum(z[band=="near",weight])-1)<1e-12,
            abs(sum(z[band=="far",weight])-1)<1e-12)
  z[order(row_id)]
}

decay_metrics <- function(values, design, block_weight=rep(1,nrow(design))) {
  stopifnot(nrow(values)==nrow(design),ncol(values)==3L,
            length(block_weight)==nrow(design),all(block_weight>=0))
  w <- design$weight*block_weight
  near <- design$band=="near"; far <- !near
  if (any(!is.finite(values)) || sum(w[near])==0 || sum(w[far])==0)
    return(setNames(rep(NA_real_,7),c("raw","corrected","reference","noise_change",
                  "raw_reference_error","corrected_reference_error","delta_squared_error")))
  means <- function(k) colSums(values[k,,drop=FALSE]*w[k])/sum(w[k])
  contrast <- means(near)-means(far)
  W <- contrast[1]; U <- contrast[2]; R <- contrast[3]
  c(raw=unname(W),corrected=unname(U),reference=unname(R),noise_change=unname(W-U),
    raw_reference_error=unname(W-R),corrected_reference_error=unname(U-R),
    delta_squared_error=unname((U-R)^2-(W-R)^2))
}

decay_prepare <- function(m) {
  E <- m$E
  list(E=E, X=E*m$xA, Y=E*m$yA, XY=E*m$xA*m$yA, C=E*m$cA,
       BI=E*m$biB, CJ=E*m$bjC, BICJ=E*m$biB*m$bjC,
       CI=E*m$biC, BJ=E*m$bjB, CIBJ=E*m$biC*m$bjB)
}

decay_values <- function(pre, donor_weight=rep(1,ncol(pre$E))) {
  sums <- lapply(pre,function(x) drop(x %*% donor_weight)); n <- sums$E
  cv <- function(xy,x,y) (xy-x*y/n)/(n-1)
  W <- cv(sums$XY,sums$X,sums$Y)
  U <- W-sums$C/n
  R <- (cv(sums$BICJ,sums$BI,sums$CJ)+cv(sums$CIBJ,sums$CI,sums$BJ))/2
  out <- cbind(raw=W,corrected=U,reference=R); out[n<2,] <- NA_real_
  out
}

decay_bootstrap <- function(pre,design,draws=2000L,seed=2026100201L) {
  blocks <- unique(design[,.(chrn,block)])
  groups <- split(seq_len(nrow(blocks)),blocks$chrn)
  bi <- match(design$block,blocks$block); D <- ncol(pre$E)
  set.seed(seed)
  ans <- vapply(seq_len(draws),function(b) {
    dw <- tabulate(sample.int(D,D,replace=TRUE),D)
    bw <- numeric(nrow(blocks))
    for (g in groups) bw[g] <- tabulate(sample.int(length(g),length(g),replace=TRUE),length(g))
    out <- decay_metrics(decay_values(pre,dw),design,bw[bi])
    if (b %% 250L==0L) cat("bootstrap",b,"/",draws,"\n")
    out
  },numeric(7))
  t(ans)
}
