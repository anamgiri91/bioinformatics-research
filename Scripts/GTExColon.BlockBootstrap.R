# Exact multiplicity implementation of whole-block resampling.
# Avoids copying and re-sorting millions of rows on every bootstrap draw.
prepare_block_stat <- function(d, kind) {
  stopifnot(nrow(d)>0L,all(is.finite(d$linkage)),!anyNA(d$block))
  block_levels <- sort(unique(d$block))
  bi <- match(d$block,block_levels); nb <- length(block_levels)
  if (kind=="mean") {
    count <- tabulate(bi,nb)
    total <- as.numeric(rowsum(d$linkage,bi,reorder=TRUE))
    return(list(nb=nb,stat=function(mult) sum(mult*total)/sum(mult*count)))
  }
  rank_layout <- function(v) {
    o <- order(v,method="radix"); v <- v[o]
    ends <- c(which(diff(v)!=0),length(v))
    group <- integer(length(v)); group[o] <- rep(seq_along(ends),diff(c(0,ends)))
    list(order=o,ends=ends,group=group)
  }
  weighted_ranks <- function(layout,w) {
    cumulative <- cumsum(w[layout$order])[layout$ends]
    count <- diff(c(0,cumulative))
    (cumulative-(count-1)/2)[layout$group]
  }
  if (kind=="spearman") {
    stopifnot(all(is.finite(d$gap)))
    lx <- rank_layout(d$linkage); ly <- rank_layout(d$gap)
    stat <- function(mult) {
      w <- mult[bi]; total <- sum(w)
      # Average ranks with ties recomputed under the replicated block counts.
      rx <- weighted_ranks(lx,w)-(total+1)/2
      ry <- weighted_ranks(ly,w)-(total+1)/2
      den <- sqrt(sum(w*rx^2)*sum(w*ry^2))
      if (den==0) return(NA_real_)
      sum(w*rx*ry)/den
    }
    return(list(nb=nb,stat=stat))
  }
  if (kind=="gain_difference") {
    stopifnot(all(is.finite(d$heldout_gain)))
    o <- order(d$linkage,method="radix"); x <- d$linkage[o]; gain <- d$heldout_gain[o]
    bix <- bi[o]
    stat <- function(mult) {
      w <- mult[bix]; cw <- cumsum(w); total <- tail(cw,1)
      # Exactly reproduces R median, including even-sized resamples and ties.
      at <- findInterval(c(ceiling(total/2),floor(total/2)+1)-1,cw)+1L
      md <- mean(x[at]); cut <- findInterval(md,x)
      low_n <- if(cut) cw[cut] else 0
      if(low_n==0 || low_n==total) return(NA_real_)
      cg <- cumsum(w*gain); low_gain <- if(cut) cg[cut] else 0
      (tail(cg,1)-low_gain)/(total-low_n)-low_gain/low_n
    }
    return(list(nb=nb,stat=stat))
  }
  stop("Unknown block statistic")
}

block_boot_exact <- function(d,kind,B=1000L,seed=20260929L,progress=TRUE) {
  prepared <- prepare_block_stat(d,kind)
  set.seed(seed)
  values <- vapply(seq_len(B),function(b) {
    mult <- tabulate(sample.int(prepared$nb,prepared$nb,replace=TRUE),prepared$nb)
    value <- prepared$stat(mult)
    if(progress && b%%100L==0L) cat("  bootstrap",kind,b,"/",B,"on",nrow(d),"pairs\n")
    value
  },0.0)
  if(all(!is.finite(values))) return(c(NA_real_,NA_real_))
  unname(quantile(values,c(.025,.975),na.rm=TRUE))
}

check_block_bootstrap <- function() {
  set.seed(52001)
  d <- data.frame(block=rep(c("1 0","1 1","2 0","2 1","3 1"),c(4,7,5,8,6)),
    linkage=sample(c(-.1,0,.1,.2,.3),30,TRUE),gap=sample(2:6,30,TRUE),heldout_gain=rnorm(30))
  blocks <- split(seq_len(nrow(d)),d$block)
  for(kind in c("mean","spearman","gain_difference")) {
    p <- prepare_block_stat(d,kind)
    for(b in 1:50) {
      draw <- sample.int(p$nb,p$nb,replace=TRUE)
      z <- d[unlist(blocks[draw],use.names=FALSE),]
      expected <- switch(kind,mean=mean(z$linkage),spearman=cor(z$linkage,z$gap,method="spearman"),
        gain_difference={md<-median(z$linkage);mean(z$heldout_gain[z$linkage>md])-mean(z$heldout_gain[z$linkage<=md])})
      actual <- p$stat(tabulate(draw,p$nb))
      stopifnot(isTRUE(all.equal(actual,expected,tolerance=1e-12)))
    }
  }
  cat("PASS: 150 exact comparisons against expanded whole-block bootstrap, including ties and medians\n")
}

if(sys.nframe()==0L) check_block_bootstrap()
