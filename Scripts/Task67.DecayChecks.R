source("Scripts/SharedReadNoise.Decay.R")
set.seed(671)
m <- lapply(seq_len(7),function(j) matrix(runif(48),8,6))
names(m)<-c("xA","yA","cA","biB","bjB","biC","bjC")
m$cA <- m$cA/50; m$E <- matrix(1,8,6);m$E[cbind(c(2,4),c(1,6))]<-0
pre <- decay_prepare(m);w <- c(1,2,0,1,1,1)
v <- decay_values(pre,w)
brute <- t(vapply(seq_len(8),function(i){
  ix <- rep(seq_len(6),w);ix<-ix[m$E[i,ix]==1]
  W<-cov(m$xA[i,ix],m$yA[i,ix]);U<-W-mean(m$cA[i,ix])
  R<-(cov(m$biB[i,ix],m$bjC[i,ix])+cov(m$biC[i,ix],m$bjB[i,ix]))/2
  c(W,U,R)
},numeric(3)))
stopifnot(max(abs(v-brute))<1e-12)
meta<-data.table(chrn=rep(1:2,each=4),pos=rep(c(1,5,50,100),2),
  gap=rep(c(5,5,180,180),2),depth_A=rep(c(7,8,8,9),2),density=5)
d<-decay_design(meta);stopifnot(nrow(d)==8L,uniqueN(d$cell)==2L)
balance<-d[,.(n=sum(weight[band=="near"]),f=sum(weight[band=="far"])),by=cell]
stopifnot(max(abs(balance$n-balance$f))<1e-12)
vals<-cbind(raw=ifelse(d$band=="near",2,1),corrected=1,reference=1)
z<-decay_metrics(vals,d);stopifnot(z["raw"]==1,z["noise_change"]==1,z["delta_squared_error"] == -1)
ord<-sample(nrow(d));stopifnot(identical(z,decay_metrics(vals[ord,],d[ord])))
vals[,2]<-vals[,1];stopifnot(decay_metrics(vals,d)["delta_squared_error"]==0)
vals[1,1]<-NA;stopifnot(all(is.na(decay_metrics(vals,d))))
stopifnot(inherits(try(decay_design(meta[gap<10]),silent=TRUE),"try-error"))
cat("PASS: independent donor-resampling covariance, matched-cell balance, known contrast/loss, order invariance, zero correction and failure handling\n")
