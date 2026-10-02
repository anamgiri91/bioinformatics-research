# Known-target simulations for linear decay/noise contrasts, not a validation
# of the squared-error endpoint interval or of all WGBS missingness designs.
suppressPackageStartupMessages({library(data.table);library(jsonlite)})
source("Scripts/SharedReadNoise.Decay.R")
B<-40L; P<-B*4L; REPS<-300L; DRAWS<-400L
scenarios<-data.table(name=c("flat_29","decay_57","sperm_like_52","no_shared_noise_52"),
  donors=c(29,57,52,52),mu=c(.5,.5,.85,.85),a=c(.1,.1,.03,.03),
  rho_near=c(.4,.65,.65,.65),rho_far=c(.4,.15,.15,.15),
  eta_near=c(.9,.9,.9,0),eta_far=c(.2,.2,.2,0),depth=c(4L,4L,2L,2L))
fwrite(scenarios,"Results/Task67_CalibrationScenarios.csv")
meta<-data.table(chrn=rep(rep(1:4,each=10),each=4),
  pos=rep(seq_len(B)*1000000L,each=4)+rep(c(1,10,100,200),B),
  gap=rep(c(5,5,175,175),B),depth_A=4,density=5)
design<-decay_design(meta); pr<-rep(seq_len(B),each=4);near<-meta$gap<10
truth<-function(s){
  t<-CJ(x=c(-1,1),y=c(-1,1),u=c(-1,1))
  noise<-function(rho,eta){
    pi<-s$mu+s$a*t$x+.01*t$u;pj<-s$mu+s$a*(rho*t$x+sqrt(1-rho^2)*t$y)+.01*t$u
    mean(eta*(pmin(pi,pj)-pi*pj))/s$depth
  }
  c(corrected=s$a^2*(s$rho_near-s$rho_far),
    noise_change=noise(s$rho_near,s$eta_near)-noise(s$rho_far,s$eta_far))
}
simulate<-function(s){
  D<-s$donors;U<-sample(c(-1,1),D,TRUE)
  X<-matrix(sample(c(-1,1),B*D,TRUE),B,D)[pr,]
  Y<-matrix(sample(c(-1,1),B*D,TRUE),B,D)[pr,]
  rho<-ifelse(near,s$rho_near,s$rho_far);eta<-ifelse(near,s$eta_near,s$eta_far)
  pi<-s$mu+s$a*X+matrix(rep(.01*U,each=P),P,D)
  pj<-s$mu+s$a*(rho*X+sqrt(1-rho^2)*Y)+matrix(rep(.01*U,each=P),P,D)
  q11<-pi*pj+eta*(pmin(pi,pj)-pi*pj);q10<-pi-q11;q01<-pj-q11
  stopifnot(all(pi>0 & pi<1),all(pj>0 & pj<1))
  one<-function(){
    n11<-matrix(rbinom(P*D,s$depth,q11),P,D)
    n10<-matrix(rbinom(P*D,s$depth-n11,q10/(1-q11)),P,D)
    n01<-matrix(rbinom(P*D,s$depth-n11-n10,q01/(1-q11-q10)),P,D)
    x<-(n10+n11)/s$depth;y<-(n01+n11)/s$depth
    list(x=x,y=y,c=(n11/s$depth-x*y)/(s$depth-1))
  }
  A<-one();BB<-one();CC<-one()
  decay_prepare(list(E=matrix(1,P,D),xA=A$x,yA=A$y,cA=A$c,
    biB=BB$x,bjB=BB$y,biC=CC$x,bjC=CC$y))
}
answer<-list()
for(k in seq_len(nrow(scenarios))){
  s<-scenarios[k];target<-truth(s);set.seed(67100L+k)
  rows<-lapply(seq_len(REPS),function(r){
    pre<-simulate(s);est<-decay_metrics(decay_values(pre),design)
    # Use independent, deterministic simulation/bootstrap seed streams.
    bt<-decay_bootstrap(pre,design,draws=DRAWS,seed=671000L+k*REPS+r)
    set.seed(672000L+k*REPS+r)
    ci<-apply(bt,2,quantile,c(.025,.975));ci<-sweep(ci,2,colMeans(bt)-est,"-")
    if(r %% 50L==0L)cat(s$name,r,"/",REPS,"\n")
    rbindlist(lapply(names(target),function(v)data.table(metric=v,target=target[[v]],
      estimate=est[[v]],lo=ci[1,v],hi=ci[2,v])))
  })
  z<-rbindlist(rows)
  answer[[k]]<-z[,.(target=first(target),mean_estimate=mean(estimate),
    bias=mean(estimate-target),coverage=mean(lo<=target & target<=hi),
    coverage_mcse=sqrt(mean(lo<=target & target<=hi)*(1-mean(lo<=target & target<=hi))/.N),
    mean_interval_width=mean(hi-lo),replicates=.N,bootstrap_draws=DRAWS),by=metric][,scenario:=s$name]
  fwrite(rbindlist(answer),"Results/Task67_DecayCalibration.csv")
}
cat("done\n");print(rbindlist(answer))
