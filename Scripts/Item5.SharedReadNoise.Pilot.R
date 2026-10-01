# Exploratory feasibility study; not part of the frozen Task 47 validation.
# All data here are simulated. No parameters are fitted to GTEx methylation.
suppressPackageStartupMessages({library(data.table); library(jsonlite)})
arg <- grep("^--file=",commandArgs(FALSE),value=TRUE)
root <- dirname(dirname(normalizePath(sub("^--file=","",arg[1L]))))
res <- file.path(root,"Results")
SEED <- 20260930L; REPS <- 20000L; PEOPLE <- 29L

# Check the read-noise moments by exhaustive enumeration of three fragments.
prob <- c(.4,.1,.1,.4) # 11, 10, 01, 00
z <- as.data.table(expand.grid(a=0:3,b=0:3,c=0:3))
z <- z[a+b+c<=3][,d:=3-a-b-c]
pr <- apply(as.matrix(z),1,function(v) dmultinom(v,prob=prob))
x <- (z$a+z$b)/3; y <- (z$a+z$c)/3
delta <- z$a/3-x*y
v <- x*(1-x)/2; cv <- delta/2
stopifnot(abs(sum(pr)-1)<1e-14,
          abs(sum(pr*v)-.25/3)<1e-14,
          abs(sum(pr*cv)-.15/3)<1e-14,
          abs(sum(pr*((x-y)^2-v-y*(1-y)/2+2*cv)))<1e-14)
cat("PASS: exact three-fragment enumeration for variance, covariance and squared disagreement\n")

cases <- data.table(
  scenario=c("constant_linked","independent_linked","independent_partial_overlap",
             "independent_no_overlap","independent_heterogeneous_depth",
             "linear_linked","linear_heterogeneous_depth","inverse_linked"),
  latent=c("constant","independent","independent","independent","independent","linear","linear","inverse"),
  depth=c("10","10","20_30","10","heterogeneous","10","heterogeneous","10"),
  overlap=c(1,1,.5,0,.5,1,.5,1),coupling=.8)
write_json(list(status="exploratory simulation; separate from frozen validation",seed=SEED,
  replicates=REPS,patients=PEOPLE,scenarios=cases,
  target="conditional latent sample covariance and mean squared difference",
  rule="retain signed corrected moments; do not clip or conceal invalid corrected correlations"),
  file.path(res,"Item5_SharedReadNoise_Protocol.json"),pretty=TRUE,auto_unbox=TRUE)
set.seed(SEED)
one <- function(k) {
  cfg <- cases[k]; len <- REPS*PEOPLE
  p <- if(cfg$latent=="constant") rep(.5,len) else runif(len,.3,.7)
  q <- switch(cfg$latent,constant=rep(.5,len),independent=runif(len,.3,.7),
              linear=.1+.8*p,inverse=1-p)
  if(cfg$depth=="heterogeneous") {
    ni <- sample(10:80,len,TRUE); nj <- sample(10:80,len,TRUE)
  } else {ni<-rep(if(cfg$depth=="20_30")20L else 10L,len)
          nj<-rep(if(cfg$depth=="20_30")30L else 10L,len)}
  shared <- floor(pmin(ni,nj)*cfg$overlap)
  stopifnot(all(shared==0 | shared>=2))
  p11 <- (1-cfg$coupling)*p*q+cfg$coupling*pmin(p,q)
  p10 <- p-p11; p01 <- q-p11
  a <- rbinom(len,shared,p11)
  b <- rbinom(len,shared-a,p10/(1-p11))
  c <- rbinom(len,shared-a-b,p01/(1-p11-p10))
  mx <- a+b+rbinom(len,ni-shared,p)
  my <- a+c+rbinom(len,nj-shared,q)
  x <- mx/ni; y <- my/nj
  vx <- x*(1-x)/(ni-1); vy <- y*(1-y)/(nj-1)
  noise_cov <- numeric(len); ok <- shared>=2
  delta <- a[ok]/shared[ok]-(a[ok]+b[ok])*(a[ok]+c[ok])/shared[ok]^2
  noise_cov[ok] <- shared[ok]^2*delta/(ni[ok]*nj[ok]*(shared[ok]-1))
  mat <- function(v) matrix(v,nrow=REPS,ncol=PEOPLE)
  sample_cov <- function(x,y) rowSums((x-rowMeans(x))*(y-rowMeans(y)))/(PEOPLE-1)
  xm <- mat(x); ym <- mat(y); pm <- mat(p); qm <- mat(q)
  observed_cov <- sample_cov(xm,ym); truth_cov <- sample_cov(pm,qm)
  corrected_cov <- observed_cov-rowMeans(mat(noise_cov))
  observed_vx <- sample_cov(xm,xm); observed_vy <- sample_cov(ym,ym)
  corrected_vx <- observed_vx-rowMeans(mat(vx)); corrected_vy <- observed_vy-rowMeans(mat(vy))
  observed_d2 <- rowMeans((xm-ym)^2); truth_d2 <- rowMeans((pm-qm)^2)
  corrected_d2 <- rowMeans(mat((x-y)^2-vx-vy+2*noise_cov))
  raw_r <- observed_cov/sqrt(observed_vx*observed_vy)
  positive_variance <- corrected_vx>0 & corrected_vy>0
  corr_r <- rep(NA_real_,REPS)
  corr_r[positive_variance] <- corrected_cov[positive_variance]/
    sqrt(corrected_vx[positive_variance]*corrected_vy[positive_variance])
  out <- data.table(scenario=cfg$scenario,replicate=seq_len(REPS),
    truth_cov,observed_cov,corrected_cov,truth_d2,observed_d2,corrected_d2,raw_r,
    corrected_nonpositive_variance=!positive_variance,
    corrected_outside_unit=positive_variance & abs(corr_r)>1)
  summary <- out[, .(replicates=.N,mean_true_cov=mean(truth_cov),
    mean_observed_cov=mean(observed_cov),mean_corrected_cov=mean(corrected_cov),
    raw_cov_bias=mean(observed_cov-truth_cov),corrected_cov_bias=mean(corrected_cov-truth_cov),
    corrected_cov_bias_mcse=sd(corrected_cov-truth_cov)/sqrt(.N),
    raw_cov_rmse=sqrt(mean((observed_cov-truth_cov)^2)),
    corrected_cov_rmse=sqrt(mean((corrected_cov-truth_cov)^2)),
    raw_d2_bias=mean(observed_d2-truth_d2),corrected_d2_bias=mean(corrected_d2-truth_d2),
    corrected_d2_bias_mcse=sd(corrected_d2-truth_d2)/sqrt(.N),
    median_observed_r=median(raw_r,na.rm=TRUE),
    share_corrected_nonpositive_variance=mean(corrected_nonpositive_variance),
    share_corrected_outside_unit=mean(corrected_outside_unit),
    share_negative_corrected_d2=mean(corrected_d2<0)),by=scenario]
  # A seeded Monte Carlo diagnostic supplements the exact algebra check.
  stopifnot(abs(summary$corrected_cov_bias)<5*summary$corrected_cov_bias_mcse,
            abs(summary$corrected_d2_bias)<5*summary$corrected_d2_bias_mcse)
  cat("Completed",cfg$scenario,"\n")
  list(summary=summary,replicates=out)
}
ans <- lapply(seq_len(nrow(cases)),one)
summary <- rbindlist(lapply(ans,`[[`,"summary"))
fwrite(summary,file.path(res,"Item5_SharedReadNoise_Summary.csv"))
fwrite(rbindlist(lapply(ans,`[[`,"replicates")),file.path(res,"Item5_SharedReadNoise_Replicates.csv.gz"))
print(summary[,.(scenario,raw_cov_bias,corrected_cov_bias,median_observed_r,
                share_corrected_nonpositive_variance,share_corrected_outside_unit)])
writeLines(capture.output(sessionInfo()),file.path(res,"Item5_SharedReadNoise_SessionInfo.txt"))
cat("done\n")
