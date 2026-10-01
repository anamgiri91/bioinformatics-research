# Independent small-count checks for the shared-fragment moment engine.
# This script uses synthetic counts only and writes no data/results files.
suppressPackageStartupMessages(library(data.table))
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- dirname(dirname(normalizePath(sub("^--file=", "", arg[1L]))))
source(file.path(root, "Scripts", "SharedReadNoise.Core.R"))

checks <- 0L
eq <- function(actual, expected, label, tolerance = 1e-11) {
  if (!isTRUE(all.equal(as.numeric(actual), as.numeric(expected),
                       tolerance = tolerance, check.attributes = FALSE)))
    stop(label, ": actual=", paste(actual, collapse = ","),
         " expected=", paste(expected, collapse = ","))
  checks <<- checks + 1L
}
yes <- function(value, label) {
  if (!isTRUE(value)) stop(label)
  checks <<- checks + 1L
}
fails <- function(expr, label) {
  errored <- tryCatch({ force(expr); FALSE }, error = function(e) TRUE)
  yes(errored, label)
}
ct <- function(Ni, Mi, Nj, Mj, n00 = 0, n01 = 0, n10 = 0, n11 = 0)
  data.table(Ni, Mi, Nj, Mj, n00, n01, n10, n11)

# Hand-computed asymmetric partial overlap: using all-read betas in place
# of shared marginals would give 1/18 rather than the correct 1/12.
partial <- ct(3, 2, 4, 2, n00 = 1, n11 = 1)
z <- srn_contributions(partial)
eq(z$noise_cov, 1/12, "partial-overlap shared covariance")
eq(z$vx, 1/9, "partial-overlap site i variance")
eq(z$vy, 1/12, "partial-overlap site j variance")
negative <- srn_contributions(ct(2, 1, 2, 1, n01 = 1, n10 = 1))
eq(negative$noise_cov, -1/4, "negative covariance must remain negative")
eq(negative$d2-negative$vx-negative$vy+2*negative$noise_cov, -1,
   "signed disagreement estimate must not be clipped")
zero <- srn_contributions(ct(3, 2, 4, 2))
eq(zero$noise_cov, 0, "r=0 exact zero covariance correction")
yes(isTRUE(zero$eligible), "r=0 remains eligible")
single <- srn_contributions(ct(2, 1, 2, 1, n11 = 1))
yes(!single$eligible, "r=1 must be unsupported")
yes(nchar(single$row_status) > 0, "r=1 needs explicit row status")
low <- srn_contributions(ct(1, 0, 2, 1))
yes(!low$eligible, "site coverage below two must be ineligible")
empty <- srn_contributions(ct(0, 0, 0, 0))
yes(!empty$eligible, "zero coverage must be ineligible")

# Count corruption is not equivalent to a merely unsupported observation.
for (field in c("Ni", "Mi", "Nj", "Mj", "n00", "n01", "n10", "n11")) {
  for (bad in c(-1, 0.5, NA_real_, Inf)) {
    d <- copy(partial); set(d, j = field, value = bad)
    fails(srn_validate_counts(d), paste("reject", field, "=", bad))
  }
}
fails(srn_validate_counts(ct(2, 3, 2, 1)), "reject methylated count above total")
fails(srn_validate_counts(ct(2, 1, 2, 1, n00 = 3)), "reject overlap above depth")
fails(srn_validate_counts(ct(3, 0, 4, 2, n00 = 1, n11 = 1)),
      "reject shared methylated count above site total")
fails(srn_validate_counts(ct(3, 3, 4, 2, n00 = 1, n11 = 1)),
      "reject nonshared methylated count above nonshared depth")
big <- srn_contributions(ct(100000L, 50000L, 100000L, 50000L,
                           n00 = 50000L, n11 = 50000L))
eq(big$noise_cov, 0.25/99999, "integer products must not overflow")

# Independent exhaustive 128-state enumeration: two shared fragments,
# one i-only fragment, and two j-only fragments. The expectation target is
# obtained from p=.3, q=.6 and q11=.1, not from the engine's algebra.
states <- rbind(c(0, 0), c(0, 1), c(1, 0), c(1, 1))
probs <- c(.2, .5, .2, .1)
events <- expand.grid(a = 1:4, b = 1:4, u = 0:1, v = 0:1, w = 0:1)
enum <- rbindlist(lapply(seq_len(nrow(events)), function(k) {
  e <- events[k, ]; shared <- states[c(e$a, e$b), , drop = FALSE]
  count <- tabulate(c(e$a, e$b), nbins = 4)
  row <- ct(3, sum(shared[, 1]) + e$u,
            4, sum(shared[, 2]) + e$v + e$w,
            count[1], count[2], count[3], count[4])
  row[, probability := probs[e$a] * probs[e$b] *
        ifelse(e$u == 1, .3, .7) * ifelse(e$v == 1, .6, .4) *
        ifelse(e$w == 1, .6, .4)]
  row
}))
z <- srn_contributions(enum)
eq(sum(enum$probability), 1, "enumeration probability mass")
eq(sum(enum$probability * z$vx), 7/100, "exact expected site i noise")
eq(sum(enum$probability * z$vy), 3/50, "exact expected site j noise")
eq(sum(enum$probability * z$noise_cov), -1/75,
   "exact expected negative shared covariance")
eq(sum(enum$probability * (z$d2-z$vx-z$vy+2*z$noise_cov)), 9/100,
   "exact expected latent squared disagreement")

# Explicit fragment expansion is intentionally independent of vectorized
# PAT parsing. Dots remain missing, and pattern rows never share a fragment.
slow_pat <- function(p, offset, positions, max_gap = 200) {
  nsite <- length(positions)
  fragments <- list()
  for (k in seq_len(nrow(p))) {
    chars <- strsplit(p$pat[k], "", fixed = TRUE)[[1L]]
    at <- p$idx[k] - offset + seq_along(chars) - 1L
    stopifnot(all(at >= 1L & at <= nsite))
    x <- rep(NA_integer_, nsite)
    x[at[chars == "C"]] <- 1L
    x[at[chars == "T"]] <- 0L
    for (replicate in seq_len(p$n[k])) fragments[[length(fragments) + 1L]] <- x
  }
  mat <- do.call(rbind, fragments)
  out <- lapply(seq_len(nsite - 1L), function(i) {
    x <- mat[, i]; y <- mat[, i + 1L]
    Ni <- sum(!is.na(x)); Nj <- sum(!is.na(y))
    if (Ni == 0L || Nj == 0L || positions[i + 1L] - positions[i] > max_gap)
      return(NULL)
    d <- ct(Ni, sum(x, na.rm = TRUE), Nj, sum(y, na.rm = TRUE),
            sum(x == 0 & y == 0, na.rm = TRUE),
            sum(x == 0 & y == 1, na.rm = TRUE),
            sum(x == 1 & y == 0, na.rm = TRUE),
            sum(x == 1 & y == 1, na.rm = TRUE))
    d[, idx := offset + i]
    d
  })
  rbindlist(out)
}

# Stream one donor at a time, then compare against R's batch covariance on
# explicitly count-eligible donors. Counts are built from individual calls.
set.seed(20261001L)
draw_counts <- function(Ni, Nj, r) {
  shared <- states[sample.int(4L, r, replace = TRUE), , drop = FALSE]
  nstate <- tabulate(shared[, 1]*2L + shared[, 2]+1L, nbins = 4L)
  ct(Ni, sum(shared[, 1])+rbinom(1, Ni-r, .4),
     Nj, sum(shared[, 2])+rbinom(1, Nj-r, .6),
     nstate[1], nstate[2], nstate[3], nstate[4])
}
by_donor <- lapply(1:25, function(k) {
  d <- rbindlist(list(draw_counts(5+k%%4, 7+k%%3, if(k%%4==0) 0 else 2),
                          draw_counts(5, 8, if(k<=8) 1 else 2),
                          draw_counts(if(k<=3) 1 else 4, 6, 0)))
  d[, pair := 1:3]
  d
})
acc <- srn_accumulator(3L)
for (d in by_donor) srn_update(acc, d, d$pair)
fit <- srn_finish(acc)
eq(fit$index, c(1, 3), "pair donor thresholds must use each pair's mask")
eq(acc$n, c(25, 17, 22), "eligible counts must exclude r=1 and low depth")
all_counts <- rbindlist(by_donor)
for (i in fit$index) {
  dd <- all_counts[pair==i]
  rr <- dd$n00+dd$n01+dd$n10+dd$n11
  use <- dd$Ni>=2 & dd$Nj>=2 & rr!=1
  dd <- dd[use]; zz <- srn_contributions(dd)
  x <- dd$Mi/dd$Ni; y <- dd$Mj/dd$Nj
  row <- fit[index==i]
  eq(row$n, length(x), "matched donor count")
  eq(row$raw_cov, cov(x, y), "Welford covariance vs stats::cov")
  eq(row$raw_var_i, var(x), "Welford variance i vs stats::var")
  eq(row$raw_var_j, var(y), "Welford variance j vs stats::var")
  eq(row$cov_corrected, cov(x,y)-mean(zz$noise_cov), "mean noise uses n denominator")
  eq(row$var_i_corrected, var(x)-mean(zz$vx), "variance correction i")
  eq(row$var_j_corrected, var(y)-mean(zz$vy), "variance correction j")
  eq(row$raw_d2, mean((x-y)^2), "raw squared disagreement")
  eq(row$d2_corrected, mean((x-y)^2-zz$vx-zz$vy+2*zz$noise_cov),
     "corrected squared disagreement uses same mask")
  eq(row$manhattan_similarity, 1-mean(abs(x-y)), "Manhattan same donor mask")
}
permuted <- srn_accumulator(3L)
for (d in rev(by_donor)) srn_update(permuted, d[c(3, 1, 2)], c(3, 1, 2))
other <- srn_finish(permuted)
for (field in c("raw_cov", "var_i_corrected", "var_j_corrected", "cov_corrected", "d2_corrected"))
  eq(other[[field]], fit[[field]], paste("donor and row permutation", field))
fails(srn_update(srn_accumulator(2L), rbind(partial, partial), c(1,1)),
      "duplicate pair index within a donor must fail")

# Small-sample exact check catches replacing mean(c) by sum(c)/(n-1).
# Two donors have constant latent p_i=p_j=1/2 and two perfectly linked
# fragments each. Integrate all nine binomial count combinations.
out <- numeric(9); weights <- numeric(9); k <- 0L
for (m1 in 0:2) for (m2 in 0:2) {
  k <- k+1L; ac <- srn_accumulator(1L)
  for (m in c(m1,m2)) srn_update(ac, ct(2,m,2,m,n00=2-m,n11=m),1L)
  out[k] <- srn_finish(ac,min_donors=2L)$cov_corrected
  weights[k] <- dbinom(m1,2,.5)*dbinom(m2,2,.5)
}
eq(sum(weights*out), 0, "exact covariance expectation across two donors")

# Invalid ratios are observable outcomes, not zero or clipped coefficients.
constant <- srn_accumulator(1L)
for (k in 1:3) srn_update(constant, ct(2,1,2,1),1L)
f <- srn_finish(constant,min_donors=2)
yes(f$var_i_corrected<0 && f$var_j_corrected<0, "negative corrected variances retained")
yes(is.na(f$rho_corrected) && f$rho_status=="nonpositive_corrected_variance",
    "two negative variances must not create a valid ratio")
outside <- srn_accumulator(1L)
for (m in c(0,2,4)) srn_update(outside,ct(4,m,4,m),1L)
f <- srn_finish(outside,min_donors=2)
eq(f$rho_corrected, 9/8, "out-of-range ratio must be retained")
yes(f$rho_status=="outside_unit_interval", "out-of-range ratio must be flagged")
exact_constant <- srn_accumulator(1L)
for (k in 1:25) srn_update(exact_constant,ct(3,2,3,2),1L)
eq(srn_finish(exact_constant)$raw_var_i, 0, "exact constant must have zero raw variance")

positions <- c(100,103,110,111,311,512)
p <- data.table(idx=c(101,101,102,104,105,106,101,102),
                pat=c("C.T","TC","CTC","TT","C","T","T","C"),
                n=c(2,1,3,2,4,2,1,1))
compare_pat <- function(p, offset, positions, label) {
  expected <- slow_pat(p, offset, positions)
  actual <- srn_pat_counts(p,offset,length(positions),positions)
  eq(actual$idx,expected$idx,paste(label,"pair universe"))
  for (field in c("Ni","Mi","Nj","Mj","n00","n01","n10","n11"))
    eq(actual[[field]],expected[[field]],paste(label,field))
  invisible(actual)
}
parsed <- compare_pat(p,100L,positions,"explicit PAT expansion")
eq(parsed$idx,c(101,102,103,104),"gap 200 included and 201 excluded")
no_bridge <- data.table(idx=c(1,2),pat=c("C.T","T"),n=c(2,2))
nb <- compare_pat(no_bridge,0,c(10,20,30),"missing-call bridge")
eq(nb$n00+nb$n01+nb$n10+nb$n11,c(0,0),"dots must not bridge adjacent reference sites")
boundary <- data.table(idx=c(1,2),pat=c("C","T"),n=c(3,4))
bb <- compare_pat(boundary,0,c(10,20),"PAT row boundary")
eq(sum(bb$n00+bb$n01+bb$n10+bb$n11),0,"separate pattern rows cannot share a fragment")
fails(srn_pat_counts(data.table(idx=1,pat="CX",n=1),0,2,c(10,20)),"bad PAT alphabet")
fails(srn_pat_counts(data.table(idx=1,pat="CC",n=0),0,2,c(10,20)),"zero PAT multiplicity")
fails(srn_pat_counts(data.table(idx=2,pat="CC",n=1),0,2,c(10,20)),"PAT crossing chromosome end")

for (replicate in 1:20) {
  starts <- sample(1:8,20,replace=TRUE)
  lengths <- vapply(starts,function(i)sample.int(min(4,9-i),1),integer(1))
  pats <- vapply(lengths,function(n)paste(sample(c("C","T","."),n,replace=TRUE),collapse=""),character(1))
  pp <- data.table(idx=starts+50,pat=pats,n=sample(1:4,20,replace=TRUE))
  compare_pat(pp,50,c(10,20,30,40,240,441,450,460),paste("random fragment expansion",replicate))
}
cat("PASS:", checks, "independent count, PAT, streaming and exact moment checks\n")
