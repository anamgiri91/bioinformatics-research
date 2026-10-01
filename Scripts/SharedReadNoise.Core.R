# Shared-fragment moment estimation. See plan_shared_read_noise.md for targets and assumptions.
# This module does not establish physical-fragment independence or novel methodology.
suppressPackageStartupMessages(library(data.table))

srn_validate_counts <- function(d) {
  fields <- c("Ni", "Mi", "Nj", "Mj", "n00", "n01", "n10", "n11")
  if (!all(fields %in% names(d))) stop("missing count fields")
  for (v in fields) {
    x <- d[[v]]
    if (!is.numeric(x) || any(!is.finite(x) | x < 0 | x != floor(x)))
      stop("invalid nonnegative integer count: ", v)
    if (any(x > 2^26)) stop("count exceeds exact-product arithmetic limit: ", v)
  }
  Ni <- as.double(d$Ni); Nj <- as.double(d$Nj)
  mi <- as.double(d$n10) + d$n11; mj <- as.double(d$n01) + d$n11
  r <- as.double(d$n00) + d$n01 + d$n10 + d$n11
  if (any(r > pmin(Ni, Nj) | d$Mi > Ni | d$Mj > Nj |
          d$Mi-mi < 0 | d$Mi-mi > Ni-r | d$Mj-mj < 0 | d$Mj-mj > Nj-r))
    stop("inconsistent site/shared/nonshared counts")
  invisible(TRUE)
}

srn_contributions <- function(d) {
  srn_validate_counts(d)
  z <- copy(as.data.table(d))
  for (v in c("Ni", "Mi", "Nj", "Mj", "n00", "n01", "n10", "n11"))
    set(z, j=v, value=as.double(z[[v]]))
  z[, r := n00+n01+n10+n11]
  z[, row_status := fifelse(Ni < 2 | Nj < 2, "low_coverage",
                           fifelse(r == 1, "single_shared_fragment", "eligible"))]
  z[, eligible := row_status == "eligible"]
  for (v in c("x", "y", "vx", "vy", "noise_cov", "d2", "l1"))
    set(z, j=v, value=rep(NA_real_, nrow(z)))
  z[eligible == TRUE, `:=`(x=Mi/Ni, y=Mj/Nj)]
  z[eligible == TRUE, `:=`(vx=x*(1-x)/(Ni-1), vy=y*(1-y)/(Nj-1),
                   noise_cov=0, d2=(x-y)^2, l1=abs(x-y))]
  z[eligible & r >= 2, noise_cov :=
      (r*n11-(n10+n11)*(n01+n11))/(Ni*Nj*(r-1))]
  z
}

srn_accumulator <- function(size) {
  stopifnot(length(size)==1L, is.finite(size), size >= 0, size==floor(size))
  z <- data.table(n=integer(size), observed=integer(size),
                  low_coverage=integer(size), single_shared=integer(size), zero_shared=integer(size))
  for (v in c("mean_x", "mean_y", "m2x", "m2y", "cxy", "sum_vx", "sum_vy",
              "sum_cov", "sum_d2", "sum_l1", "sum_Ni", "sum_Nj", "sum_r"))
    set(z, j=v, value=numeric(size))
  z
}

# Each call adds one donor. The caller must enforce unique donor IDs across calls.
# indices map this donor's rows to reference-adjacent pairs in one chromosome.
srn_update <- function(acc, d, indices) {
  if (length(indices)!=nrow(d) || anyNA(indices) || anyDuplicated(indices) ||
      any(indices < 1 | indices > nrow(acc) | indices != floor(indices)))
    stop("invalid or duplicate accumulator indices")
  indices <- as.integer(indices)
  z <- srn_contributions(d)
  add <- function(v, ix, values) set(acc, ix, v, acc[[v]][ix]+values)
  add("observed", indices, 1L)
  add("low_coverage", indices, as.integer(z$row_status=="low_coverage"))
  add("single_shared", indices, as.integer(z$row_status=="single_shared_fragment"))
  ok <- which(z$eligible)
  if (!length(ok)) return(invisible(acc))
  ix <- indices[ok]; z <- z[ok]
  nn <- acc$n[ix]+1L
  dx <- z$x-acc$mean_x[ix]; dy <- z$y-acc$mean_y[ix]
  mx <- acc$mean_x[ix]+dx/nn; my <- acc$mean_y[ix]+dy/nn
  # Welford updates preserve exact constants and avoid subtracting large squares.
  add("m2x", ix, dx*(z$x-mx)); add("m2y", ix, dy*(z$y-my))
  add("cxy", ix, dx*(z$y-my))
  set(acc, ix, "n", nn); set(acc, ix, "mean_x", mx); set(acc, ix, "mean_y", my)
  fields <- c(vx="sum_vx", vy="sum_vy", noise_cov="sum_cov", d2="sum_d2",
              l1="sum_l1", Ni="sum_Ni", Nj="sum_Nj", r="sum_r")
  for (source in names(fields)) add(fields[[source]], ix, z[[source]])
  add("zero_shared", ix, as.integer(z$r==0))
  invisible(acc)
}

srn_finish <- function(acc, min_donors=20L) {
  stopifnot(length(min_donors)==1L, min_donors >= 2, min_donors==floor(min_donors))
  ix <- which(acc$n >= min_donors)
  a <- acc[ix]; n <- a$n
  z <- data.table(index=ix, n=n, observed=a$observed,
    low_coverage=a$low_coverage, single_shared=a$single_shared, zero_shared=a$zero_shared,
    mean_i=a$mean_x, mean_j=a$mean_y,
    raw_cov=a$cxy/(n-1), raw_var_i=a$m2x/(n-1), raw_var_j=a$m2y/(n-1),
    mean_noise_cov=a$sum_cov/n, mean_noise_var_i=a$sum_vx/n, mean_noise_var_j=a$sum_vy/n,
    raw_d2=a$sum_d2/n, manhattan_similarity=1-a$sum_l1/n,
    mean_Ni=a$sum_Ni/n, mean_Nj=a$sum_Nj/n, mean_shared=a$sum_r/n)
  z[, `:=`(cov_corrected=raw_cov-mean_noise_cov,
            var_i_corrected=raw_var_i-mean_noise_var_i,
            var_j_corrected=raw_var_j-mean_noise_var_j,
            d2_corrected=raw_d2-mean_noise_var_i-mean_noise_var_j+2*mean_noise_cov)]
  z[, agreement_corrected := 1-d2_corrected]
  z[, `:=`(raw_pearson=NA_real_, rho_variance_only=NA_real_, rho_corrected=NA_real_,
            rho_status="nonpositive_corrected_variance")]
  z[raw_var_i > 0 & raw_var_j > 0, raw_pearson := raw_cov/sqrt(raw_var_i*raw_var_j)]
  z[var_i_corrected > 0 & var_j_corrected > 0,
    `:=`(rho_variance_only=raw_cov/sqrt(var_i_corrected*var_j_corrected),
          rho_corrected=cov_corrected/sqrt(var_i_corrected*var_j_corrected), rho_status="valid")]
  z[rho_status=="valid" & abs(rho_corrected)>1, rho_status := "outside_unit_interval"]
  z[, variance_only_status := fifelse(is.na(rho_variance_only), "nonpositive_corrected_variance",
                              fifelse(abs(rho_variance_only)>1, "outside_unit_interval", "valid"))]
  z
}

# Arithmetic reconstruction from PAT records; multiplicity is a counted unit.
# Whether that unit is an independent physical fragment is a separate provenance gate.
srn_pat_counts <- function(p, offset, ncpg, positions, max_gap=200L) {
  required <- c("idx", "pat", "n")
  if (!all(required %in% names(p)) || length(positions)!=ncpg) stop("invalid PAT/map schema")
  if (length(max_gap)!=1L || !is.finite(max_gap) || max_gap < 0 ||
      any(!is.finite(positions) | positions < 1 | positions!=floor(positions)) ||
      any(diff(positions) <= 0)) stop("invalid genomic positions or maximum gap")
  p <- as.data.table(p)
  if (!nrow(p)) return(data.table(idx=integer(),Ni=numeric(),Mi=numeric(),Nj=numeric(),Mj=numeric(),
                                  n00=numeric(),n01=numeric(),n10=numeric(),n11=numeric()))
  if (anyNA(p[, ..required]) || any(p$n <= 0 | p$n!=floor(p$n) | !is.finite(p$n)) ||
      any(p$idx!=floor(p$idx))) stop("invalid PAT record")
  L <- nchar(p$pat, type="bytes")
  if (any(L < 1 | p$idx <= offset | p$idx+L-1 > offset+ncpg)) stop("PAT outside chromosome CpG map")
  codes <- utf8ToInt(paste0(p$pat, collapse=""))
  if (any(!codes %in% c(46L,67L,84L))) stop("invalid PAT alphabet")
  row <- rep.int(seq_len(nrow(p)), L); off <- sequence(L)-1L
  called <- codes!=46L; methylated <- codes==67L
  k <- which(called)
  sites <- data.table(idx=p$idx[row[k]]+off[k], n=as.double(p$n[row[k]]), m=methylated[k])[
    , .(N=sum(n), M=sum(n[m])), keyby=idx]
  nxt <- match(sites$idx+1L, sites$idx)
  keep <- which(!is.na(nxt) & sites$idx-offset < ncpg)
  rank <- sites$idx[keep]-offset
  keep <- keep[positions[rank+1L]-positions[rank] <= max_gap]
  pairs <- data.table(idx=sites$idx[keep], Ni=sites$N[keep], Mi=sites$M[keep],
                       Nj=sites$N[nxt[keep]], Mj=sites$M[nxt[keep]])
  for (v in c("n00","n01","n10","n11")) set(pairs, j=v, value=numeric(nrow(pairs)))
  K <- length(codes)
  k <- if (K >= 2L) which(c(row[-1L]==row[-K],FALSE) & called & c(called[-1L],FALSE)) else integer()
  if (length(k) && nrow(pairs)) {
    joint <- data.table(idx=p$idx[row[k]]+off[k], n=as.double(p$n[row[k]]),
                        state=as.integer(methylated[k])*2L+as.integer(methylated[k+1L]))[
      , .(n00=sum(n[state==0L]), n01=sum(n[state==1L]),
          n10=sum(n[state==2L]), n11=sum(n[state==3L])), keyby=idx]
    at <- match(pairs$idx, joint$idx); ok <- which(!is.na(at))
    for (v in c("n00","n01","n10","n11")) set(pairs, ok, v, joint[[v]][at[ok]])
  }
  srn_validate_counts(pairs)
  attr(pairs,"site_counts") <- sites
  pairs
}
