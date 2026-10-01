# ================================================================
# Task 57 - Does the joint donor-and-block bootstrap give honest
#           intervals for the three-way benchmark endpoint?
# Date: Oct 1, 2026
#
# Plan section 10 says the joint bootstrap "must be checked in
# development simulations with donor and spatial dependence" before
# the external test. This is that check.
#
# SIMULATED GENOMES
#   2,000 pairs in 100 blocks of 20 (5 chromosomes x 20 blocks). Each
#   block has its own latent correlation (0 to 0.9), spread, mean and
#   within-fragment coupling (0 to 1), shared by its pairs: spatial
#   dependence. Each donor has a global methylation shift on all pairs:
#   donor dependence. Depth per site 6 + NegBin(mean 15, size 3);
#   85% of the smaller depth is shared. Fragments are split into A, B
#   and C exactly as in the benchmark (multinomial 1/3 by state).
#   Each replicate draws new donors AND new block parameters, so the
#   target is the average endpoint over both.
#
# SCENARIOS (fixed here)
#   coupled_29  coupling 0 to 1 by block, 29 donors (GTEx colon size)
#   coupled_57  the same with 57 donors (the external cohort size)
#   null_57     no coupling: shared reads add no covariance noise, so
#               the correction cannot help; the endpoint is >= 0
#   300 replicates each, 400 bootstrap draws per replicate.
#
# CHECKS
#   Coverage of the 95% interval for the average endpoint (Monte
#   Carlo SE about 1.3%), and in the null scenario how often the upper
#   end falls below 0 (a false claim that the correction helps).
#
# SECOND VERSION (after the first run showed the plan's percentile
#   interval covering only 38% to 82%): the same draws are also turned
#   into a bias-corrected percentile interval (shift by the mean of the
#   draws minus the estimate) and a normal interval (estimate +- 1.96
#   bootstrap SE), and a donor-only bootstrap is added. The squared-error
#   endpoint gains variance when donors are resampled, which shifts the
#   plain percentile interval upward.
#
# OUTPUT: Results/Task57_BootstrapCalibration.csv
# ================================================================

suppressPackageStartupMessages({ library(data.table) })
pick_dir <- function(...) { for (d in c(...)) if (dir.exists(d)) return(d)
  stop("no candidate directory exists") }
repo_dir <- pick_dir("/mmfs1/home/wln26/Experiments.Outlier.July31.2026",
                     path.expand("~/Desktop/bioinformatics-research"))
res <- file.path(repo_dir, "Results")
NP <- 2000L; NBLK <- 100L; NCHR <- 5L; REPS <- 300L; NB <- 400L; MIN_DON <- 20L

split3 <- function(x) {                      # multinomial split of counts into A, B, C
  a <- rbinom(length(x), x, 1/3); b <- rbinom(length(x), x - a, 1/2); list(a, b, x - a - b)
}
one_genome <- function(D, coupled) {
  blk <- rep(seq_len(NBLK), each = NP / NBLK); chr <- (blk - 1L) %/% (NBLK / NCHR) + 1L
  rho <- runif(NBLK, 0, 0.9)[blk]; sdl <- runif(NBLK, 0.3, 1.2)[blk]; mu <- rnorm(NBLK, 0, 1.5)[blk]
  eta <- if (coupled) runif(NBLK, 0, 1)[blk] else rep(0, NP)
  u <- rnorm(D, 0, 0.3)                          # donor-wide shift
  m <- NP * D; pr <- rep(seq_len(NP), times = D); dn <- rep(seq_len(D), each = NP)
  z1 <- rnorm(m); z2 <- rho[pr] * z1 + sqrt(1 - rho[pr]^2) * rnorm(m)
  pi <- plogis(mu[pr] + u[dn] + sdl[pr] * z1); pj <- plogis(mu[pr] + u[dn] + sdl[pr] * z2)
  Ni <- 6L + rnbinom(m, mu = 15, size = 3); Nj <- 6L + rnbinom(m, mu = 15, size = 3)
  r <- as.integer(round(0.85 * pmin(Ni, Nj)))
  q11 <- (1 - eta[pr]) * pi * pj + eta[pr] * pmin(pi, pj); q10 <- pmax(pi - q11, 0); q01 <- pmax(pj - q11, 0)
  n11 <- rbinom(m, r, q11); rem <- r - n11
  n10 <- rbinom(m, rem, ifelse(1 - q11 > 0, q10 / (1 - q11), 0)); rem <- rem - n10
  n01 <- rbinom(m, rem, ifelse(1 - q11 - q10 > 0, pmin(1, q01 / (1 - q11 - q10)), 0)); n00 <- rem - n01
  Mio <- rbinom(m, Ni - r, pi); Uio <- Ni - r - Mio; Mjo <- rbinom(m, Nj - r, pj); Ujo <- Nj - r - Mjo
  S <- lapply(list(n00 = n00, n01 = n01, n10 = n10, n11 = n11, Mio = Mio, Uio = Uio, Mjo = Mjo, Ujo = Ujo), split3)
  part <- function(k) {                          # counts of partition k (1 = A, 2 = B, 3 = C)
    g <- function(v) S[[v]][[k]]
    Mi <- g("n10") + g("n11") + g("Mio"); Ni_ <- g("n00") + g("n01") + g("n10") + g("n11") + g("Mio") + g("Uio")
    Mj <- g("n01") + g("n11") + g("Mjo"); Nj_ <- g("n00") + g("n01") + g("n10") + g("n11") + g("Mjo") + g("Ujo")
    list(Ni = Ni_, Mi = Mi, Nj = Nj_, Mj = Mj, n00 = g("n00"), n01 = g("n01"), n10 = g("n10"), n11 = g("n11"))
  }
  A <- part(1); B <- part(2); C <- part(3); rA <- A$n00 + A$n01 + A$n10 + A$n11
  E <- A$Ni >= 2 & A$Nj >= 2 & B$Ni >= 2 & B$Nj >= 2 & C$Ni >= 2 & C$Nj >= 2 & (rA == 0 | rA >= 2)
  mat <- function(v) matrix(ifelse(E, v, 0), NP, D)
  xA <- A$Mi / pmax(A$Ni, 1); yA <- A$Mj / pmax(A$Nj, 1)
  cA <- ifelse(rA >= 2, (rA * A$n11 - (A$n10 + A$n11) * (A$n01 + A$n11)) / (pmax(A$Ni, 1) * pmax(A$Nj, 1) * pmax(rA - 1, 1)), 0)
  biB <- B$Mi / pmax(B$Ni, 1); bjB <- B$Mj / pmax(B$Nj, 1); biC <- C$Mi / pmax(C$Ni, 1); bjC <- C$Mj / pmax(C$Nj, 1)
  list(E = matrix(as.numeric(E), NP, D), EX = mat(xA), EY = mat(yA), EXY = mat(xA * yA), EC = mat(cA),
       EBI = mat(biB), ECJ = mat(bjC), EBIC = mat(biB * bjC), ECI = mat(biC), EBJ = mat(bjB), ECIB = mat(biC * bjB),
       chr = chr, blk = blk)
}
dC_w <- function(g, w) {
  n <- drop(g$E %*% w); f <- function(M) drop(M %*% w)
  cv <- function(sxy, sx, sy) (sxy - sx * sy / n) / (n - 1)
  W <- cv(f(g$EXY), f(g$EX), f(g$EY)); U <- W - f(g$EC) / n
  R <- (cv(f(g$EBIC), f(g$EBI), f(g$ECJ)) + cv(f(g$ECIB), f(g$ECI), f(g$EBJ))) / 2
  d <- (U - R)^2 - (W - R)^2; d[n < 2] <- NA_real_; d
}
run_scenario <- function(name, D, coupled, seed) {
  set.seed(seed)
  out <- t(vapply(seq_len(REPS), function(r) {
    g <- one_genome(D, coupled)
    base_ok <- rowSums(g$E) >= MIN_DON                 # fixed pair mask, as in the benchmark
    est <- mean(dC_w(g, rep(1, D))[base_ok], na.rm = TRUE)
    blocks <- unique(data.table(chr = g$chr, blk = g$blk))
    bt <- vapply(seq_len(NB), function(b) {
      w <- tabulate(sample.int(D, D, replace = TRUE), D)
      bm <- blocks[, .(blk = sample(blk, .N, replace = TRUE)), by = chr][, .N, by = blk]
      bw <- bm$N[match(g$blk, bm$blk)]; bw[is.na(bw)] <- 0; bw[!base_ok] <- 0
      d <- dC_w(g, w); k <- !is.na(d) & bw > 0
      kk <- !is.na(d) & base_ok                         # donor-only: every pair once
      c(sum(bw[k] * d[k]) / sum(bw[k]), mean(d[kk]))
    }, numeric(2))
    q <- function(v) unname(quantile(v, c(0.025, 0.975)))
    pj <- q(bt[1, ]); bj <- mean(bt[1, ]) - est; pd <- q(bt[2, ]); bd <- mean(bt[2, ]) - est
    c(est = est, lo = pj[1], hi = pj[2], bc_lo = pj[1] - bj, bc_hi = pj[2] - bj,
      n_lo = est - 1.96 * sd(bt[1, ]), n_hi = est + 1.96 * sd(bt[1, ]),
      d_lo = pd[1] - bd, d_hi = pd[2] - bd, pairs = sum(base_ok))
  }, numeric(10)))
  target <- mean(out[, "est"])
  cov_ <- function(l, h) round(mean(out[, l] <= target & target <= out[, h]), 4)
  rbindlist(lapply(list(c("percentile (plan)", "lo", "hi"), c("bias-corrected percentile", "bc_lo", "bc_hi"),
                       c("normal, bootstrap SE", "n_lo", "n_hi"), c("donor-only, bias-corrected", "d_lo", "d_hi")), function(m)
    data.table(scenario = name, interval = m[1], donors = D, replicates = REPS, bootstrap_draws = NB,
               mean_pairs = round(mean(out[, "pairs"])), target_mean_endpoint = signif(target, 4),
               sd_endpoint = signif(sd(out[, "est"]), 4), coverage_95 = cov_(m[2], m[3]),
               coverage_mcse = round(sqrt(0.95 * 0.05 / REPS), 4),
               share_upper_below_0 = round(mean(out[, m[3]] < 0), 4),
               mean_interval_width = signif(mean(out[, m[3]] - out[, m[2]]), 4))))
}
t0 <- Sys.time()
cal <- rbindlist(list(run_scenario("coupled_29", 29L, TRUE, 20261004L),
                      run_scenario("coupled_57", 57L, TRUE, 20261005L),
                      run_scenario("null_57", 57L, FALSE, 20261006L)))
fwrite(cal, file.path(res, "Task57_BootstrapCalibration.csv"))
cat("calibration in", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n"); print(cal)
cat("\ndone\n")
