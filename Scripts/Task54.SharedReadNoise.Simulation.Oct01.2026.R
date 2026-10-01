# ================================================================
# Task 54 - Shared-read noise: the simulation grid (plan phase D)
# Date: Oct 1, 2026
#
# Follows plan_shared_read_noise.md section 7. The 40 scenarios and the
# 9-scenario calibration subset were declared in
# Results/Task54_SimulationScenarios.json before anything ran; the
# script stops if that file has changed (SHA-256 below).
#
# WHAT IS SIMULATED
#   2,000 cohorts per scenario. Each donor has true methylation p_i, p_j
#   at two close CpGs (drawn from the declared latent structure), and a
#   set of DNA fragments: r fragments call both sites (their joint
#   state follows the declared coupling eta), the rest call one site.
#   Fragments are simulated one by one, so the violations (PCR copies,
#   conversion errors, state-dependent call loss, cell mixtures) act
#   on real fragments. Counts are then built exactly as from .pat files.
#
# WHAT IS MEASURED (plan section 7, "Assess")
#   Covariance and squared disagreement, raw and corrected, against
#   the sampled donors' true moments (conditional target) and the
#   population moments, separately: bias with Monte Carlo SE, RMSE.
#   Correlation: errors only where the true correlation is defined,
#   with every failure counted. No-overlap check: with no shared
#   fragments the covariance correction must be exactly 0.
#   Calibration subset: 95% whole-donor percentile bootstrap intervals
#   (2,000 draws) for the population covariance and disagreement:
#   coverage, and false positives where the true covariance is 0.
#
# The fast vectorised estimator here is checked against the tested
# Scripts/SharedReadNoise.Core.R on real simulated cohorts.
#
# OUTPUTS (Results/)
#   Task54_SimulationSummary.csv, Task54_Calibration.csv,
#   Task54_PopulationTargets.csv, Fig40_SimulationBias.png
# ================================================================

suppressPackageStartupMessages({ library(data.table); library(jsonlite); library(matrixStats); library(ggplot2) })

pick_dir <- function(...) { for (d in c(...)) if (dir.exists(d)) return(d)
  stop("no candidate directory exists") }
repo_dir <- pick_dir("/mmfs1/home/wln26/Experiments.Outlier.July31.2026",
                     path.expand("~/Desktop/bioinformatics-research"))
res <- file.path(repo_dir, "Results")
source(file.path(repo_dir, "Scripts", "SharedReadNoise.Core.R"))
SPEC_FILE <- file.path(res, "Task54_SimulationScenarios.json")
SPEC_SHA <- "6033b772555a38f305fc5b067b2fee8125f2164cd1fd44eba44207260eff585a"
stopifnot(digest::digest(file = SPEC_FILE, algo = "sha256") == SPEC_SHA)
spec <- fromJSON(SPEC_FILE)
SC <- as.data.table(spec$scenarios); R <- spec$cohorts_per_scenario; B <- spec$bootstrap_draws
MIN_DON <- spec$min_eligible_donors
cat("declared scenario list verified:", nrow(SC), "scenarios,", R, "cohorts each\n")

# ----------------------------------------------------------------
# latent structures
# ----------------------------------------------------------------
LATENT <- list(
  constant      = list(type = "constant"),
  independent   = list(type = "ln", mu = c(0, 0), sd = 1, rho = 0),
  positive      = list(type = "ln", mu = c(0, 0), sd = 1, rho = 0.7),
  inverse       = list(type = "ln", mu = c(0, 0), sd = 1, rho = -0.7),
  low_spread    = list(type = "ln", mu = c(0, 0), sd = 0.15, rho = 0.7),
  extreme_means = list(type = "ln", mu = c(-3, -3), sd = 0.5, rho = 0.7),
  mean_offset   = list(type = "ln", mu = c(-0.5, 0.5), sd = 1, rho = 0.9))
draw_latent <- function(lat, m) {                     # m donors -> list(pi, pj)
  if (lat$type == "constant") return(list(pi = rep(0.5, m), pj = rep(0.5, m)))
  z1 <- rnorm(m); z2 <- lat$rho * z1 + sqrt(1 - lat$rho^2) * rnorm(m)
  list(pi = plogis(lat$mu[1] + lat$sd * z1), pj = plogis(lat$mu[2] + lat$sd * z2))
}
draw_mixture <- function(m) { phi <- rbeta(m, 2, 5); list(phi = phi, pi = 0.1 + 0.8 * phi, pj = 0.1 + 0.8 * phi) }

set.seed(spec$master_seed)
pop <- rbindlist(lapply(c(names(LATENT), "cell_mixture"), function(nm) {
  d <- if (nm == "cell_mixture") draw_mixture(4e6) else draw_latent(LATENT[[nm]], 4e6)
  data.table(latent = nm, cov = cov(d$pi, d$pj), var_i = var(d$pi), var_j = var(d$pj), d2 = mean((d$pi - d$pj)^2),
             rho = if (var(d$pi) > 0 && var(d$pj) > 0) cor(d$pi, d$pj) else NA_real_)
}))
pop[latent == "constant", `:=`(cov = 0, var_i = 0, var_j = 0, d2 = 0)]
pop[latent == "independent", `:=`(cov = 0, rho = 0)]         # exactly 0 by construction
fwrite(pop, file.path(res, "Task54_PopulationTargets.csv"))
cat("\npopulation targets:\n"); print(pop)

# ----------------------------------------------------------------
# one scenario: fragments -> counts
# ----------------------------------------------------------------
depths <- function(dep, m, p_i, p_j, violation) {
  if (violation == "informative_coverage")
    return(list(Ni = 4L + rpois(m, 12 * p_i), Nj = 4L + rpois(m, 12 * p_j)))
  if (dep == "heterogeneous") { d <- 2L + rnbinom(m, mu = 10, size = 2); return(list(Ni = d, Nj = d)) }
  if (dep == "5/30") return(list(Ni = rep(5L, m), Nj = rep(30L, m)))
  d <- as.integer(dep); list(Ni = rep(d, m), Nj = rep(d, m))
}
simulate <- function(sc, Rn) {
  n <- sc$donors; m <- Rn * n
  if (sc$violation == "cell_mixture") lat <- draw_mixture(m) else lat <- draw_latent(LATENT[[sc$latent]], m)
  pi <- lat$pi; pj <- lat$pj
  dd <- depths(sc$depth, m, pi, pj, sc$violation)
  share <- if (sc$violation == "informative_coverage") 1 else sc$shared
  r <- as.integer(round(share * pmin(dd$Ni, dd$Nj)))
  si <- dd$Ni - r; sj <- dd$Nj - r
  cell <- rep.int(seq_len(m), r + si + sj)
  type <- rep.int(rep(1:3, m), as.vector(rbind(r, si, sj)))
  # per-fragment marginal probabilities (shared fragments may be shifted)
  ps_i <- pi; ps_j <- pj
  if (sc$violation == "nonrepresentative_shared") { ps_i <- plogis(qlogis(pi) + 0.5); ps_j <- plogis(qlogis(pj) + 0.5) }
  eta <- sc$eta
  qb <- if (eta >= 0) pmin(ps_i, ps_j) else pmax(0, ps_i + ps_j - 1)
  q11 <- (1 - abs(eta)) * ps_i * ps_j + abs(eta) * qb
  q10 <- pmax(ps_i - q11, 0); q01 <- pmax(ps_j - q11, 0)
  F <- length(cell); u <- runif(F); ci <- cj <- rep(NA_integer_, F)
  if (sc$violation == "cell_mixture") {
    a <- runif(F) < lat$phi[cell]                       # the fragment's cell type
    pc <- ifelse(a, 0.9, 0.1)
    ci[type != 3L] <- as.integer(runif(sum(type != 3L)) < pc[type != 3L])
    cj[type != 2L] <- as.integer(runif(sum(type != 2L)) < pc[type != 2L])
  } else {
    sh <- type == 1L; c_ <- cell[sh]; us <- u[sh]
    s11 <- us < q11[c_]; s10 <- !s11 & us < q11[c_] + q10[c_]; s01 <- !s11 & !s10 & us < q11[c_] + q10[c_] + q01[c_]
    ci[sh] <- as.integer(s11 | s10); cj[sh] <- as.integer(s11 | s01)
    io <- type == 2L; ci[io] <- as.integer(u[io] < pi[cell[io]])
    jo <- type == 3L; cj[jo] <- as.integer(u[jo] < pj[cell[jo]])
  }
  if (sc$violation == "conversion_errors") {
    fl <- function(x) { k <- !is.na(x); f <- runif(length(x)) < 0.02; x[k & f] <- 1L - x[k & f]; x }
    ci <- fl(ci); cj <- fl(cj)
  }
  if (sc$violation == "state_dependent_loss") {
    lose <- function(x) { k <- !is.na(x); p <- ifelse(x %in% 1L, 0.20, 0.05); x[k & runif(length(x)) < p] <- NA_integer_; x }
    ci <- lose(ci); cj <- lose(cj)
  }
  if (sc$violation == "pcr_duplicates") {
    cp <- 1L + rpois(F, 0.3); cell <- rep.int(cell, cp); ci <- rep.int(ci, cp); cj <- rep.int(cj, cp)
  }
  fr <- data.table(cell, ci, cj)
  agg <- fr[, .(Ni = sum(!is.na(ci)), Mi = sum(ci, na.rm = TRUE), Nj = sum(!is.na(cj)), Mj = sum(cj, na.rm = TRUE),
                n00 = sum(ci %in% 0L & cj %in% 0L), n01 = sum(ci %in% 0L & cj %in% 1L),
                n10 = sum(ci %in% 1L & cj %in% 0L), n11 = sum(ci %in% 1L & cj %in% 1L)), keyby = cell]
  full <- data.table(cell = seq_len(m))[agg, on = "cell", names(agg)[-1] := mget(paste0("i.", names(agg)[-1]))]
  for (v in names(agg)[-1]) set(full, which(is.na(full[[v]])), v, 0L)
  mat <- function(v) matrix(as.numeric(v), Rn, n, byrow = TRUE)
  list(pi = mat(pi), pj = mat(pj), Ni = mat(full$Ni), Mi = mat(full$Mi), Nj = mat(full$Nj), Mj = mat(full$Mj),
       n00 = mat(full$n00), n01 = mat(full$n01), n10 = mat(full$n10), n11 = mat(full$n11), counts = full)
}

# ----------------------------------------------------------------
# the estimator, vectorised over cohorts (rows) and donors (columns)
# ----------------------------------------------------------------
estimate <- function(s) {
  r <- s$n00 + s$n01 + s$n10 + s$n11
  E <- s$Ni >= 2 & s$Nj >= 2 & r != 1
  x <- ifelse(E, s$Mi / s$Ni, 0); y <- ifelse(E, s$Mj / s$Nj, 0)
  vx <- ifelse(E, x * (1 - x) / (s$Ni - 1), 0); vy <- ifelse(E, y * (1 - y) / (s$Nj - 1), 0)
  cc <- ifelse(E & r >= 2, (r * s$n11 - (s$n10 + s$n11) * (s$n01 + s$n11)) / (s$Ni * s$Nj * pmax(r - 1, 1)), 0)
  ne <- rowSums(E); ok <- ne >= MIN_DON
  mv <- function(a, b) (rowSums(E * a * b) - rowSums(E * a) * rowSums(E * b) / ne) / (ne - 1)
  pi <- ifelse(E, s$pi, 0); pj <- ifelse(E, s$pj, 0)
  o <- data.table(ok, ne, unsupported_r1 = rowSums(s$Ni >= 2 & s$Nj >= 2 & r == 1), low_cov = rowSums(s$Ni < 2 | s$Nj < 2),
                  raw_cov = mv(x, y), raw_vi = mv(x, x), raw_vj = mv(y, y),
                  mc = rowSums(E * cc) / ne, mvx = rowSums(E * vx) / ne, mvy = rowSums(E * vy) / ne,
                  raw_d2 = rowSums(E * (x - y)^2) / ne,
                  T = mv(pi, pj), Tii = mv(pi, pi), Tjj = mv(pj, pj), Q = rowSums(E * (pi - pj)^2) / ne)
  o[, `:=`(cor_cov = raw_cov - mc, cor_vi = raw_vi - mvx, cor_vj = raw_vj - mvy,
           cor_d2 = raw_d2 - mvx - mvy + 2 * mc)]
  o[, `:=`(raw_rho = fifelse(raw_vi > 0 & raw_vj > 0, raw_cov / sqrt(pmax(raw_vi * raw_vj, 0)), NA_real_),
           cor_rho = fifelse(cor_vi > 0 & cor_vj > 0, cor_cov / sqrt(pmax(cor_vi * cor_vj, 0)), NA_real_),
           T_rho = fifelse(Tii > 0 & Tjj > 0, T / sqrt(pmax(Tii * Tjj, 0)), NA_real_))]
  o[, status := fcase(is.na(cor_rho), "nonpositive_variance", abs(cor_rho) > 1, "outside_unit_interval", default = "valid")]
  o
}

# check the fast estimator against the tested core code on real cohorts
check_core <- function(sc) {
  set.seed(7); s <- simulate(sc, 3L); o <- estimate(s)
  for (k in 1:3) {
    d <- data.table(Ni = s$Ni[k, ], Mi = s$Mi[k, ], Nj = s$Nj[k, ], Mj = s$Mj[k, ],
                    n00 = s$n00[k, ], n01 = s$n01[k, ], n10 = s$n10[k, ], n11 = s$n11[k, ])
    acc <- srn_accumulator(1L)
    for (q in seq_len(nrow(d))) srn_update(acc, d[q], 1L)
    f <- srn_finish(acc, min_donors = MIN_DON)
    if (!o$ok[k]) { stopifnot(nrow(f) == 0); next }
    stopifnot(isTRUE(all.equal(c(f$raw_cov, f$cov_corrected, f$var_i_corrected, f$d2_corrected),
                               c(o$raw_cov[k], o$cor_cov[k], o$cor_vi[k], o$cor_d2[k]), tolerance = 1e-10)))
  }
  TRUE
}
for (nm in c("base", "eta_-1.0", "depth_heterogeneous", "violation_pcr_duplicates")) stopifnot(check_core(SC[scenario == nm]))
cat("\nfast estimator matches SharedReadNoise.Core.R on simulated cohorts\n")

# ----------------------------------------------------------------
# run every declared scenario
# ----------------------------------------------------------------
mse <- function(a) mean(a^2); se <- function(a) sd(a) / sqrt(length(a))
t0 <- Sys.time(); keep <- list()
summ <- rbindlist(lapply(seq_len(nrow(SC)), function(i) {
  sc <- SC[i]; set.seed(spec$master_seed + i)
  s <- simulate(sc, R); o <- estimate(s)
  pp <- pop[latent == (if (sc$violation == "cell_mixture") "cell_mixture" else sc$latent)]
  g <- o[ok == TRUE]
  if (sc$scenario %in% spec$calibration_subset) keep[[sc$scenario]] <<- list(s = s, o = o, pp = pp)
  rho_def <- g[!is.na(T_rho)]
  out <- data.table(scenario = sc$scenario, group = sc$group, cohorts_ok = nrow(g),
    donor_rows_unsupported_r1 = round(mean(o$unsupported_r1) / sc$donors, 4),
    donor_rows_low_coverage = round(mean(o$low_cov) / sc$donors, 4),
    true_cov_population = signif(pp$cov, 4),
    bias_raw_cov_cond = signif(mean(g$raw_cov - g$T), 4), bias_cor_cov_cond = signif(mean(g$cor_cov - g$T), 4),
    mcse_cor_cov_cond = signif(se(g$cor_cov - g$T), 3),
    rmse_raw_cov_cond = signif(sqrt(mse(g$raw_cov - g$T)), 4), rmse_cor_cov_cond = signif(sqrt(mse(g$cor_cov - g$T)), 4),
    bias_raw_cov_pop = signif(mean(g$raw_cov - pp$cov), 4), bias_cor_cov_pop = signif(mean(g$cor_cov - pp$cov), 4),
    mcse_cor_cov_pop = signif(se(g$cor_cov - pp$cov), 3),
    rmse_raw_cov_pop = signif(sqrt(mse(g$raw_cov - pp$cov)), 4), rmse_cor_cov_pop = signif(sqrt(mse(g$cor_cov - pp$cov)), 4),
    bias_raw_d2_cond = signif(mean(g$raw_d2 - g$Q), 4), bias_cor_d2_cond = signif(mean(g$cor_d2 - g$Q), 4),
    rmse_raw_d2_cond = signif(sqrt(mse(g$raw_d2 - g$Q)), 4), rmse_cor_d2_cond = signif(sqrt(mse(g$cor_d2 - g$Q)), 4),
    share_cor_d2_negative = round(mean(g$cor_d2 < 0), 4),
    share_rho_valid = round(mean(g$status == "valid"), 4), share_rho_outside = round(mean(g$status == "outside_unit_interval"), 4),
    share_rho_nonpositive_var = round(mean(g$status == "nonpositive_variance"), 4),
    rmse_raw_rho_cond = if (nrow(rho_def)) signif(sqrt(mean((rho_def$raw_rho - rho_def$T_rho)^2, na.rm = TRUE)), 4) else NA_real_,
    rmse_cor_rho_cond_valid = if (nrow(rho_def[status == "valid"])) signif(sqrt(mean((rho_def[status == "valid"]$cor_rho - rho_def[status == "valid"]$T_rho)^2)), 4) else NA_real_,
    max_abs_noise_cov = signif(max(abs(g$mc)), 4))
  cat(sprintf("  %2d/%d %-34s ok %4d  cov bias raw %+.2e  corrected %+.2e (%.1f min)\n", i, nrow(SC), sc$scenario,
              nrow(g), out$bias_raw_cov_cond, out$bias_cor_cov_cond, as.numeric(Sys.time() - t0, units = "mins")))
  out
}))
fwrite(summ, file.path(res, "Task54_SimulationSummary.csv"))
nov <- summ[scenario == "shared_0.0"]
stopifnot(nrow(nov) == 1)
cat("\nno-overlap check (shared = 0): largest |mean noise covariance| =", nov$max_abs_noise_cov, "\n")
stopifnot(nov$max_abs_noise_cov == 0)

# ----------------------------------------------------------------
# calibration subset: whole-donor percentile bootstrap
# ----------------------------------------------------------------
boot_one <- function(s, o, k) {
  r <- s$n00[k, ] + s$n01[k, ] + s$n10[k, ] + s$n11[k, ]
  e <- s$Ni[k, ] >= 2 & s$Nj[k, ] >= 2 & r != 1
  x <- (s$Mi[k, ] / s$Ni[k, ])[e]; y <- (s$Mj[k, ] / s$Nj[k, ])[e]
  Ni <- s$Ni[k, e]; Nj <- s$Nj[k, e]; re <- r[e]
  n11 <- s$n11[k, e]; n10 <- s$n10[k, e]; n01 <- s$n01[k, e]
  cc <- ifelse(re >= 2, (re * n11 - (n10 + n11) * (n01 + n11)) / (Ni * Nj * pmax(re - 1, 1)), 0)
  vx <- x * (1 - x) / (Ni - 1); vy <- y * (1 - y) / (Nj - 1)
  m <- length(x); id <- matrix(sample.int(m, B * m, replace = TRUE), B, m)
  S <- function(v) rowSums(matrix(v[id], B, m))
  sx <- S(x); sy <- S(y); rc <- (S(x * y) - sx * sy / m) / (m - 1)
  mc <- S(cc) / m; d2 <- S((x - y)^2) / m; cd2 <- d2 - S(vx) / m - S(vy) / m + 2 * mc
  q <- function(v) unname(quantile(v, c(0.025, 0.975)))
  c(q(rc), q(rc - mc), q(d2), q(cd2))
}
calib <- rbindlist(lapply(spec$calibration_subset, function(nm) {
  z <- keep[[nm]]; set.seed(spec$master_seed + 1000L + match(nm, spec$calibration_subset))
  okr <- which(z$o$ok)
  iv <- t(vapply(okr, function(k) boot_one(z$s, z$o, k), numeric(8)))
  cov_in <- function(lo, hi, t) mean(lo <= t & t <= hi)
  bin_se <- function(p) sqrt(p * (1 - p) / length(okr))
  cr <- cov_in(iv[, 1], iv[, 2], z$pp$cov); cc_ <- cov_in(iv[, 3], iv[, 4], z$pp$cov)
  dr <- cov_in(iv[, 5], iv[, 6], z$pp$d2); dc <- cov_in(iv[, 7], iv[, 8], z$pp$d2)
  data.table(scenario = nm, cohorts = length(okr), true_cov = signif(z$pp$cov, 4),
             coverage_raw_cov = round(cr, 4), coverage_cor_cov = round(cc_, 4), mcse_cor_cov = round(bin_se(cc_), 4),
             coverage_raw_d2 = round(dr, 4), coverage_cor_d2 = round(dc, 4), mcse_cor_d2 = round(bin_se(dc), 4),
             false_pos_raw_cov = if (z$pp$cov == 0) round(mean(iv[, 1] > 0 | iv[, 2] < 0), 4) else NA_real_,
             false_pos_cor_cov = if (z$pp$cov == 0) round(mean(iv[, 3] > 0 | iv[, 4] < 0), 4) else NA_real_)
}))
fwrite(calib, file.path(res, "Task54_Calibration.csv"))
cat("\nCalibration (95% whole-donor bootstrap intervals for population targets):\n"); print(calib)
cat("\nSummary (bias against the sampled donors' true covariance):\n")
print(summ[, .(scenario, cohorts_ok, bias_raw_cov_cond, bias_cor_cov_cond, mcse_cor_cov_cond,
               rmse_raw_cov_cond, rmse_cor_cov_cond, share_rho_valid)], nrows = 60)

# ----------------------------------------------------------------
# figure: covariance bias, raw and corrected, every scenario
# ----------------------------------------------------------------
SURF <- "#fcfcfb"; INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"; GRID <- "#e1e0d9"
f40 <- melt(summ[, .(scenario = factor(scenario, levels = rev(scenario)), group,
                     Observed = bias_raw_cov_cond, Corrected = bias_cor_cov_cond)],
            id.vars = c("scenario", "group"), variable.name = "estimator", value.name = "bias")
f40[, group := factor(group, levels = c("core", "interaction", "violation"),
                      labels = c("Core grid", "Interactions", "Assumption violations"))]
p40 <- ggplot(f40, aes(bias, scenario, colour = estimator)) +
  geom_vline(xintercept = 0, colour = MUTED, linewidth = 0.4) +
  geom_point(size = 2.4) +
  scale_colour_manual(values = c(Observed = "#52514e", Corrected = "#2a78d6"), name = NULL) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  labs(title = "Simulation: bias of the covariance across donors, observed and corrected",
       subtitle = paste0("2,000 cohorts per scenario, bias against the sampled donors' true covariance. ",
                         "Scenarios declared before the run (Task54_SimulationScenarios.json)."),
       x = "Mean error (estimate minus true covariance)", y = NULL,
       caption = "Source: Results/Task54_SimulationSummary.csv (Task 54)") +
  theme_minimal(base_size = 10) +
  theme(plot.background = element_rect(fill = SURF, colour = NA), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = GRID, linewidth = 0.3),
        strip.text.y = element_text(angle = 0, hjust = 0, face = "bold", colour = INK),
        plot.title = element_text(colour = INK, face = "bold", size = 12.5), plot.title.position = "plot",
        plot.subtitle = element_text(colour = INK2, size = 9.5), plot.caption = element_text(colour = MUTED, hjust = 0),
        plot.caption.position = "plot", legend.position = "top", legend.justification = "left",
        axis.text = element_text(colour = INK2))
ggsave(file.path(res, "Fig40_SimulationBias.png"), p40, width = 10, height = 10, dpi = 200, bg = SURF)
cat("\ndone\n")
