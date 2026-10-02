# Fixed n=18 benchmark following the Task42 candidate rejection.
# Synthetic beta ratios, not sequencing reads or inferred biological truth.
# Run: Rscript Scripts/Task43.Item5.ControlledBenchmark.Sep29.2026.R
# Optional smoke run (separate filenames): --smoke
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- dirname(dirname(normalizePath(sub("^--file=", "", file_arg[1L]))))
source(file.path(root, "Scripts", "Task41.RRBS.Item5.Feasibility.Sep28.2026.R"))

distance_kernel <- function(a, unbiased = FALSE) {
  n <- ncol(a)
  t(vapply(seq_len(nrow(a)), function(i) {
    z <- abs(outer(a[i, ], a[i, ], "-"))
    if (unbiased) {
      z <- z - rowSums(z)/(n-2) - rep(colSums(z)/(n-2), each = n) +
        sum(z)/((n-1)*(n-2))
      diag(z) <- 0
    } else {
      z <- z - rowMeans(z) - rep(colMeans(z), each = n) + mean(z)
    }
    as.vector(z)
  }, numeric(n*n)))
}

kernel_cor <- function(a, b, den, unbiased = FALSE) {
  z <- rowSums(a*b)
  if (unbiased) return(z/(18*15)) # unnormalized unbiased squared dCov
  ans <- sqrt(pmax(z, 0)/den)
  ans[den == 0] <- NA_real_
  ans
}

methods <- c("ordinal_similarity", "linear_kappa", "quadratic_kappa", "pearson_beta",
             "pearson_levels", "spearman_beta", "spearman_levels", "xi_levels_mean",
             "dcor_beta", "dcor_levels", "unbiased_dcov2_levels")
agreement_methods <- c("ordinal_similarity", "linear_kappa", "quadratic_kappa")

prepare_pair <- function(a, b) {
  x <- levels_beta(a); y <- levels_beta(b)
  px <- sapply(1:10, function(k) rowMeans(x == k))
  py <- sapply(1:10, function(k) rowMeans(y == k))
  d0 <- rowSums((px %*% abs(outer(1:10, 1:10, "-"))) * py)
  ord_exp <- rowMeans(x^2) + rowMeans(y^2) - 2*rowMeans(x)*rowMeans(y)
  kxa <- distance_kernel(a); kya <- distance_kernel(b)
  kxl <- distance_kernel(x); kyl <- distance_kernel(y)
  kxu <- distance_kernel(x, TRUE); kyu <- distance_kernel(y, TRUE)
  list(a = a, b = b, x = x, y = y, d0 = d0, ord_exp = ord_exp,
       ra = rowRanks(a, ties.method = "average", preserveShape = TRUE),
       rb = rowRanks(b, ties.method = "average", preserveShape = TRUE),
       rx = rowRanks(x, ties.method = "average", preserveShape = TRUE),
       ry = rowRanks(y, ties.method = "average", preserveShape = TRUE),
       kxa = kxa, kya = kya, kxl = kxl, kyl = kyl, kxu = kxu, kyu = kyu,
       dena = sqrt(rowSums(kxa^2)*rowSums(kya^2)),
       denl = sqrt(rowSums(kxl^2)*rowSums(kyl^2)))
}

score_prepared <- function(q, perms = NULL) {
  nr <- nrow(q$a); n <- ncol(q$a)
  b <- q$b; y <- q$y; rb <- q$rb; ry <- q$ry
  kya <- q$kya; kyl <- q$kyl; kyu <- q$kyu
  if (!is.null(perms)) {
    ix <- cbind(rep(seq_len(nr), each = n), as.vector(t(perms)))
    reorder <- function(v) matrix(v[ix], nr, n, byrow = TRUE)
    b <- reorder(b); y <- reorder(y); rb <- reorder(rb); ry <- reorder(ry)
    # Every replicate receives its own independent patient permutation.
    # All methods use those same permutations for fair Monte Carlo comparison.
    ki <- t(vapply(seq_len(nr), function(i) {
      as.vector(outer(perms[i, ], (perms[i, ] - 1L)*n, "+"))
    }, integer(n*n)))
    ixk <- cbind(rep(seq_len(nr), each = n*n), as.vector(t(ki)))
    reorder_k <- function(v) matrix(v[ixk], nr, n*n, byrow = TRUE)
    kya <- reorder_k(kya); kyl <- reorder_k(kyl); kyu <- reorder_k(kyu)
  }
  d <- rowMeans(abs(q$x-y)); dl2 <- rowMeans((q$x-y)^2)
  xi1 <- vapply(seq_len(nr), function(i) xi_tie_average(q$x[i, ], y[i, ]), 0.0)
  xi2 <- vapply(seq_len(nr), function(i) xi_tie_average(y[i, ], q$x[i, ]), 0.0)
  ans <- cbind(ordinal_similarity = 1-d/9,
    linear_kappa = ifelse(q$d0 > 0, 1-d/q$d0, NA_real_),
    quadratic_kappa = ifelse(q$ord_exp > 0, 1-dl2/q$ord_exp, NA_real_),
    pearson_beta = abs(row_correlation(q$a, b)),
    pearson_levels = abs(row_correlation(q$x, y)),
    spearman_beta = abs(row_correlation(q$ra, rb)),
    spearman_levels = abs(row_correlation(q$rx, ry)),
    xi_levels_mean = (xi1+xi2)/2,
    dcor_beta = kernel_cor(q$kxa, kya, q$dena),
    dcor_levels = kernel_cor(q$kxl, kyl, q$denl),
    unbiased_dcov2_levels = kernel_cor(q$kxu, kyu, q$denl, TRUE))
  stopifnot(identical(colnames(ans), methods))
  ans
}

# Distances checked against independent explicit row/column formulas.
test_a <- rbind(seq(0, 1, length.out = 18), rep(.1, 18), c(1, rep(0, 17)))
ka <- distance_kernel(test_a); ku <- distance_kernel(test_a, TRUE)
stopifnot(abs(kernel_cor(ka, ka, rowSums(ka^2))[1] - 1) < 1e-12,
          is.na(kernel_cor(ka, ka, rowSums(ka^2))[2]),
          max(abs(ku[3, ])) < 1e-12)
for (i in 1:3) {
  z <- abs(outer(test_a[i, ], test_a[i, ], "-")); n <- 18L
  for (j in 1:n) for (k in 1:n) {
    expected <- if (j == k) 0 else z[j,k] - sum(z[j,])/(n-2) - sum(z[,k])/(n-2) + sum(z)/((n-1)*(n-2))
    stopifnot(abs(ku[i, j+(k-1)*n]-expected) < 1e-12)
  }
}
set.seed(43001)
ta <- matrix(runif(3*18), 3); tb <- matrix(runif(3*18), 3)
tp <- t(replicate(3, sample.int(18)))
tbp <- t(vapply(1:3, function(i) tb[i, tp[i,]], numeric(18)))
cached <- score_prepared(prepare_pair(ta, tb), tp)
fresh <- score_prepared(prepare_pair(ta, tbp))
stopifnot(max(abs(cached-fresh), na.rm = TRUE) < 1e-12,
          identical(is.na(cached), is.na(fresh)))
cat("PASS: distance kernels, U-centering, singleton degeneracy, independent row permutations\n")

smoke <- "--smoke" %in% commandArgs(TRUE)
nr <- if (smoke) 30L else 400L
B <- if (smoke) 19L else 199L
n <- 18L; alpha <- .05; seed <- 20260929L
prefix <- if (smoke) "Task43_Smoke_" else "Task43_"
out <- file.path(root, "Results")
write_tab <- function(x, name) fwrite(x, file.path(out, paste0(prefix, name, ".csv")))
clamp <- function(x, lo=0, hi=1) pmin(pmax(x, lo), hi)
scenario <- data.table(
  scenario = c("independent_continuous", "independent_sparse", "independent_bin_boundary",
               "positive_linear", "offset_linear", "inverse_linear", "nonlinear_u",
               "within_low_bin", "shared_patient_factor"),
  truth = c(rep("independent_null", 3), rep("dependent", 5), "dependent_common_cause"),
  description = c("Independent Uniform(.05,.95) beta ratios",
    "Independent .05/.95 ratios, high probability .10 at each site",
    "Independent N(.10,.001^2), clipped to [0,1]",
    "X Uniform(.10,.90); Y=.1+.8X+N(0,.08^2), clipped",
    "X Uniform(.05,.25); Y=.65+X+N(0,.02^2), clipped",
    "X Uniform(.10,.90); Y=1-X+N(0,.08^2), clipped",
    "X Uniform(.02,.98); Y=4(X-.5)^2+N(0,.06^2), clipped",
    "X Uniform(.005,.085); Y=X+N(0,.008^2), clipped to [0,.099]",
    "Independent errors: X=.5+.25Z+e1, Y=.5+.25Z+e2; Z~N(0,1), e~N(0,.04^2), clipped"))
write_tab(scenario, "ProtocolScenarios")
jsonlite::write_json(list(seed=seed, n=n, replicates_per_scenario=nr, permutations=B,
  alpha=alpha, frozen_before_simulation=TRUE,
  denominators="Rejections/all generated replicates and rejections/defined separately",
  methods=methods, agreement_target="Positive paired similarity, one-sided",
  dependence_target="Any dependence; Pearson/Spearman absolute, remaining statistics upper-tail",
  U_statistic="Unbiased squared dCov, not normalized, negative values retained; degenerate zero is not omitted",
  uncertainty="Exact binomial 95% intervals across independent synthetic replicates",
  no_tuning="No fitted method, threshold optimization, external cohort access or invented read depths",
  script_md5=unname(tools::md5sum(sub("^--file=", "", file_arg[1L])))),
  file.path(out, paste0(prefix, "Protocol.json")), pretty=TRUE, auto_unbox=TRUE)

generate_case <- function(name) {
  draw <- function(expr) matrix(expr, nr, n)
  if (name == "independent_continuous") {
    a <- draw(runif(nr*n,.05,.95)); b <- draw(runif(nr*n,.05,.95))
  } else if (name == "independent_sparse") {
    a <- .05+.9*draw(rbinom(nr*n,1,.1)); b <- .05+.9*draw(rbinom(nr*n,1,.1))
  } else if (name == "independent_bin_boundary") {
    a <- clamp(draw(rnorm(nr*n,.10,.001))); b <- clamp(draw(rnorm(nr*n,.10,.001)))
  } else if (name == "shared_patient_factor") {
    z <- draw(rnorm(nr*n)); a <- clamp(.5+.25*z+draw(rnorm(nr*n,0,.04)))
    b <- clamp(.5+.25*z+draw(rnorm(nr*n,0,.04)))
  } else if (name == "within_low_bin") {
    a <- draw(runif(nr*n,.005,.085)); b <- clamp(a+draw(rnorm(nr*n,0,.008)),0,.099)
  } else if (name == "offset_linear") {
    a <- draw(runif(nr*n,.05,.25)); b <- clamp(.65+a+draw(rnorm(nr*n,0,.02)))
  } else if (name == "nonlinear_u") {
    a <- draw(runif(nr*n,.02,.98)); b <- clamp(4*(a-.5)^2+draw(rnorm(nr*n,0,.06)))
  } else {
    a <- draw(runif(nr*n,.1,.9))
    if (name == "positive_linear") b <- clamp(.1+.8*a+draw(rnorm(nr*n,0,.08)))
    else b <- clamp(1-a+draw(rnorm(nr*n,0,.08)))
  }
  list(a=a,b=b)
}

all_results <- list(); all_summary <- list(); started <- proc.time()[3]
for (s in seq_len(nrow(scenario))) {
  set.seed(seed+s*1000L)
  name <- scenario$scenario[s]
  data <- generate_case(name); q <- prepare_pair(data$a, data$b)
  obs <- score_prepared(q); defined <- is.finite(obs)
  exceed <- matrix(0L, nr, length(methods))
  cat("Scenario", s, "of", nrow(scenario), ":", name, "\n")
  for (p in seq_len(B)) {
    perms <- t(replicate(nr, sample.int(n)))
    values <- score_prepared(q, perms)
    stopifnot(identical(is.finite(values), defined))
    # Scaled tolerance handles small U statistics without rounding to zero.
    tol <- 1e-12*pmax(1, abs(obs))
    comparison <- values >= obs - tol
    comparison[!defined] <- FALSE
    exceed <- exceed + comparison
  }
  pval <- (1+exceed)/(B+1); pval[!defined] <- NA_real_
  result <- rbindlist(lapply(seq_along(methods), function(j) data.table(
    scenario=name, truth=scenario$truth[s], replicate=seq_len(nr), method=methods[j],
    statistic=obs[,j], defined=defined[,j], permutation_p=pval[,j],
    rejected=is.finite(pval[,j]) & pval[,j] <= alpha)))
  summary <- result[, {
    nd <- sum(defined); hits <- sum(rejected); ci_all <- binom.test(hits,.N)$conf.int
    ci_def <- if (nd) binom.test(hits,nd)$conf.int else c(NA_real_,NA_real_)
    list(generated=.N, defined=nd, rejections=hits, evaluable_share=nd/.N,
         operational_rejection_rate=hits/.N, operational_ci_low=ci_all[1], operational_ci_high=ci_all[2],
         rejection_rate_given_defined=if(nd) hits/nd else NA_real_,
         conditional_ci_low=ci_def[1], conditional_ci_high=ci_def[2])
  }, by=.(scenario,truth,method)]
  summary[, target := fifelse(method %in% agreement_methods, "positive_similarity", "any_dependence")]
  all_results[[s]] <- result; all_summary[[s]] <- summary
  write_tab(rbindlist(all_summary), "BenchmarkSummary")
  cat("Completed", name, "; elapsed", round(proc.time()[3]-started,1), "seconds\n")
}
write_tab(rbindlist(all_results), "BenchmarkReplicates")
summary <- rbindlist(all_summary)

# Sparse coincidence: finite sample evidence cannot be inferred from score=1.
stress <- data.table(n=18, rare_patients=c(1,2),
  exact_one_sided_overlap_p=c(1/choose(18,1),1/choose(18,2)),
  identical_profiles_dcor=c(1,1),
  uncentered_note=c("Perfect score, p=1/18; U-centered distance kernel is zero",
                   "Two shared rare samples; p=1/153 for maximum-overlap statistic"))
write_tab(stress, "SparseEvidence")

if (!smoke) {
  suppressPackageStartupMessages(library(ggplot2))
  shown <- summary[method %in% c("pearson_beta","spearman_levels","xi_levels_mean","dcor_beta",
                                "dcor_levels","ordinal_similarity")]
  scenario_order <- scenario$scenario
  shown[, scenario := factor(scenario, levels=scenario_order)]
  shown[, method := factor(method, levels=c("pearson_beta","spearman_levels","xi_levels_mean",
    "dcor_beta","dcor_levels","ordinal_similarity"), labels=c("Pearson beta","Spearman levels",
    "Xi levels","dCor beta","dCor levels","Ordinal agreement"))]
  g <- ggplot(shown, aes(method,operational_rejection_rate,color=method))+
    geom_hline(yintercept=.05,linetype="dashed",color="grey60")+
    geom_errorbar(aes(ymin=operational_ci_low,ymax=operational_ci_high),width=.2)+geom_point(size=2)+
    facet_wrap(~scenario,ncol=3)+scale_y_continuous(limits=c(0,1),breaks=c(0,.5,1))+
    labs(title="Item 5: what 18 samples can detect",subtitle="400 synthetic replicates per scenario; 199 independent permutations each",
      x=NULL,y="Rejections / all replicates (95% binomial interval)",
      caption="First three panels are independence nulls. Undefined scores count as no call.\nOrdinal agreement tests positive similarity; other shown methods test dependence. Common cause is real dependence, not a false positive.")+
    theme_bw(base_size=10)+theme(legend.position="none",axis.text.x=element_text(angle=55,hjust=1),
      plot.title=element_text(face="bold"),strip.background=element_rect(fill="#eef3f8"))
  ggsave(file.path(out,"Task43_ControlledBenchmark.png"),g,width=12,height=9,dpi=180)
}
capture.output(sessionInfo(),file=file.path(out,paste0(prefix,"SessionInfo.txt")))
cat("Completed benchmark in",round(proc.time()[3]-started,1),"seconds\n")
