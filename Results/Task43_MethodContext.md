# Task 43: established dependence estimators and a defensible benchmark

Reviewed 2026-09-29. This note informs a controlled simulation at 18 patients; it does not report simulation results or claim a new coefficient.

**Established methods already handle discrete, tied dependence.** Distance covariance is a suitable additional comparator. For bounded decile variables its moment assumptions hold automatically; ties do not require random breaking. The population coefficient is zero exactly under independence. Ordinary sample distance correlation is a nonnegative descriptive statistic and can be appreciably positive under the finite-sample null; its numerical magnitude alone is not a significance test. [Székely, Rizzo and Bakirov, *Measuring and testing dependence by correlation of distances*, Annals of Statistics 35 (2007), 2769–2794](https://doi.org/10.1214/009053607000000505).

**Distinguish ordinary distance correlation from unbiased squared distance covariance.** For `a_ij=|x_i-x_j|`, U-centering sets the diagonal to zero and, off diagonal, uses

```text
A_ij = a_ij - row_sum_i/(n-2) - col_sum_j/(n-2)
             + total_sum/((n-1)*(n-2)).
U_dcov2 = sum_(i != j) A_ij B_ij / (n*(n-3)).
```

For iid observations and `n>3`, `U_dcov2` is unbiased for population squared distance covariance. It can be negative; this does not indicate inverse association. The normalized ratio is not thereby an unbiased estimator of population distance correlation. Crucially, a U-centered matrix vanishes when at least `n-1` sample values are identical. Thus 17 identical values plus one distinct value at `n=18` gives zero U-distance variance even though ordinary distance correlation is defined. Retain an explicit degeneracy flag. These properties and the estimator appear in Definition 2 and Proposition 1 of [Székely and Rizzo, *Partial Distance Correlation with Methods for Dissimilarities*, Annals of Statistics 42 (2014), 2382–2412](https://doi.org/10.1214/14-AOS1255); [author manuscript](https://arxiv.org/abs/1310.2926).

**Xi is a complementary estimator, not a guaranteed small-sample improvement.** Published analyses establish deficient local power against certain alternatives relative to established rank tests. Those asymptotic comparisons do not prove the ordering of methods for the particular tied distributions at `n=18`; the proposed simulation can assess that restricted question. [Shi, Drton and Han, *On the power of Chatterjee's rank correlation*, Biometrika 109 (2022), 317–333](https://doi.org/10.1093/biomet/asab028).

**Finite-sample permutation validity does not require a large n.** Permutation inference requires exchangeability under the chosen null. Use the same patient-label permutations for every method, preserve ties, count equality in the upper tail, and use `(1 + exceedances)/(B + 1)` for independent uniform Monte Carlo permutations including the observed configuration through the correction. Complete enumeration uses its exact tail proportion instead. Never report zero p-values. A small permutation budget limits p-value resolution and cannot support arbitrary genome-wide significance thresholds. [Phipson and Smyth, *Permutation P-values Should Never Be Zero: Calculating Exact P-values When Permutations Are Randomly Drawn*, Statistical Applications in Genetics and Molecular Biology 9 (2010), article 39](https://doi.org/10.2202/1544-6115.1585); [author manuscript](https://gksmyth.github.io/pubs/PermPValuesPreprint.pdf).

Practical benchmark decisions, based on these properties:

- Compare calibrated rejection rates under independence before comparing power. Report Monte Carlo uncertainty and the rate of undefined/degenerate scores, not just the average score. Use two-sided absolute Pearson/Spearman for general association and upper-tail xi/distance statistics. Keep agreement performance as a separate target from dependence detection.
- Include both ordinary distance correlation and `U_dcov2`; do not use the absolute value of `U_dcov2` as a dependence test. Retain negative values rather than clipping them before calibration. Normalization by positive marginal distance variances cannot improve a given pair's permutation ranking because those variances remain fixed across permutations.
- Include balanced and highly imbalanced marginals, constants, singleton changes, monotone/inverse patterns, nonlinear patterns, shifted supports, and bin-boundary perturbations. Compare methods on the same representation; a raw-beta versus decile comparison measures information loss as well as estimator choice.
- The 18-patient rare-high example is a useful resolution limit: if both vectors contain one high and 17 low values, alignment of the two highs has probability `1/18` under random pairing. Even a maximal observed association then has one-sided exact p-value `1/18`, above 0.05. This is a limitation of that marginal configuration, not proof that all 18-patient analyses are impossible.

**What the benchmark can establish.** It can select a useful established score for specified RRBS-like scenarios and identify a reproducible failure requiring further work. It cannot establish novelty, universal superiority, or biological validation. A new formula would need a precise unmet target, a demonstrated distinction from existing estimators, justified mathematical properties, and evaluation on separately specified scenarios and independent data. If an existing estimator performs adequately, adopting it is a successful outcome for item 5. Rescaling weighted kappa does not create additional permutation information.
