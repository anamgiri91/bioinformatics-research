# Independent review for the n = 18 benchmark and report section 5

Reviewed 2026-09-29. Task 41 is the feasibility analysis and Task 42 is the tested ordinal candidate. This note proposes checks and interpretations; it contains no invented benchmark results. Existing scripts and `report.md` were inspected without editing them.

## Minimal fixed simulation design

Generate independent patients within every stochastic replicate, and independently generate replicates. Keep n = 18, ten-bin boundaries, noise scales, alpha, and permutation count fixed before reviewing results. The following six scenarios are sufficient to cover distinct questions:

| Scenario | Suggested fixed data-generating mechanism | Interpretation |
|---|---|---|
| Continuous independence | X and Y independently uniform on [0.05, 0.95] | Genuine independence null with broad support. |
| Sparse independence | X = 0.05 + 0.9 A, Y = 0.05 + 0.9 B; A and B independently Bernoulli(0.1) | Genuine independence null with ties, constant samples, and coarse permutation p-values. |
| Positive monotone | X uniform on [0.05, 0.95]; Y = clip(0.1 + 0.8 X + 0.06 epsilon), epsilon standard normal | Positive association; appropriate signal for one-sided agreement tests and two-sided correlations. |
| Inverse monotone | Same X; Y = clip(0.9 - 0.8 X + 0.06 epsilon) | Dependence with disagreement. One-sided similarity tests target a different alternative. |
| Nonlinear U shape | Same X; Y = clip(0.05 + 0.85 (2 X - 1)^2 + 0.04 epsilon) | Nonmonotone dependence. Finite-n power must be measured, not inferred from xi's asymptotic consistency. |
| Shared patient factor | Z independently Bernoulli(0.5) across patients; X = clip(0.15 + 0.7 Z + 0.04 epsilon_x), Y similarly with independent epsilon_y | True marginal dependence but conditional independence given Z. This is a confounder illustration, never an independence-null false-positive scenario. |

Here `clip` truncates to [0,1]; all noise terms and latent draws are independent unless their sharing is stated. These numbers are suggestions for a fixed design, not parameters selected to favor a score. If an existing implementation has already fixed other sensible values, retain those rather than retuning against outcomes.

Add identical constants, different constants, and one shared high patient as deterministic evaluability/influence examples, outside the null-size and power denominators. A one-high profile can also be generated stochastically, but conditioning on exactly one high is a separate stress design and must be labelled.

## Calibration and fair comparisons

- Report generated replicates, defined replicates, undefined share, rejection rate among defined replicates, and rejection count divided by all replicates. The latter is operational yield, not calibration conditional on being evaluable. Constants are not failed or nonsignificant numerical correlations; their statistic is undefined.
- Use the plus-one Monte Carlo permutation p-value and inclusive exceedance for ties. Keep the same permutation for all statistics within a replicate to support paired comparisons.
- State that Pearson/Spearman are two-sided while similarity, kappa, and xi are one-sided. A low similarity-test rejection rate for inverse dependence is not proof that the method is defective; it answers a different question.
- Raw Pearson versus discretized xi/Spearman mixes the effect of quantization with the effect of the statistic. A Pearson-on-levels comparator is inexpensive and makes this explicit. Quadratic kappa is a positive marginal rescaling of Pearson on those same levels, so their one-sided fixed-margin permutation ordering agrees.
- Sparse-null tests can be conservative because of discreteness. If both profiles have one high among 18 patients, the smallest exact one-sided p-value is 1/18, above 0.05. A unit score is therefore not strong evidence in that case.
- Give Monte Carlo uncertainty. At 400 independent evaluable null replicates and a true rejection rate of 0.05, the approximate standard error is 0.0109. A few percentage points of difference do not establish a reliable ordering. Use the number actually defined for each metric when constructing intervals.
- Prespecified synthetic scenarios establish behavior under their generators, not biological validity, read-depth robustness, patient exchangeability in the real cohort, genome-wide error control, or novelty.

## Material qualifications for `report.md` section 5

1. **Permutation inference:** the current statement that every aggregate p-value equals 0.005 and therefore neighbours “really do” line up is too strong. The observed means exceed shuffled means descriptively, but shared patient composition/batch can produce this result. In addition, the two matrices used by the common shuffle contain overlapping CpGs across edges. They are not disjoint independent blocks under a simple all-edge null, so the aggregate permutation statistic has no automatically justified genome-wide p-value. Per-pair tests still require sample exchangeability. A shuffle diagnostic does not show excess proximity-specific coupling; that needs an appropriate nonneighbor comparison.
2. **Distance:** the within-R-state comparison currently cites Task 39, whose complete-site pairs can bridge missing sites. Label this as the older pair population, or recompute it on Task 41 true input-row adjacency. Describe the displayed medians as a descriptive distance trend; the large-gap bins contain few pairs. “All scores fall” should not imply a universal monotone law or include the lone 205-bp pair.
3. **Xi:** replace “catches any pattern” with its intended general-dependence interpretation and note n = 18 limitations. Xi need not describe sites going up and down together; inverse and nonlinear dependence can also give positive xi. When one level vector is constant, the *symmetric mean* xi is undefined; constant predictor/variable response has a defined directional tie-average of zero.
4. **Score definitions:** level similarity uses the mean **absolute** level difference. Exact agreement equals one minus **normalized** Hamming distance. Quadratic-weighted kappa equals Lin concordance applied to the same discrete-level vectors, not generally the raw-beta concordance score.
5. **Denominators:** change “median, all 32,664 pairs” to “median among defined pairs,” with the coverage column retained. “Works” should mean either numerically defined or demonstrably useful, not conflate the two.
6. **Candidate interpretation:** the numerical equivalence of C* to a marginal rescaling of weighted kappa is strong grounds against a new-correlation claim. Its identical p-values apply to the specified within-pair fixed-margin permutation tests; cross-pair ranking can change. Do not broaden “adds nothing to testing” beyond that scope.
7. **Recommendation:** similarity S measures closeness, while signed Spearman measures monotone association and xi measures more general dependence. Their joint use is a useful definition of the task output. It does not imply that xi and Spearman are interchangeable or that all observed association is local biological coupling.
8. **Provenance:** identify adjacency as consecutive rows in the supplied CpG file unless reference-genome completeness has been independently checked. The analysis preserves coordinate origin but does not verify whether it is zero- or one-based.

## Implementation review after the script was built

Read-only inspection of `Scripts/Task43.Item5.ControlledBenchmark.Sep29.2026.R` found no correctness blocker. The implementation extends the proposed design to nine fixed scenarios, including a bin-boundary null, an offset signal, and a signal entirely inside one level.

The cached distance matrices use correct double centering or U-centering. In R's column-major storage, index `p_i + n*(p_j-1)` correctly applies the same patient permutation to both axes of a distance matrix. Each synthetic replicate receives its own independent patient permutation; all methods share that permutation within the replicate. Centering commutes with this permutation, so caching the kernels is valid. The biased distance-correlation normalization is correct, and the unbiased squared distance-covariance divisor is `n*(n-3) = 18*15` for this explicitly fixed-n benchmark.

The script checks finite-statistic masks under every permutation, retains negative unbiased statistics, gives zero U-statistics a defined no-call test rather than dropping them, uses inclusive plus-one permutation p-values, and separates rejection/defined from rejection/all. Its binomial intervals use independent synthetic replicates. The shared-factor scenario is correctly labelled marginally dependent, not an independence null. These implementation findings concern the inspected code; completion counts and results must come from the completed run artifacts.
