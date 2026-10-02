# Task 41: what the earlier CpG correlation work contributes

Date: 2026-09-28. Scope: contextual source review for supervisor item 5 in
[Task 39](../Scripts/Task39.RRBS.StatesAndNeighbours.Sep28.2026.R): discretize
RRBS beta values into ten levels and compare consecutive CpGs, including
Chatterjee's xi. This note does not rerun or independently validate the historical
experiments. It does not resume the separate outlier Branch B or its paused team.

## Feasibility and the correct question

The requested measurement is feasible on the 18-column chr22 matrix.
Task 39 has already produced 38,563 consecutive-complete-site pair scores in
approximately 13 seconds for the whole script, according to its saved
[transcript](../Scripts/Task39.RRBS.StatesAndNeighbours.Sep28.2026.Rout).
This establishes computational feasibility, not statistical novelty.

Two questions must remain visible: **Do paired samples have similar methylation
levels?** and **Do levels co-vary across samples?** Equal, constant vectors provide
perfect agreement but no variation from which to estimate correlation. A vector
and a shifted copy can have perfect Pearson correlation while disagreeing in
absolute methylation. A reversed vector can have strong dependence and poor
agreement. No single ordering of these measures makes all three cases equivalent.

The most transferable earlier implementation is
[src/similarity.py](<../correlation for cg sites/src/similarity.py>):

\[
\operatorname{mean}[(x-y)^2]=(\bar x-\bar y)^2+(s_x-s_y)^2
 +2s_xs_y(1-r),
\]

where the SDs use divisor n. The final
term can be written `2*(sx*sy-cov(x,y))`, which stays defined when a vector is
constant. This identity explains how Euclidean discrepancy combines differences
in level, spread and paired pattern. The source uses level and spread beside
dependence rather than treating them as evidence of an association by themselves.

For this task a useful primary agreement record is exact-bin agreement,
`1 - mean(abs(level_x-level_y))/9`, and raw-beta MAE/RMSE. Report raw concordance,
signed Pearson/Spearman and the two directional xi values alongside it. The
normalized ordinal L1 is an interpretable rescaling of an established distance,
not a new statistical coefficient. Raw-beta distances are needed because ten-bin
agreement loses within-bin information and changes abruptly at a boundary.

## Transferable lessons actually tested in the older study

| Saved finding | Consequence for the RRBS analysis |
| --- | --- |
| Robust weighting discarded rare switched samples and sharply reduced power in the recorded simulations. | Do not automatically delete a CpG pair's unusual patient. Report leave-one-patient-out influence, and retain rare patterns as possible signal with uncertainty. |
| Shared purity/batch shifts created dependence; no tested historical candidate passed the entire contract. | A large correlation or xi is an observed association, not proof of direct biological coupling or novelty. Do not transplant a tumor-purity adjustment to these samples without relevant metadata. |
| Shared median imputation caused serious false-positive inflation in the RS simulations. | Do not replace missing RRBS values with a common level. Score complete sites as a transparent primary set and use genuinely adjacent original rows with pairwise-complete samples as a separate missingness sensitivity analysis. |
| Rare-switch analytic Pearson/dCor tails were badly miscalibrated in tested settings. | A descriptive score is feasible without claiming a calibrated p-value. Use n=18, tie-preserving sample-permutation diagnostics if inference is examined; historical large-n calibration does not transfer. |
| Some screens discarded every pair in a difficult case. | Report denominators and undefined fractions per measure. No score is not evidence of no association or a passing result. |

Sources: [September 16 benchmark](<../correlation for cg sites/results/benchmark/2026-09-16/FINDINGS.md>),
[analytic-tail checks](<../correlation for cg sites/results/benchmark/tail_calibration/README.md>),
[sharp-null checks](<../correlation for cg sites/results/benchmark/tail_calibration_sharp/README.md>),
and [formula-search log, E4–E5](<../correlation for cg sites/docs/formula_search_log.md>).
Their numerical results concern their stated simulations at n=155/466/776 and
must not be quoted as measured n=18 RRBS performance.

## Small-sample and tie limitations

Eighteen observations assigned to ten levels necessarily create ties. One patient
changes the exact agreement fraction by 1/18, approximately 5.56 percentage
points. A nominal 2% rare subtype has expected count 0.36 in 18 observations;
the earlier 2% switch experiments therefore provide no rare-subtype validation
for this matrix. A 10-by-10 level table has at most 18 occupied cells, making
unqualified large-sample contingency-table approximations implausible here.

Chatterjee's paper explicitly permits discrete observations, gives the tie-aware
formula, and suggests averaging over all permissible predictor tie orders
(Remark 8). Xi is directional; maximum-direction symmetrization is also discussed
(Remark 1). For untied responses its finite-sample maximum is `(n-2)/(n+1)`:
16/19 at n=18, so xi is not a unit-scaled agreement fraction. Its constant-response
denominator is zero. A constant predictor with varying response is a different
case. Source: [Chatterjee, original paper, pp. 2–4](https://arxiv.org/pdf/1909.10140).

Exact averaging removes arbitrary random tie-order noise, while leaving sampling
uncertainty and information loss from discretization. For a predictor tie block
of size m with response ranks r, the expected within-block rank-jump total is
`sum_{u,v in block}|r_u-r_v|/m`. At a boundary between successive tie blocks it
is their mean cross-block rank difference. Summing these terms gives an exact
conditional average of the usual xi numerator. Constant predictor/nonconstant
response then gives zero. Distinguish `max(E[xi_xy], E[xi_yx])` from
`E[max(xi_xy,xi_yx)]`; these generally differ. This is an implementation of a
published suggestion, not a novel coefficient.

Keep negative finite-sample xi estimates rather than silently clipping them.
Report valid n, unique levels, constant flags and leave-one-out influence.
Bin-resolution sensitivity should use fixed prespecified alternatives and the
raw-beta measures. A sample-permutation diagnostic assesses dependence under an
exchangeability assumption; it cannot identify or correct shared batch effects.
Adjacent pairs share CpGs, so pair counts are not counts of independent
replicates. Genome-wide error control is a further analysis, not an automatic
consequence of running permutations on a subset.

## Task 39 details to correct or qualify

- The score table describes consecutive **complete** CpGs, which can skip
  intervening incomplete sites. Preserve genomic gap and adjacency in the
  original file; do not describe all such pairs as original-row neighbors.
- Its header says correlations, kappa and xi are undefined if either vector is
  constant. The code is more nuanced. Quadratic weighted kappa is undefined
  when both vectors are the same constant; with one constant and another
  varying vector it is zero. Raw Pearson can remain defined when a discretized
  vector is constant. Xi needs directional handling as described above.
- Task 39's `xi` is a maximum of two randomly tie-broken directional estimates.
  Seeded reproducibility does not quantify tie-order variation. A deterministic
  tie-averaged replacement must be named separately, preserving the old result.
- Its [same-state AUC table](Task39_RRBS_MeasureVsState_AUC.csv) reports different
  coverage: Pearson 94.17%, Spearman/xi 76.19%, and kappa 91.61% within 1 kb.
  Raw AUC rankings on different valid sets are not a fair universal ranking.
  States are derived from the same measurements and are not independent truth.
  This is a descriptive relation to state labels, not validation of a superior
  similarity measure. If compared, also report a common valid subset.
- No read-depth columns exist in this input. Ratio precision and technical
  reliability cannot be recovered from the ratios alone, and counts must not
  be inferred from rounded fractions. Neighbor differences also contain biology.

## What is and is not accomplished

A reproducible table with defined pair membership, complementary agreement and
dependence scores, tie/constant handling, missingness/bin/patient sensitivity and
a small n=18 diagnostic benchmark can accomplish **item 5's requested
measurement and comparison**. It cannot by itself prove an original method,
superiority for biological discovery, or a validated methylation mechanism.

The correlation folder's [AGENTS.md](<../correlation for cg sites/AGENTS.md>),
[roadmap](<../correlation for cg sites/docs/RESEARCH_ROADMAP.md>),
[plan](<../correlation for cg sites/docs/RESEARCH_PLAN.md>),
[novelty review](<../correlation for cg sites/docs/OUTLIER_NOVELTY.md>),
[CLAUDE.md](<../correlation for cg sites/CLAUDE.md>), and
[method guidelines](<../correlation for cg sites/skills.md>) concern the separate
per-patient outlier branch. Its pending choices, frozen gates and patient
restrictions remain intact. The [archived pairwise notes](<../correlation for cg sites/docs/LEGACY_CO_METHYLATION_NOTES.md>)
record that no historical pairwise candidate passed the full contract. Neither
branch establishes that the present task is novel, and this scoped analysis
does not reopen either branch's research decisions.

## Addendum: defining and testing a proposed alignment score

The user subsequently clarified that item 5 includes considering a defined
similarity score or a new correlation. A mathematically explicit candidate and
its counterexamples therefore belong in the deliverable. Its novelty and
usefulness must be assessed separately from successfully implementing it.

For paired ordinal levels x and y on n complete samples, define

\[
D=\frac1n\sum_i|x_i-y_i|,\quad
D_0=\frac1{n^2}\sum_{i,j}|x_i-y_j|,
\]

\[
D_{\min}=\frac1n\sum_i|x_{(i)}-y_{(i)}|,\quad
D_{\max}=\frac1n\sum_i|x_{(i)}-y_{(n+1-i)}|.
\]

Sorted matching minimizes, and reversed matching maximizes, the total absolute
distance over permutations. Consequently
`D = Dmin + (D-Dmin)` decomposes paired disagreement into empirical marginal
distribution mismatch and excess due to the particular patient pairing. Dmin
is the one-dimensional empirical Wasserstein-1 distance; it does not use the
original patient correspondence.

The corrected candidate is

\[
C^*=\frac{D_0-D}
{\max(D_0-D_{\min},\;D_{\max}-D_0)}.
\]

Use NA when its denominator is zero. This convention distinguishes absence of
alignment information from an observed score near the chance reference. C* is
symmetric, deterministic with ties, bounded in [-1,1], and has exactly zero mean
over all permutations of patient correspondence conditional on its empirical
margins. One endpoint need not reach magnitude one. These properties follow
because D0 is the permutation mean of D and the denominator is fixed throughout
that permutation distribution. They do not establish an unbiased population
estimator, confounder adjustment, a general independence characterization, or
small-sample significance.

This construction closely overlaps established agreement methodology. Linear
weighted kappa is `kL = 1-D/D0`; therefore `C* = kL*D0/M`, where M is the stated
denominator. Fixed-margin maximum normalization is already studied by
[Rapallo, 2022](https://link.springer.com/article/10.1007/s11336-022-09844-y).
Piecewise min/chance/max scaling of nominal agreement appears in
[Safak's 2020/2021 preprint](https://arxiv.org/abs/2006.12904).
C* should be described as a proposed normalization/adaptation of weighted
agreement. No claim of a new statistical principle follows.

Three decisive tests constrain the interpretation:

1. **Sparse support:** both profiles contain seventeen 1s and one 2. Shared
   spikes give C*=1, but their exact one-sided permutation probability is 1/18.
   Nonoverlapping spikes give C*=-1/17. Removing the shared spike makes both
   profiles constant and C* undefined. A score of one is therefore not strong
   evidence with 18 patients. Report the denominator and leave-one-out behavior.
2. **Separated ranges:** let x repeat levels 1,2,3 and y=x+7 repeat 8,9,10. All
   pairings have the same D, so Dmin=D0=Dmax and C* is undefined, although
   Pearson correlation is one. This is a limitation of absolute-distance
   agreement for detecting shifted patterns, consistent with the earlier SAED
   failure. Constants are another denominator-zero case.
3. **Rejected piecewise denominator:** using D0-Dmin for positive scores and
   Dmax-D0 for negative scores maps the sparse example to +1 with probability
   1/18 and -1 with probability 17/18. Its permutation mean is -8/9, despite
   being zero at D=D0. The single-denominator correction avoids this bias.

For fixed margins C* is monotone in -D and yields the same one-sided permutation
p-value as paired L1 or linear weighted kappa. It adds interpretation and changes
comparisons between different margins; it adds no new within-pair ordering
information. Bounded absolute differences limit each patient's raw contribution,
but a small denominator can amplify influence in C*. Ties, n=18 and absent read
counts do not justify claiming robustness to technical error. Keep raw L1/MAE,
Dmin, excess mismatch, opportunity M, standard comparators and directional xi
beside the candidate. Implementing and testing this explicit score addresses
the clarified task while leaving scientific novelty unproved.
