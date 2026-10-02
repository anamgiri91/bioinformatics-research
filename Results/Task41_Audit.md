# Task 39 audit for supervisor item 5

Audit date: 2026-09-28. Existing `Scripts/Task39.RRBS.StatesAndNeighbours.Sep28.2026.R` and `Results/Task39*` were inspected without modifying them. Counts below were independently recomputed from `Data/NN.hg38.18P.forw.chr22.w.header.txt`; the AUC comparison uses the existing rounded pair export.

**Item 5 is computationally feasible.** The file supports ten-level paired methylation similarity for 32,664 genuinely consecutive file CpG pairs observed in all 18 samples. Task 39 already implements most requested descriptive measures correctly. Its principal limitations are the definition of adjacency, unstable randomized xi at n = 18, and treating state agreement as validation of a similarity measure.

## Data and adjacency

- Input: 634,646 unique, sorted chr22 positions, 18 sample columns; observed values are finite and in [0,1]. No duplicate chromosome-position keys were found.
- 503,265 sites have all 18 values missing; 38,564 have none missing. Missingness is therefore a substantial selection issue, not a small cleanup step.
- Task 39 creates 38,563 pairs after discarding incomplete sites. **5,899 pairs skip at least one input CpG**; among those, median skipped sites = 39 and maximum = 37,611. Of these skipped pairs, 2,523 remain within 1 kb.
- `adjacent_in_file` is computed correctly in memory for this single chromosome but omitted from the pair export. The supplement should retain row indices and both adjacency definitions.
- Of the 32,664 truly adjacent complete pairs, 31,464 are at most 50 bp apart; 1,116 are 50–100 bp apart; 83 are 100–200 bp apart; one is 200–500 bp apart. **Every Task 39 complete-site pair above 500 bp skips input CpGs.** Its long-distance profile cannot be described as a profile of genuinely consecutive CpGs.
- A sensitivity analysis can retain original adjacency and use jointly observed samples. At minimum overlap 6, 10, 12, 15, and 18, there are respectively 70,273, 61,585, 57,592, 50,233, and 32,664 eligible pairs. These are descriptive sample-size tradeoffs, not evidence that a particular cutoff is optimal.

## Measures and constant vectors

The ten-level transformation implements the requested intervals, including beta = 1 in level 10. State classification matches the documented L/H/M/HM/LM/R rules. Exact agreement, ordinal mean absolute difference, beta RMS difference, rank correlation, and the quadratic-weighted kappa formula are internally consistent.

Two comments need correction: `L2beta` is RMS difference (Euclidean distance divided by sqrt(n)); and kappa is not undefined whenever either vector is constant. Of all 38,563 complete-site pairs:

| Discrete vectors | Pairs | Kappa undefined | Spearman/xi undefined |
|---|---:|---:|---:|
| Both variable | 29,929 | 0 | 0 |
| Exactly one constant | 5,672 | 0 | 5,672 |
| Both constant | 2,962 | 2,957 | 2,962 |

Kappa is zero with exactly one constant vector, or two different constants; it is undefined for identical constants because both observed and expected squared disagreement are zero. Discretisation makes 5,798 complete sites constant, versus 1,271 on raw beta. Constants must remain available for identity/distance descriptions even when correlation is uninformative. `undefined_corr` in Task 39 summaries measures only missing Pearson values; it does not describe undefined Spearman, xi, or kappa.

`agree1` means within one discrete level, not an exact raw-beta difference threshold of 0.1. Identity-by-state is an analogy to genotype agreement, not a genotype-derived estimator here. No coverage-weighted/binomial uncertainty analysis is possible from this ratio-only input because methylated and total read counts are absent.

## Chatterjee xi

Task 39 uses the correct tied-response denominator and uniformly randomized predictor ties. Taking the maximum of both directions is permitted in the original method, but it remains a different finite-sample statistic requiring its own calibration. The original paper explicitly discusses averaging over possible tie orders and warns that finite-sample xi need not reach one. These facts support a deterministic tie-averaged supplement; this stabilization alone is not a new correlation method. [Chatterjee, *A new coefficient of correlation*, sections 1–2](https://arxiv.org/pdf/1909.10140).

An independent diagnostic sampled 500 complete-site pairs with two nonconstant level vectors (R seed 20260929), then recalculated the current symmetric maximum 50 times per pair. The median within-pair SD caused solely by tie breaking was **0.1094**, and the median observed range was **0.4448**. After those repetitions, 50 response permutations per pair gave pooled mean directional xi of -0.00027 and 0.00082, but mean maximum **0.08780**. These are finite Monte Carlo diagnostics, not discovery tests or exact null expectations.

For tie groups g ordered by x, the exact expected rank-difference sum is

`sum_g (2 / |g|) sum_{i<j in g} |r_i-r_j| + sum_adjacent(g,h) mean_{i in g,j in h} |r_i-r_j|`.

Independent exhaustive enumeration of 12, 48, 720, and one allowed tie order in four small examples agreed with that formula to within 8.4e-17. A constant predictor and variable response yields tie-averaged xi = 0; a constant response is undefined. Both directional values should be exported and labelled; retain a flag indicating whether both margins vary. Averaging over ties before taking a directional maximum does not eliminate maximum-selection bias.

Xi measures dependence rather than numerical identity. For 18 balanced binary observations, both an identical pair and its reversed levels give xi = 0.8889. Even an identical untied 18-value vector gives only 0.8421. Thus `1 - xi` should not be presented as a methylation distance or compared directly against identity agreement as though they measured the same property.

## Validation and bounded next steps

The same-state AUC uses labels derived from the same 18 methylation values used to compute the scores. It is a descriptive association with a coarse rule, not independent biological validation. Different metrics also use different subsets because their undefined cases differ. Restricting every metric to the same 26,809 pairs within 1 kb changes the ranking context:

| Measure | Task 39 own-subset AUC | Common-subset AUC |
|---|---:|---:|
| Exact agreement | 0.5803 | 0.5078 |
| Negative level L1 | 0.6208 | 0.5420 |
| Negative beta RMS difference | 0.6535 | 0.5783 |
| Weighted kappa | 0.6839 | 0.7168 |
| Pearson | 0.6722 | 0.6900 |
| Spearman | 0.6572 | 0.6572 |
| Xi maximum | 0.6286 | 0.6286 |

A defensible completion of item 5 should use true file adjacency as the primary population; retain complete-site adjacency as a labelled sensitivity analysis; report exact agreement and `1 - mean(abs(level_x-level_y))/9` as direct similarity; retain beta RMS difference and correlations as complementary descriptions; and calibrate the exact implemented xi statistic by sample-label permutation with fixed margins. Partial-overlap outputs must retain their overlap counts. No independent-pair standard errors are justified because adjacent pairs share CpGs, and no independent biological replication or read-depth uncertainty is available in this file. Method superiority, discovery claims, and biological novelty require validation beyond state AUC and this chromosome.
