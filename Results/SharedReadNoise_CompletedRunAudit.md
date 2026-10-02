# Completed-run audit: blood and plasma shared-read correction

Date: 2026-10-01. This is a post-result audit, not a new locked test.

## What was verified

- All 63 entries in the Task 58 SHA-256 manifest and all 54 entries in the
  Task 64 manifest match the current files. These 117 entries can include
  the same file in both manifests; they are not 117 distinct files.
- The saved Task 59 and Task 65 pair-level estimates reproduce the primary
  CSV point estimates, relative error changes, per-pair improvement shares,
  disagreement endpoints, and the deterministic bootstrap-subset point
  estimates at the precision saved in those CSVs.
- Covariance loss was independently recomputed as `(estimate-reference)^2`
  from each pair's saved covariance estimates, rather than trusting its
  saved loss column. The identity `corrected = raw - mean_noise` holds.
- The pair keys are unique, gaps are in (0,200] bp and all retained pairs
  have at least 20 eligible donors. The plasma short-gap endpoint and
  pair count reproduce the saved sensitivity analysis.
- The existing independent estimator checks pass: 300 count, PAT,
  streaming and exact-moment checks in `SharedReadNoise.Checks.R`.

This audit reads saved estimates. It does not independently reconstruct
counts from raw alignments, verify molecule identities from read names,
rerun read splitting or regenerate the 2,000 bootstrap draws. Hash matches
establish current integrity against a recorded lock, not public timestamping
or a guarantee that every assumption in the locked protocol holds.

## Results supported by the saved outputs

| Quantity | Blood RRBS, Task 59 | Plasma capture, Task 65 |
|---|---|---|
| Cohort donors | 57 | 46 |
| Eligible pairs | 4,408,849 | 995,729 |
| Mean squared-error difference | -9.343e-6 | -3.498e-5 |
| Relative error change against the reference | -47.86% | -48.33% |
| Saved 95% interval | [-1.099e-5, -8.219e-6] | [-4.167e-5, -2.949e-5] |
| Eligible pairs used for bootstrap | 176,146 | 39,792 |
| Recorded primary criterion | Met | Met |

The intervals use the locked CpG-index-modulo-25 subset, whereas the full
point estimates average all eligible pairs. The full and subset point
estimates are close, but that does not make the intervals full-pair
bootstraps. The intervals subtract the bootstrap mean shift from the
percentile endpoints; the protocol calls this "bias-corrected percentile."
It is not a BCa interval. Task 57's calibration used three simulation
scenarios, 300 replicates and 400 bootstrap draws per replicate, with the
target endpoint estimated from those same simulation replicates. Broader
coverage and model-violation calibration remain outstanding.

## The plasma reference is not proven fragment-independent

Task 64 section 7 explicitly says some mates remain separate mHap lines.
Those lines can be allocated to different partitions. Consequently:

1. Shared physical-fragment covariance may be absent from the reconstructed
   joint counts.
2. B and C can share physical-fragment noise despite using disjoint lines.
3. The errors in A and in the B/C reference need not be conditionally
   independent, so the usual unbiased risk-difference identity is not
   guaranteed for this input.

The at-most-40-bp sensitivity still meets its criterion, but it neither
measures the residual cross-partition covariance nor proves it is zero.
Report the measured error reduction and the declared limitation together.
Resolving this requires physical fragment identity from alignments or
another verified molecule-preserving data source; aggregated unlinked
records cannot recover which mates belonged together.

The mathematical split target is the latent covariance of the sampled
eligible donors. Even with an ideal fragment split, a risk comparison about
a population covariance requires additional reasoning and sampling
assumptions, as specified in the plan.

## Distance interpretation

The plasma (150,200] bp bin has 1,548 pairs. Every one has zero represented
overlap in A and zero estimated correction. This establishes what the
available records contain, not the absence of a physical fragment covering
both CpGs. Do not attribute that zero to a typical cfDNA fragment length.

From the (0,10] to (150,200] bins, the plasma mean observed covariance
declines 84.8%, the corrected mean 42.4%, and the B/C mean 47.3%. The first
bin has 525,074 pairs. These are descriptive between-bin differences;
depth, represented overlap, composition and eligibility vary. No causal
effect of fragment length or a controlled biological decay rate has been
identified. The new figure shows these descriptive profiles without
uncertainty intervals.

## Correlation, ranking and novelty

Plasma's correlation reference is defined for 19.64% of pairs. Among 23,098
pairs with all three estimators and the reference valid, Pearson's MSE is
0.221, versus 0.350 for variance-only and 0.399 for full correction. This
noisy ratio reference and validity selection limit latent-correlation
claims. There is no evidence of a practical correlation improvement in
these comparisons. Raw disagreement and Manhattan agreement remain the
strongest tested rankers for the agreement endpoint.

The existing [prior-art audit](Task60_PriorArt.md) treats the subtraction
and count estimator as applications of known identities. This output audit
does not independently establish novelty. A possible contribution is the
CpG-specific empirical measurement of the bias and its effect on standard
analyses; an improved item 5 similarity score has not been demonstrated.

## Reproduce the audit

From the repository root:

```sh
Rscript Scripts/SharedReadNoise.Checks.R
Rscript Scripts/SharedReadNoise.CompletedRunAudit.R
```

The second command writes only files with prefix
`Results/SharedReadNoise_CompletedRunAudit`: full-precision aggregate and
distance tables, a descriptive figure and a manifest of audited inputs.
The frozen analysis scripts, locks and result tables are unchanged.

Next: recover or independently verify fragment identity; prespecify a
standard-analysis endpoint; then select and lock an independent WGBS
evaluation with design-specific calibration. Do not retune on plasma and
call a repeat analysis independent validation.
