# GTEx colon: completed continuation of the frozen validation

Completed 2026-09-30. The frozen F+ score and Task 47 protocol retain their
recorded SHA-256 hashes. This report covers the public-cohort validation and
the inherited **adjacent-pair** shared-read analysis, with deviations stated
below. It does not claim that every pair within 200 bp was tested.

The strongest result is a repeated failure of genomic distance weighting to
improve the planned prediction-versus-agreement comparison. Read-depth
weighting passes its planned rule, with a small absolute improvement. Shared
read linkage is associated with prediction gain, but that association does
not establish biological covariance between people or methodological novelty.

## Cohort and endpoints

- GSE233417, 29 distinct GTEx colon donors; hg19 coordinates.
- 1,839,097 autosomal CpGs with at least 10 called reads in every donor.
- 1,574,581 truly consecutive genomic CpG pairs; median gap 11 bp, range 2–531 bp.
- Ten patient folds, seed 20260929. Scores use training patients; endpoints use
  held-out patients. Top-decile selection splits ties fractionally.
- Endpoints: absolute beta disagreement and bidirectional absolute-error gain
  from a clipped linear prediction relative to the training-mean prediction.
- The halves are position halves within each chromosome. Passing requires the
  full set and both halves to meet the recorded rule.

Fold-by-error-level cells share data. Their win counts are the frozen decision
rule, **not independent trials or a binomial significance test**.

## Public-cohort results

| Hypothesis | Result | Evidence |
|---|---|---|
| H1: distance weighting helps F | Fails | F wins 3/60 cells in the full set and 3/60 in each half. |
| H2: F beats the planned New C | Passes | New C wins only 7/60 in every subset. |
| H3: dCor improves on shrunk Spearman inside New C | Passes | Shrunk Spearman wins 0/58 in every subset. |
| H4: sign guard has negligible cost | Passes | Largest reported fold-mean endpoint difference is 0.00001; threshold is 0.0005. |
| H5: spread guard costs gain | Passes | Guarded F wins 0/60 in every subset. |
| H6: flagged pairs behave worse at matched spread | Passes after annotation correction | Nonunique, poor-mapping and SINE flags each satisfy the rule overall and in both halves. |
| H7: weighting S by read depth helps | Passes | Weighted F+ wins 45/60 overall, 44/60 first-half and 47/60 second-half cells. |

H7's fold-mean gain advantage at the six matched agreement levels is roughly
0.00001–0.00008 beta units overall, from the five-decimal summary. It meets the
protocol but should not be described as a large practical improvement. H3
compares the components **inside New C**; it does not show universal superiority
of distance correlation over Spearman or Chatterjee's xi.

For H6, there are 146,998 nonunique, 358,242 poor-mapping and 284,514 SINE pairs;
414 pairs have unassignable flags. The flagged-minus-unflagged agreement-error
differences, matched on 20 spread bins, are respectively +0.00563, +0.00692 and
+0.00712; prediction-gain differences are −0.00421, −0.00472 and −0.00519.

[Hypothesis table](Task51_Hypotheses.csv),
[matched comparisons](Task51_Matched.csv),
[reliability details](Task51_Reliability_SpreadMatched.csv),
[frontier figure](Fig37_GTExColon_Frontier.png).

## Shared-read results

Linkage is computed within a donor from fragments calling both adjacent CpGs:

\[
L=P_{\rm observed}(\text{same state})-
  \{\hat p_i\hat p_j+(1-\hat p_i)(1-\hat p_j)\}
 =2(\hat p_{11}-\hat p_i\hat p_j).
\]

Each donor needs at least 10 jointly called fragments. The cache then averages
linkage equally over qualifying donors for each pair; the donor set can vary
between pairs. The analysis contains 4,323,033 adjacent pairs within 200 bp.

**H8 fails its implemented monotone-distance test.** Mean linkage is positive:
0.0362, with a 1-Mb-block bootstrap interval [0.0355, 0.0370]. But the across-pair
Spearman correlation with gap is **+0.0265**, interval [0.0229, 0.0301], rather
than negative. Both halves show the same signs. The mean by distance bin peaks
at 10–20 bp (0.0396) and declines to 0.0273 at 150–200 bp; this binned mean
pattern does not override the failed pair-level rank criterion. The protocol
says "falls with distance" without naming a statistic; Spearman was the
operationalization in the inherited Task 52 script, not a statistic explicitly
named in the frozen Markdown. That distinction limits a strict preregistration
claim for this part.

**H9 passes the recorded direction rule.** Among 911,729 low-spread pairs with
linkage, the above-median group has mean held-out prediction gain 0.00185,
versus −0.00013 in the at-or-below-median group. The difference is **0.00198**,
block interval [0.00192, 0.00204]. Differences in the halves are 0.00197 and
0.00198, with positive intervals. "Below median" in the inherited figure
includes ties at the median. Exploratorily, linkage and held-out gain have
Spearman correlation 0.5528 across 1,560,416 eligible cohort pairs.

[H8 details](Task52_H8_Detail.csv), [H9 details](Task52_H9_Detail.csv),
[distance bins](Task52_H8_ByDistance.csv),
[shared-read figure](Fig38_SharedReads.png).

## Interpretation and limits

1. **H9 is an association with cross-validated model outcomes.** Linkage and the
   low-spread grouping use all donors. Pair grouping is therefore not fully
   cross-fitted. This is not a prospective test selecting pairs solely from
   other people's read linkage.
2. **Shared fragments couple measurement errors.** Within-fragment concordance
   can produce covariance of estimated betas even when latent methylation
   proportions do not covary across people. H9 does not separate that from
   biological between-person covariance. Predicting observed beta and
   estimating biological covariance are different targets.
3. **Uncertainty is conditional on this cohort and assay.** Genomic block
   resampling does not resample the 29 donors. Millions of pairs do not replace
   independent people or independent cohorts. RRBS and complete-coverage
   selection further restrict generalization.
4. **Adjacent-pair scope.** The inherited build covers consecutive genomic CpGs.
   The frozen text "each CpG pair within 200 bp" is broader if interpreted
   literally. Nonadjacent shared-read pairs remain outside this analysis.
5. **Mapping and annotation.** PAT counts map much better to the correct genome
   CpG indices than shifted indices, but do not exactly reproduce the separate
   GEO beta product. hg19 pairs were annotated through two-base CG-dyad liftOver
   to hg38 tracks, a deviation from the data-choice addendum's hg19-track plan.

## What this changes for item 5

The TCGA validation already found H1 failures (0/57 and 5/52 wins) and H2–H5
passes. GTEx now repeats those directions with sequencing data. The earlier
outlier-support primary test remains a failure under its frozen criterion:
effect +0.1401 with patient-bootstrap interval [−0.0077, 0.1906]. It has not
been relabeled a success by the newer cohort results.

For further development, a transparent candidate is

\[
S_{ij}^{\rm depth}=1-
 \frac{\sum_k \min(N_{ik},N_{jk})|\beta_{ik}-\beta_{jk}|}
      {\sum_k \min(N_{ik},N_{jk})},\qquad
F_{ij}^{\rm next}=S_{ij}^{\rm depth}\{\alpha+(1-\alpha)E_{ij}^{+}\}.
\]

Here alpha is one common agreement/dependence tradeoff rather than a function
of gap. **This joint combination has not been validated:** the frozen study
tested removal of distance and addition of read-depth weighting separately.
The equation uses established operations, so it is a candidate adaptation,
not evidence of a new correlation coefficient. Alpha would need to be fixed
on development data and checked on a new cohort.

A more specific methodological target is distinguishing biological covariance
from covariance introduced by overlapping sequenced fragments. A separate
[derivation and simulation](Item5_SharedReadNoise_Design.md) demonstrates both
feasibility and important finite-sample failures. That work is exploratory and
does not modify the frozen F+. Any proposed advance now needs comparison with
measurement-error methods and validation against independent fragments or
technical replicates, not only prediction of the same noisy beta measurements.

## Corrections, checks and reproduction

- The exact-constant linear-predictor correction changes no pair-level gain
  here; legacy and corrected hypothesis tables are identical.
- A multi-chromosome BigBed import initially returned only chr1 annotations.
  The final script reads each chromosome separately. This changed H6's verdict;
  the initial failure and fix are documented in
  [Task51_H6_AnnotationFix.md](Task51_H6_AnnotationFix.md).
- The bootstrap optimization passed 150 comparisons with explicit resampled
  block expansion, including ties and median changes.
- Final checks verify frozen hashes, all 3,780 fold records, aggregation,
  legacy equivalence, annotation counts and H8/H9 table consistency.

Commands, from the repository root (Task 50's completed cache is required):

```sh
Rscript Scripts/GTExColon.Continuation.Checks.R
Rscript Scripts/Task51.GTExColon.FrozenTests.Sep29.2026.R
Rscript Scripts/GTExColon.LegacyEndpointSummary.R
Rscript Scripts/Task52.GTExColon.SharedReads.Sep29.2026.R
Rscript Scripts/Item5.SharedReadNoise.Pilot.R
Rscript Scripts/GTExColon.Validation.Checks.R
```

Task 52's exact rank bootstraps can be expensive on millions of pairs; the
completed tables need not be recomputed for reporting. Source and output
hashes are recorded in `GTExColon_Validation_Manifest.json`. See also the
[full deviation log](GTExColon_Continuation_Deviations.md).
