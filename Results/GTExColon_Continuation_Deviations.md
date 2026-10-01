# Continuation audit: frozen sections 5 and 6

The Task 47 score and protocol SHA-256 hashes match their recorded values.
The GTEx colon build completed before this continuation: 1,839,097 autosomal
CpGs have at least 10 called reads in each of the 29 donors. The existing
Task 50 cache and count reconstruction are retained. No score coefficient,
cohort, coverage threshold, fold seed or success threshold has been changed.

## Corrections made before running Task 51

1. **Constant predictors in the held-out straight-line model.** The untested
   script inferred constants from computed variance being zero. A repeated
   value of 1/3 across 26 training samples produced a tiny positive variance
   and slope -1.538462 against an unrelated varying response. A mathematically
   constant predictor must give the intercept-only prediction. The corrected
   code detects exact constants by their range. Both legacy and corrected
   held-out gains are retained in fold tables and pair outcomes, with separate
   hypothesis summaries. This corrects the prediction endpoint implementation;
   the frozen score file is untouched.
2. **Reverse-chain coordinate mapping.** Annotation now lifts the two-base CG
   dyad, requiring a single two-base target interval. The leftmost base is the
   target CpG C even for reverse-oriented mappings. Lifting only the source C
   can instead yield the target G. Pairs with either unannotatable endpoint
   remain unknown for H6; they are not silently classified as unflagged.
3. **Annotation-build discrepancy in the inherited plan.** The cohort-choice
   addendum mentions hg19 annotation tracks, whereas the already written
   Task 51 script and downloaded tracks use hg38 via liftOver. This continuation
   retains those track definitions and makes the coordinate conversion explicit.
   This is an annotation adaptation/deviation, not an exact implementation of
   the addendum's statement. It affects H6 only; the pair gaps and all scores
   still use the original hg19 coordinates.
4. **Computational changes.** Unused large cache components are released;
   annotation imports are restricted to query loci; top-decile boundaries are
   shared across endpoints without changing tie weights. Fold-level tables are
   written as progress checkpoints. These changes do not tune the method.

## Input audit

`GTExColon.Continuation.Checks.R` verifies the frozen hashes, count-matrix
dimensions, unique donors, methylated <= total counts and the >=10 threshold,
and every retained site's genome-index position. On the first donor and the
donor with the lowest original beta-file agreement, the correct map was
compared with CpG-index offsets -1, +1 and +100. For the retained cohort sites,
correct-map mean absolute beta differences are approximately 0.00309 and
0.00266, versus approximately 0.059–0.067 for one-CpG shifts. Among intermediate
GEO betas, exact agreement is approximately 64.3% and 57.5%, versus 4.5–6.9%
for one-CpG shifts. Thus the map is strongly supported, but the PAT-derived
counts are not an exact reproduction of the separate GEO beta product.
The raw comparisons are in `GTExColon_MappingAudit.csv`.

The official PAT specification defines one character per consecutive genomic
CpG, a genome-wide CpG index and a multiplicity count; chromosome order follows
that index rather than UCSC lexical sorting. Source:
[wgbstools PAT format](https://github.com/nloyfer/wgbs_tools/blob/master/docs/pat_format.md).

## Section 6 interpretation to retain

The inherited build stores linkage for **adjacent genome CpGs**, matching the
user's handoff description. The protocol's phrase "each CpG pair within 200 bp"
is broader if interpreted literally; the reported H8/H9 scope must explicitly
say adjacent pairs. It does not cover every nonadjacent pair spanned by a read.

H9's inherited calculation relates cross-validated prediction outcomes to
linkage averaged over all donors and a low-spread label computed from all
donors. That is an association with held-out model outcomes, not fully
cross-fitted selection of pairs using other donors' linkage. The distinction
must remain explicit. It cannot establish causality or exclude correlated
sampling noise from reads shared by two sites.

The frozen section 6 success rules use point estimates. The inherited script
adds block-bootstrap intervals as descriptive uncertainty. Any computational
optimization of those intervals must be checked against expansion of sampled
blocks on toy inputs and documented, without changing the frozen verdict rule.

## Completed-run findings and additional implementation correction

The constant-predictor fix changed **zero** of the 1,574,581 pair-level gains
on this cohort. Legacy and corrected H1–H7 summaries are identical.

The first Task 51 run's multi-chromosome BigBed query returned annotations for
chr1 only, incorrectly marking all other chromosomes as nonunique. The final
script queries one chromosome at a time. This correction is recorded separately
in [Task51_H6_AnnotationFix.md](Task51_H6_AnnotationFix.md), including its timing
after the initial result. The final count is 146,998 flagged pairs; all three
H6 flags now satisfy the frozen rule in both halves and overall. The erroneous
initial H6 summary did not pass overall. This affects H6 only. The first
reliability table was not retained as an original file. The user's September 30
handoff supplies its three affected rows; they are transcribed and labeled as
such in the annotation-fix note. The final result and failure mechanism are
recorded alongside those historical values.

The exact block-bootstrap implementation was checked against explicit block
expansion for 150 toy resamples with unequal block sizes, ties and recomputed
medians. It uses multiplicity weights to reproduce the same statistic; Spearman
ranks are recalculated under each resample's multiplicities. Task 52 completed
with 1,000 resamples. These genomic-block intervals do not resample donors and
do not represent population-level uncertainty across new patients.

No frozen file was modified. Full results and remaining scope limitations are
in [GTExColon_Validation_Report.md](GTExColon_Validation_Report.md).
