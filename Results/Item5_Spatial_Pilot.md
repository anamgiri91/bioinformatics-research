# Item 5: a concrete Manhattan–dependence–genomic-distance pilot

Date: 2026-09-29. This continues the requested candidate investigation after
reviewing the [primary literature](Item5_Manhattan_Spatial_PriorArt.md).

## Candidate and decision

We implemented the exploratory affinity

\[
S_{ij}=1-\frac{1}{n}\sum_k|\beta_{ik}-\beta_{jk}|,\qquad
A_{ij}(\ell)=e^{-|p_i-p_j|/\ell}\frac{S_{ij}+\operatorname{dCor}(\beta_i,\beta_j)}{2}.
\]

The primary length scale is 100 bp; 50 and 200 bp are sensitivity settings.
These choices were recorded before computing this pilot's outcomes and were
not fitted. The data had already been examined in earlier work, so this is
internal exploratory validation, not a preregistration or independent cohort.
The quality multiplier is one: the ratio-only input supplies no read depths
or conversion-quality estimates. A constant profile leaves dCor and the
composite undefined; its Manhattan agreement remains available separately.

**Decision: do not promote this equal-weight combination as an improvement
for predicting absolute agreement.** On the stated held-patient endpoint,
Manhattan similarity alone performed better. Spatial weighting improved the
mean result of the unweighted blend slightly, but did not recover the
performance of the simpler agreement score. This does not test every possible
combination, dependence target, region definition, or biological endpoint.

The ingredients have published precedents. This experiment is a concrete
adaptation to our RRBS matrix, not evidence of a new statistical principle.

## Data and endpoint

The input is `Data/NN.hg38.18P.forw.chr22.w.header.txt`. We used all 32,664 pairs
that are consecutive in the original input and complete in all 18 patients.
Incomplete intervening sites were not bridged. The gap range is 2–205 bp,
with median 8 bp; these pairs cannot establish kilobase-scale behavior.

Each patient was held out once. Training on the other 17 patients, each method
selected the highest-scoring 10% of eligible pairs. The outcome is their mean
absolute beta difference in the held-out patient, so lower is better. This
measures agreement transfer to another patient. It does not measure discovery
of nonlinear dependence, direct co-regulation, or technical accuracy.

All methods used the same within-fold eligible population: both training
profiles must vary, leaving 30,095–30,694 pairs per fold. Eligibility and scores
do not use the held-out patient's beta values. The complete-case universe
does depend on observed availability across all patients, and limits
generalization to sparsely observed sites. Genomic position breaks score ties
deterministically; this especially affects the gap-only baseline and can
concentrate tied selections in earlier regions of the chromosome.

## Results

Values below are the equally weighted mean of the 18 patient-fold means.

| Ranking score | Held-out mean absolute beta difference | Median selected gap across folds |
|---|---:|---:|
| Manhattan similarity alone | 0.01635 | 5 bp |
| Distance correlation alone | 0.05518 | 6 bp |
| Equal-weight mixture, no genomic weighting | 0.04272 | 6 bp |
| Mixture × genomic weight, 50 bp scale | 0.04032 | 3 bp |
| Mixture × genomic weight, **100 bp scale** | **0.03990** | **4 bp** |
| Mixture × genomic weight, 200 bp scale | 0.04007 | 4 bp |
| Genomic proximity alone | 0.04056 | 2 bp |

Patient folds share training data, and neighboring pairs share CpGs. These
rows are descriptive comparisons, without an independence-based standard error
or significance claim. The endpoint directly concerns absolute agreement,
whereas dCor also rewards inverse/shifted/nonlinear relationships; its worse
performance here is not evidence that dCor is a worse general dependence test.

With all 18 patients, raw dCor and symmetric raw-beta xi are defined for
30,724 pairs; symmetric ten-level xi is defined for 24,726. Raw xi was added
to the saved pair scores as an important missing comparator from the earlier
simulation study; this pilot does not rerun that simulation with raw xi.

The 100-bp weight retains 73.97% of the unweighted blend's top 10% of pairs,
with Spearman rank correlation 0.86475 across all 30,724 defined pairs. Thus
spatial weighting materially changes prioritization even in this short-gap
population. Changed ranking alone is not improved validity.

## Mathematical check and validation

For a fixed pair and any fixed positive spatial weight w,
`w*T_permuted >= w*T_observed` exactly when
`T_permuted >= T_observed`. Its within-pair permutation p-value is unchanged.
We checked the complete 18 possible alignments of two profiles each containing
one high and 17 low values. The aligned profiles have exact p=1/18 at every
tested genomic gap, despite their maximal unweighted agreement/dependence.

Implementation checks passed for the distance-correlation formula against an
independently expanded expression, symmetry, identical profiles, constant
handling, original-row adjacency, score bounds and spatial permutation
invariance. No read counts, biological labels or quality measurements were
inferred from beta ratios.

## Implication for further method design

The question must determine how factors are combined. Absolute agreement,
patient-pattern dependence and spatial region membership are different
targets. For an agreement score, the present result favors keeping Manhattan
similarity as the base and testing whether context explains its residual
errors. For a dependence or regional score, use an endpoint that actually
rewards those properties, and include raw xi/dCor and simpler spatial methods.

A context-relative score could compare observed similarity with the expected
similarity at the same gap, baseline methylation, dispersion and local density.
This is a proposed next experiment, not an implemented or proven novel method.
Any fitted context model needs validation on separate genomic blocks; its
parameters must not be selected by inspecting the validation outcomes.
Build-matched reference CpG density and tissue-specific annotations are needed
before claiming to test those biological factors.

## Reproducible files

- [Runnable R script](../Scripts/Item5.SpatialAffinity.Pilot.R)
- [Protocol](Item5_Spatial_Protocol.json)
- [All pair scores, compressed CSV](Item5_Spatial_PairScores.csv.gz)
- [Held-patient results](Item5_Spatial_HeldPatientFolds.csv)
- [Method summary](Item5_Spatial_HeldPatientSummary.csv)
- [Rank sensitivity](Item5_Spatial_RankSensitivity.csv)
- [Exact permutation check](Item5_Spatial_PermutationInvariance.csv)
- [Run log](Item5_Spatial_Run.log) and [R session](Item5_Spatial_SessionInfo.txt)
