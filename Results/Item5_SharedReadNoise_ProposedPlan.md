# Proposed plan: shared-fragment correction of between-person CpG covariance

Date: 2026-10-01. Status: a concrete proposal for novelty review, not a frozen
protocol or a completed real-data validation. It extends the
[checked moment derivation](Item5_SharedReadNoise_Design.md) and should be read
with the [prior-art assessment](Item5_SharedReadNoise_PriorArt.md).

The proposal investigates whether using observed shared-fragment joint states
improves estimation of latent between-person CpG covariance and agreement.
It does not establish a new bounded correlation coefficient. The first
implementation will expose unstable estimates rather than silently regularize
them. A constrained latent model would be a separately specified extension.

Review clarification, 2026-10-01: the fragment-split risk comparison below
targets the latent sample covariance of the eligible donors, conditional on
their proportions and the sampling design. It does not automatically compare
MSE for a superpopulation covariance. See the
[draft audit and derivation](Item5_SharedReadNoise_DraftAudit.md).

## 1. Target and initial scope

For adjacent reference-genome CpGs i and j, estimate the covariance across
donors of their underlying methylation proportions p_ik and p_jk. Also estimate
their average squared disagreement. The target is bulk-sample co-methylation;
it does not identify direct regulation or eliminate cell-mixture effects.

Start with autosomal, truly consecutive reference CpGs separated by at most
200 bp. This is a declared starting scope, not a biological boundary.
Genomic distance is a diagnostic stratum, not a multiplier in the estimator.
Use continuous beta values; ten-level discretization is a comparison only.

## 2. Rebuild the necessary sufficient counts

For every donor k and pair, retain:

- N_ik and N_jk: numbers of fragments calling each site;
- M_ik and M_jk: their methylated counts, giving b_ik = M_ik/N_ik;
- r_k: number of fragments calling both sites;
- n_00,k, n_01,k, n_10,k, n_11,k: joint states on those r_k fragments;
- coordinates, donor ID, coverage, and available mapping/QC annotations.

The first digit indexes site i. Thus p11hat = n11/r, pihatR = (n10+n11)/r,
and pjhatR = (n01+n11)/r. These are shared-subset proportions. Do not substitute
all-read site betas into the joint covariance expression.

Count a physical fragment once at a site. Overlapping paired ends must be
resolved before counting, with an explicit duplicate policy. Verify what a PAT
multiplicity represents against its producing pipeline; a compressed pattern
is not itself a unique molecule. If fragment identity or provenance cannot be
verified, the three-way fragment validation below cannot be claimed from that
file. Current pooled linkage summaries alone are insufficient for the method.

The development cohort is RRBS. Do not apply generic coordinate deduplication
to it: [Bismark's guidance](https://felixkrueger.github.io/Bismark/usage/deduplication/)
specifically advises against that operation for RRBS and related libraries.
Choose molecule/UMI handling from the actual library design and record residual
dependence when distinct molecules cannot be resolved. Counts N and r refer to
fragments with usable calls, not merely fragments geometrically spanning sites.

Initial count eligibility: both N values >= 2; r = 0 or r >= 2. Mark r = 1 as
unsupported for the moment correction, not zero. Require at least 20 eligible
donors per pair and apply the identical donor set to every comparator. Report
eligibility fractions. Pairwise results concern those eligible donors; coverage
selection can limit generalization.

## 3. Exact estimator

Under conditionally independent fragments, representative shared/nonshared
subsets, and independent donors, estimate within-donor sampling terms by

\[
\widehat v_{ik}=\frac{b_{ik}(1-b_{ik})}{N_{ik}-1},
\qquad
\widehat c_k=\frac{r_k^2}{N_{ik}N_{jk}(r_k-1)}
 \left(\widehat p_{11,k}-\widehat p^R_{i,k}\widehat p^R_{j,k}\right).
\]

Set chat to zero for r = 0 under this model. The second expression applies only
for r >= 2. Let s_ij and s_i^2 be sample covariance and variance across the same
n eligible donors, using denominator n-1. Return

\[
\widehat C_{ij}=s_{ij}-\frac1n\sum_k\widehat c_k,
\qquad
\widehat V_i=s_i^2-\frac1n\sum_k\widehat v_{ik}.
\]

For agreement, return

\[
\widehat D^2_{ij}=\frac1n\sum_k
\left[(b_{ik}-b_{jk})^2-\widehat v_{ik}-\widehat v_{jk}
 +2\widehat c_k\right],
\qquad \widehat S_{ij}=1-\widehat D^2_{ij}.
\]

The population target S lies in [0,1]; its finite-sample moment estimate need
not. This corrects squared disagreement, not Manhattan absolute disagreement.
The covariance correction can be positive or negative according to shared
states; it is not a fixed penalty for proximity or for sharing reads.

The secondary ratio is rhohat = Chat/sqrt(Vihat Vjhat). Return an explicit
invalid-status code if either variance is nonpositive or the ratio exceeds
[-1,1]. Preserve the raw ingredients and diagnostic ratio. Do not replace
invalid values by zero, clip them, or hide their frequency. Even a valid ratio
is not automatically unbiased or reliable.

Use 2,000 whole-donor bootstrap draws for initial 95% percentile intervals for
C and D^2; assess their actual coverage in simulation before interpreting them
as calibrated. The first version does not claim a calibrated ratio interval.
Do not generate such an interval by discarding invalid bootstrap ratios.

## 4. Comparators and endpoints

Compare like quantities:

- Covariance estimation: raw sample covariance versus full correction.
- Correlation: raw Pearson, variance-only correction s_ij/sqrt(Vihat Vjhat),
  and full shared-fragment correction, with all failure rates reported.
- Agreement: raw mean squared disagreement versus corrected D^2.
- Pair ranking: Manhattan similarity, Spearman, distance correlation,
  Chatterjee's xi, and the existing frozen scores, using a common held-out
  endpoint. These are not estimators of the same covariance target.

Add an existing measurement-error estimator where its assumptions and inputs
can be matched. Specify that adaptation before external evaluation; do not
label an unimplemented comparator as tested. Variance-only correction leaves
the covariance numerator unchanged and is therefore not a distinct covariance
baseline.

## 5. Development and controlled tests

Use already examined GTEx colon data for development. Rebuild per-donor counts,
verify count identities and mapping, and quantify correction magnitude against
coverage, shared-fragment fraction, distance, and methylation spread. Keep
donor means equal-weighted in the primary target; depth weighting changes it.

Extend the existing simulations across donor count (20, 29, 50, 100), depth
(including low and unequal depth), zero/partial/full overlap, low/zero latent
variance, and positive/zero/negative latent covariance. Choose joint-state
probabilities within their feasible bounds. Include separate misspecification
experiments for dependent fragments, nonrepresentative overlap, calling errors,
and varying sample composition.

Record covariance and disagreement bias, mean squared error, interval coverage,
and correlation failure rates. Check the no-overlap covariance correction is
exactly zero. Quantify Monte Carlo uncertainty. Improved bias alone is not a
pass: the correction may increase mean squared error at low coverage.

Before external evaluation, freeze the simulator grid, software, sample/pair
eligibility, diagnostic bins, estimator, and endpoints in a versioned protocol.
Any failed calibration requires development changes and a new protocol version.

## 6. External validation without using the same fragments twice

Select an independent cohort from metadata before examining its pairwise
results. Proposed eligibility is at least 30 donors from one tissue, with raw
fragment data or independently prepared replicate libraries and adequate
coverage. The accession has not yet been selected. Its identity and processing
rules must be fixed before evaluation. If no qualifying dataset is available,
report the external-validation requirement as unmet.

For fragment-level validation, randomly assign each physical fragment within
each donor to A, B, or C with equal probability and a fixed seed. Keep all CpGs
and both read ends of that fragment in its assigned partition; do not split
individual CpG calls independently. Fit both raw and corrected estimators on A.
Construct the reference from B and C:

\[
C^{ref}_{ij}=\tfrac12\{s(b_i^B,b_j^C)+s(b_i^C,b_j^B)\}.
\]

Under the independent-fragment model this reference has no same-fragment
sampling covariance. Define T_ij = s(p_i, p_j), the latent sample covariance
of these same eligible donors. Conditional on their proportions and a
representative sampling design, E[Cref | p, design] = T_ij. It is a noisy
reference, not ground truth. For squared agreement use

\[
D^{2,ref}_{ij}=\frac1n\sum_k
(b_{ik}^{B}-b_{jk}^{B})(b_{ik}^{C}-b_{jk}^{C}).
\]

Use the same eligible donors for A and the reference. Require N >= 2 at both
sites in each partition and r_A = 0 or r_A >= 2. Retain the 20-donor minimum.
Report the loss of pairs from partitioning. Conditional independence of A from
B/C is essential: it lets paired squared-error comparisons cancel reference
noise in expectation for the conditional, finite-cohort risk comparison.
This cancellation does not
automatically hold for risk about a fixed population covariance: both A and
the reference contain the same randomly sampled donors. Independent library
replicates provide a stronger check against some library-specific artifacts.
Splitting a library cannot eliminate
all shared systematic errors.

Primary external endpoint:

\[
\Delta_{MSE}=\operatorname{mean}_{pairs}
[(\widehat C^A-C^{ref})^2-(s^A_{ij}-C^{ref})^2].
\]

Negative values favor correction. Obtain uncertainty by resampling donors and
genomic blocks, preserving the same resampling and pairs for both methods.
Use fixed 1-Mb blocks within chromosomes and 2,000 joint bootstrap replicates
as the initial specification; verify this procedure in development simulations.
Do not treat millions of adjacent pairs as independent observations.

The prespecified primary improvement criterion is an upper 95% confidence
limit for Delta_MSE below zero. Report coverage/overlap/spread strata regardless
of the result. Agreement improvement, pair rankings, and correlation stability
are secondary outcomes and cannot rescue a failed primary endpoint. Passing
this test establishes improved covariance estimation in this setting, not a
universally reliable correlation coefficient.

## 7. Exact novelty question to investigate

“Has a method used per-donor joint methylation counts on overlapping DNA
fragments to estimate and remove the heteroscedastic shared-sampling covariance
from across-donor CpG beta covariance, with corresponding variance and squared
agreement corrections and evaluation against independently sampled fragment
measurements?”

Equivalent likelihood methods or general measurement-error estimators
specialized to these counts may anticipate the proposal even if their
terminology differs. The benchmark design and a new application alone do not
establish a new estimator. The previously completed work is a derivation and
simulation pilot; the count rebuild, robust inference, and external evaluation
above remain proposed work.
