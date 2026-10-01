# Review of the shared-fragment covariance research draft

2026-10-01. This reviews the user-supplied document beginning “Deconfounding
Inter-Donor CpG Co-Methylation via Shared-Fragment Joint-State Latent Covariance
Estimation.” It is a scientific review, not evidence of implementation or a
novelty determination. The prospective plan has been clarified accordingly;
the frozen Task 47 analyses have not changed.

The core moment formulas agree with the checked derivation under its sampling
assumptions. The draft's conclusions about novelty, distribution-free validity,
preprocessing, and population-level validation are stronger than the evidence.

## 1. A directly relevant Manhattan precedent

Nazer et al.'s MPCI paper was published on 24 March 2026. It explicitly uses
weighted Manhattan similarities across sequencing reads and consecutive CpG
sites. This is direct prior art for the broad idea of a Manhattan-based CpG
similarity score. Its stated target is methylation-pattern consistency, not
our between-donor sampling-covariance correction. That distinction leaves the
narrower question open; it does not prove the correction novel.
[Primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC13035127/).

The draft should also retain the closer measurement-error references in the
[previous assessment](Item5_SharedReadNoise_PriorArt.md): MethylPCA,
Buonaccorsi et al., epiG, dSOMNiBUS, and MeasurementError.cor. Broad categories
such as “single-cell methylation mapping” are not individual comparators whose
absence of a feature can be asserted for an entire field.

## 2. Correlated measurement-error correction is established

The statement that classical measurement-error approaches universally assume
independent feature errors is incorrect. Saccenti et al. explicitly study
correlated errors and correction; MeasurementError.cor documents separate
estimates for error correlation and latent correlation.
[Saccenti et al.](https://www.nature.com/articles/s41598-019-57247-4),
[official package manual](https://bioconductor.org/packages/release/bioc/manuals/MeasurementError.cor/man/MeasurementError.cor.pdf).

Variance-only correction is one restricted baseline, not a characterization of
the whole measurement-error literature. A shared-fragment-specific estimator
must be compared with general correlated-error methods where applicable.

## 3. What the correction actually is

Let s^R_ij,k be the ordinary unbiased sample covariance of the two binary calls
over the r_k shared fragments. Then

\[
s^R_{ij,k}=\frac{r_k}{r_k-1}
(\widehat p_{11,k}-\widehat p^R_{i,k}\widehat p^R_{j,k}),
\qquad
\widehat c_k=\frac{r_k}{N_{ik}N_{jk}}s^R_{ij,k}.
\]

This simplification follows directly from the definition of sample covariance.
It identifies the proposal as an overlap-scaled estimate of error covariance.
Its possible contribution is a useful specialization and demonstrated behavior;
the displayed formula alone does not establish a new statistical principle.

The sign is not necessarily positive. Replace “inflates” and “systematically
overestimates” by “adds a sampling-covariance term, whose sign depends on the
shared-fragment joint states.” Positive and negative within-fragment coupling
should both be included in controlled tests. The existing pilot primarily uses
positive within-fragment coupling even when between-donor covariance is negative.

The result is also not assumption-free. It requires representative fragment
subsets, independent fragment sampling conditional on donor parameters and
coverage, and the specified binary-call model. The corrected covariance is
unbiased conditional on that model; the normalized correlation ratio is not
therefore unbiased.

## 4. Correct the three-way validation identity

Write P for the latent proportions of the n eligible donors and condition on
the relevant coverage/design information as well. Define

\[
T(P)=\frac1{n-1}\sum_k(p_{ik}-\bar p_i)(p_{jk}-\bar p_j).
\]

This is the latent covariance of the sampled donors. It is distinct from the
fixed superpopulation parameter theta = Cov(p_i,p_j).

Let U = corrected covariance from A, W = raw covariance from A, and R = the
B/C reference. Under the fragment model, E[R | P, design] = T(P). If R is
conditionally independent of U and W, expansion of squared errors gives

\[
E[(U-R)^2-(W-R)^2\mid P,design]
=E[(U-T)^2-(W-T)^2\mid P,design].
\]

This is the valid conditional risk identity. For a population parameter theta,
the corresponding unconditional identity instead includes an extra term:

\[
\Delta_{split}-\Delta_{population}
=-2E[(U-W)(R-\theta)].
\]

For this moment correction, let h(P,design) be the true average within-donor
sampling covariance in A. Then E[U-W | P,design] = -h, so the difference is

\[
2E[h\{T(P)-\theta\}].
\]

For representative iid donor sampling with E[T] = theta, this reduces to
2 Cov(h,T), which need not vanish: sampling covariance can itself depend on
donor methylation proportions. Thus independence of physical fragments does
not make the reference independent of donor-sampling fluctuations.

An exact finite enumeration verified the distinction. In a deliberately small
example, there are two iid donors; each has p_i = p_j = 0.1 with probability
0.8 or 0.5 with probability 0.2. Each partition has two independent fragments,
with perfectly concordant calls at the two sites. Enumerating four latent donor
configurations and 3^6 possible partition counts gives:

| Risk difference: corrected minus raw | Exact-enumeration value |
|---|---:|
| Against the B/C reference | 0.0023085 |
| Against latent covariance T of these donors | 0.0023085 |
| Against population covariance theta | 0.0010797 |
| Difference between reference and population risk differences | 0.0012288 |

The last number equals 2 Cov(h,T). This is an algebra check, not a proposed
two-donor application. For the same iid latent law and depth, the extra term
is 0.0024576/n, so it also persists for the proposed n >= 20 range.

Retain the primary benchmark as a conditional, finite-cohort measurement-error
risk comparison. Evaluate superpopulation covariance risk in simulations where
theta is known, or design a separate reference with independently sampled
donors from a suitably matched population. A new external cohort is still
necessary to assess generalization of a locked procedure. It does not by
itself change which parameter the within-cohort split benchmark evaluates.

Disjoint partitions are independent only under the generative fragment model
and appropriate conditioning. Randomly splitting a fixed, already observed
library does not create independent libraries or remove all systematic errors.

## 5. Correct the assay and counting rules

- The GTEx colon dataset in this project is RRBS, not WGBS.
  [GEO GSE233417](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233417).
- Generic coordinate deduplication must not be imposed on RRBS. Bismark
  explicitly advises against it for RRBS and related enrichment libraries;
  shared alignment positions can correspond to valid observations. Use a
  library-specific policy, exploiting UMIs when available, and document
  residual molecule dependence.
  [Official guidance](https://felixkrueger.github.io/Bismark/usage/deduplication/).
- N and r count fragments with usable methylation calls. Merely spanning a
  CpG or both CpGs is insufficient if a call is missing or excluded.
- Counting overlapping mates once is necessary; physically merging sequence
  files is not the only way to achieve that counting rule.
- Use shared-subset marginals in the unbiased covariance formula. If shared
  and nonshared fragments actually sample different methylation populations,
  the proposed latent-target model can fail even with that formula implemented
  correctly. The shared-marginal choice does not repair selection bias.
- A minimum of 20 donors is an eligibility choice, not a guarantee of numerical
  stability or adequate power. Equal weighting likewise does not fix donor
  ascertainment or informative coverage selection.

## 6. Wording and status corrections

Use “signed moment estimate that may violate parameter bounds,” rather than
claiming these estimates have mathematically unbounded ranges. Inputs are
bounded and counts are constrained; the latent population variance is in
[0,1/4], covariance in [-1/4,1/4], and squared disagreement in [0,1]. Distinguish
these parameter bounds from sample-moment bounds and estimator status codes.

The draft's table saying three-way validation has been performed must say
“proposed.” Likewise, achieving unbiased biological inference has not been
demonstrated in real data. The checked pilot uses 20,000 Monte Carlo cohorts per
scenario; a proposed 2,000-cohort expanded grid is a separate choice from the
2,000 resamples proposed for uncertainty estimation.

The joint donor/genomic bootstrap and 1-Mb block size remain proposed and must
be calibrated. They do not automatically handle every source of dependence.
Each pair should be assigned to one block using a fixed rule and retained as
one observational unit across comparisons.

A more precise title is **“Correcting Shared-Fragment Sampling Error in
Inter-Donor CpG Covariance.”** Avoid suggesting removal of biological
confounding. Replace the novelty conclusion with:

> We propose a count-based specialization of correlated measurement-error
> correction for adjacent CpGs. Its equivalence to prior estimators, practical
> benefit, interval calibration, and independent-cohort performance remain to
> be established.

The bibliography needs paper-level citations attached to individual claims.
Generic search-result labels, unrelated datasets, and broad field summaries
cannot establish absence of an existing estimator.
