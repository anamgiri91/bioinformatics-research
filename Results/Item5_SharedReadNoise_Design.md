# Item 5 follow-up: separate shared-fragment noise from between-person covariance

Exploratory derivation and simulation, 2026-09-30. This is separate from the
frozen Task 47 score and its validation. **A novel correlation coefficient has
not been established.** The completed GTEx tests motivate a more precise target:
estimate agreement and covariance of latent methylation proportions while
accounting for fragments that contribute to both observed CpG betas.

The existing H9 result cannot determine whether its association arises from
between-person biology, coupled sampling error, or both. This simulation
demonstrates a possible mechanism, not an estimate of its contribution in GTEx.

## Established work that bounds the novelty claim

- Read-level CpG coupling and methylation haplotype blocks are established:
  [Guo et al., Nature Genetics (2017)](https://doi.org/10.1038/ng.3805) defined
  tightly coupled CpG blocks and a methylation haplotype load summary.
- Separating biological variability from measurement error in methylation is
  established. [A statistical method for excluding non-variable CpG sites
  (2010)](https://pmc.ncbi.nlm.nih.gov/articles/PMC2876131/) uses technical
  replicates to study signal relative to measurement variance.
- Correction for both correlated and uncorrelated errors is established:
  [Saccenti et al., Scientific Reports (2020)](https://doi.org/10.1038/s41598-019-57247-4)
  discusses correlation inflation/attenuation, correction, and estimates
  outside the correlation bounds.
- Count-based latent co-expression inference also has genomic precedents:
  [CS-CORE, Nature Communications (2023)](https://doi.org/10.1038/s41467-023-40503-7),
  with [official implementation](https://changsubiostats.github.io/CS-CORE/).
  Its RNA-count setting differs from overlapping binary methylation fragments.

The following equations specialize moment correction to shared binary fragments.
That is a mathematical derivation for this project, not a priority claim.
The targeted search did not establish whether this exact implementation has
already appeared. The potential contribution would need to be a demonstrated,
useful treatment of this sampling structure with valid uncertainty and an
independent benchmark, rather than the general act of subtracting noise.

## Observation model and identifiable moments

For CpGs i and j in donor k, let N_ik and N_jk be the numbers of called
fragments, r_k the number calling both, and b_ik and b_jk the observed betas.
The latent site proportions are p_ik and p_jk. Fragments are conditionally
independent; shared fragments have binary within-fragment covariance sigma_k.
All fragments covering a site must be representative of the same site
proportion. Counts and overlap are conditioned on.

Under these assumptions,

\[
v_{ik}=\operatorname{Var}(b_{ik}\mid p,N)
 =\frac{p_{ik}(1-p_{ik})}{N_{ik}},\qquad
c_k=\operatorname{Cov}(b_{ik},b_{jk}\mid p,N,r)
 =\frac{r_k\sigma_k}{N_{ik}N_{jk}}.
\]

For N>1 and r>=2, unbiased estimators are

\[
\widehat v_{ik}=\frac{b_{ik}(1-b_{ik})}{N_{ik}-1},\qquad
\widehat c_k=\frac{r_k^2}{N_{ik}N_{jk}(r_k-1)}
 (\widehat p_{11,k}-\widehat p_{i,k}^{R}\widehat p_{j,k}^{R}).
\]

The superscript R means proportions computed from the **same shared-fragment
subset**, not from all fragments at each site. The second formula follows
because the shared-fragment covariance with denominator r has expectation
(r−1)/r times sigma. When all fragments are shared, N_i=N_j=r, it reduces to
L/[2(r−1)] using the frozen within-donor linkage L. When r=0, the shared-fragment
term is zero under the model. With r=1 it cannot be estimated by this formula;
it must be flagged or modeled, not silently set to zero.

Let s_ij be the usual across-donor sample covariance with denominator n−1.
The moment estimates of latent covariance and variance are

\[
\widehat C_{ij}=s_{ij}-\frac1n\sum_k\widehat c_k,\qquad
\widehat V_i=s_i^2-\frac1n\sum_k\widehat v_{ik}.
\]

Conditional on the donor-specific latent proportions, their expectations are
the corresponding sample covariance/variance of those latent proportions.
This follows by expanding the centered sample covariance: independent,
mean-zero errors across donors contribute the mean within-donor error
covariance. The statement allows heterogeneous depths and error variances.

An agreement target with an analogous exact moment correction is

\[
\widehat D_{ij}^{\,2}=\frac1n\sum_k
 \{(b_{ik}-b_{jk})^2-\widehat v_{ik}-\widehat v_{jk}
                         +2\widehat c_k\}.
\]

Its expectation is n^−1 sum_k (p_ik−p_jk)^2. This is **squared** disagreement;
it is not a noise correction for Manhattan absolute disagreement. Taking a
square root or clipping negatives changes the expectation. Likewise, the
ratio C/sqrt(V_i V_j) is not an unbiased correlation estimator merely because
its moment ingredients are unbiased.

## Feasibility checks actually completed

`Scripts/Item5.SharedReadNoise.Pilot.R` enumerates every four-category outcome
for three shared fragments and checks the noise-variance, noise-covariance and
squared-disagreement identities to numerical precision. It then simulates
20,000 cohorts of 29 people in each of eight scenarios (seed 20260930):
constant, independent, positive linear and inverse latent proportions;
complete, partial and absent read overlap; equal and heterogeneous depths.
The within-fragment coupling parameter is 0.8 for these illustrative scenarios.
It was not estimated from GTEx, and the simulation is not calibrated to its
observed linkage distribution.

| Scenario | Raw covariance bias | Corrected covariance bias | Median observed Pearson r |
|---|---:|---:|---:|
| Constant latent proportions; 10 shared reads | 0.020034 | 0.000033 | 0.8073 |
| Independent latent proportions; 10 shared reads | 0.014638 | −0.000039 | 0.4048 |
| Independent; partial overlap, depths 20 and 30 | 0.002429 | −0.000014 | 0.1067 |
| Independent; no overlapping reads | −0.000003 | −0.000003 | approximately 0 |
| Independent; heterogeneous depths | 0.001457 | 0.000015 | 0.0737 |
| Positive linear latent relation; 10 shared reads | 0.018396 | 0.000052 | 0.8413 |
| Positive linear relation; heterogeneous depths | 0.001787 | −0.000019 | 0.7053 |
| Inverse latent relation; 10 shared reads | 0.013098 | 0.000048 | −0.0089 |

All corrected covariance and squared-disagreement biases are within five
Monte Carlo standard errors of zero. These are bias diagnostics, not Type I
error tests or evidence of biological performance.

The limitation is substantial: with independent latent proportions and ten
shared reads, 13.845% of cohorts have a nonpositive corrected variance and a
further 11.045% have a corrected correlation outside [−1,1]. With truly constant
latent proportions, a biological correlation is undefined, although the
observed correlation can be high. Boundary cases with a true correlation of
±1 also often produce out-of-bounds corrected ratios. Reporting only the valid
ratios, clipping them, or assigning zero to undefined cases would conceal
these failures and introduce additional bias.

[Full summary and Monte Carlo errors](Item5_SharedReadNoise_Summary.csv),
[simulation protocol](Item5_SharedReadNoise_Protocol.json),
[replicate output](Item5_SharedReadNoise_Replicates.csv.gz),
[run log](Item5_SharedReadNoise_Run.log).

## What a credible next method would require

1. Retain **per-donor** N_i, N_j, r, n_11, n_10 and n_01 from the PAT files. The
   current pooled-linkage cache loses the information needed for this correction;
   pooled mean linkage divided by pooled mean depth is not a valid substitute.
2. Estimate a constrained latent covariance model, or report signed moment
   estimates with uncertainty and explicit non-identifiability flags. Generic
   regularization is itself established and needs comparison, not a new name.
3. Use read counts and linkage from training donors only for pair selection.
   Validate on fragment-disjoint beta estimates or independent library
   replicates so same-fragment noise cannot support both sides of the endpoint.
   Splitting fragments must be done consistently across all CpGs on a fragment.
4. Compare raw Pearson/dCor/xi, depth-weighted agreement, variance-only
   correction, full shared-fragment correction and a binomial hierarchical
   model. Evaluate estimation error, calibration, missingness and power at
   matched actual false-positive rates.
5. Include violations: PCR dependence, cell mixture, nonrepresentative shared
   fragments, heterogeneous calling errors and variable donor composition.
   Cell-mixture covariance remains biological covariance at the bulk-sample
   level; this correction does not identify direct regulatory coupling.

The GTEx cohort has now been examined. Any subsequent tuning on it is
development, and an independent cohort would be needed for validation of a
new score. The completed result here is an explicit, checked mechanism and a
bounded research target. A publishable new method and its novelty remain open.
