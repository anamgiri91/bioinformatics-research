# Item 5 plan: correcting shared-fragment sampling error in inter-donor CpG covariance

Updated: 2026-10-01. **Status: development plan, carried out through phase H
on 2026-10-01 (Tasks 53 to 60).** The locked external test (Task 59) met its
primary criterion in one external cohort from the same study. The method is
still not established as new, and the section 6 comparisons of correlation,
pair ranking and existing error models are not done. Results are in
[results.md](results.md) section 16. The outlier-analysis roadmap stays in
[plan.md](plan.md); this file is the separate plan for the shared-read noise
work.

The research question is whether per-donor joint methylation states on shared
DNA fragments improve estimation of latent CpG covariance and agreement across
donors, beyond raw estimates and existing measurement-error methods.

The primary deliverable is an audited count pipeline, explicit moment
estimators, diagnostics, and a calibrated evaluation. A useful new application
and a new statistical estimator are different possible outcomes. Neither is
assumed in advance.

## 1. Current evidence and remaining work

| Component | Current status | Implication for this plan |
|---|---|---|
| Frozen score and validation | Task 47 protocol and subsequent analyses already exist | Preserve the original hypotheses and results; this proposal is separate |
| GTEx colon cohort | GSE233417, RRBS, 29 donors; already examined | Use for development, not as untouched validation of a new method |
| Shared-fragment moment derivation | Checked by exact enumeration and a simulation pilot | Supports the algebra under its assumptions, not real-data utility |
| Existing simulation pilot | Eight scenarios, 20,000 simulated cohorts of 29 donors per scenario | Extend to negative within-fragment coupling, calibration, and misspecification |
| Correlation stability | The pilot found nonpositive corrected variances and out-of-range ratios | Retain explicit failures; do not present a ready-made bounded coefficient |
| Per-donor joint-count pipeline | Not implemented and audited for this proposal | Reconstruct the required counts; pooled linkage is insufficient |
| Independent cohort and three-way benchmark | Not selected or completed | Select, freeze, and evaluate prospectively |
| Novelty | Exact equivalence to prior methods remains unresolved | Use potential-contribution language |

Relevant project records:

- [Original feasibility input](Data/NN.hg38.18P.forw.chr22.w.header.txt) and
  [earlier CpG correlation work](correlation%20for%20cg%20sites/).
- [Frozen protocol](Results/Task47_FrozenProtocol.md) and
  [completed GTEx validation report](Results/GTExColon_Validation_Report.md).
- [Moment derivation and pilot](Results/Item5_SharedReadNoise_Design.md),
  [pilot script](Scripts/Item5.SharedReadNoise.Pilot.R), and
  [pilot summary](Results/Item5_SharedReadNoise_Summary.csv).
- [Prior-art assessment](Results/Item5_SharedReadNoise_PriorArt.md) and
  [mathematical and methodological audit](Results/Item5_SharedReadNoise_DraftAudit.md).

Beta-only tables can support comparisons but cannot supply shared-fragment
counts. The existing `Data/gtex_colon_rrbs/cache_cohort.rds` and pooled-linkage
outputs must not be assumed to retain everything this method requires.

## 2. Define the targets before implementing the estimator

For pair (i,j), let I_ij be its eligible donor set and n its size. Let p_ik and
p_jk be donor k's underlying bulk methylation proportions. All formulas below
use the same eligible donors and equal donor weights.

### Primary measurement target: the sampled donors

Define the latent sample covariance and variance with denominator n-1:

\[
T_{ij}=\frac{1}{n-1}\sum_{k\in I_{ij}}
(p_{ik}-\bar p_i)(p_{jk}-\bar p_j),\qquad
T_{ii}=\frac{1}{n-1}\sum_{k\in I_{ij}}(p_{ik}-\bar p_i)^2.
\]

Define latent squared disagreement and agreement with denominator n:

\[
Q_{ij}=\frac1n\sum_{k\in I_{ij}}(p_{ik}-p_{jk})^2,
\qquad A_{ij}=1-Q_{ij}.
\]

The three-way fragment benchmark evaluates estimation error about T_ij and
Q_ij, conditional on these donors and the sampling design. The underlying
agreement A_ij lies in [0,1]. A finite-sample estimate may leave that interval.

### Separate population target

The population quantities are theta_ij = Cov(p_i,p_j), Var(p_i), and
E[(p_i-p_j)^2]. Connecting sample targets to these parameters requires a
specified population and representative donor sampling. Coverage-based
eligibility can change the represented population. Equal donor weighting does
not repair selection bias.

The primary split benchmark is not automatically an MSE comparison about
theta_ij. Evaluate population risk separately in simulations with known
population parameters. Any empirical claim about population risk needs an
additional justified design, such as an independent-donor reference from a
matched population.

### Biological interpretation and initial genomic scope

- Target bulk-tissue co-methylation. Cell-mixture covariance remains part of
  the bulk target; direct regulatory coupling is not identified.
- Start with autosomal, consecutive reference-genome CpGs separated by at
  most 200 bp. Consecutive measured sites are not necessarily consecutive
  reference CpGs.
- Use continuous beta values. Ten-level discretization is a secondary
  comparison and is not used to estimate the sampling correction.
- Use genomic distance for diagnostic strata, not as a score multiplier.
- Preserve each dataset's genome build and coordinate conventions. Record
  any annotation liftover; do not calculate gaps between mixed builds.

## 3. Sampling model and its limits

Condition on the donor proportions, coverage, and overlap design. The moment
identities require:

1. Fragments are independent draws within a donor, after appropriate molecule
   handling. Donors are independent for the initial inferential model.
2. Shared and nonshared fragments at a site represent the same underlying
   methylation proportion. Coverage or calling selection must not invalidate
   this representativeness assumption.
3. Calls are binary observations of methylation under the primary model.
   Conversion errors, genotype effects, and other misclassification require
   separate assessment; subtracting sampling covariance does not fix them.
4. The joint state distribution of a shared fragment can have positive,
   zero, or negative covariance. Physical overlap alone does not determine
   the direction of the correction.

Without an additional error model, persistent calling errors can make the
latent target a call probability rather than true methylation probability.
Describe any such interpretation explicitly. These are conditional moment
results, not assumption-free guarantees.

## 4. Build and audit the fragment counts

For every donor and CpG pair retain this record:

| Field | Definition |
|---|---|
| Donor ID, genome build, chromosome, positions | Stable pair and sample identifiers |
| N_ik, N_jk | Fragments with usable calls at each site |
| M_ik, M_jk | Methylated calls among those fragments |
| r_k | Fragments with usable calls at both sites |
| n_00,k, n_01,k, n_10,k, n_11,k | Joint calls on exactly the r_k shared fragments |
| Provenance and QC | Source file, library, processing version, available mapping/base-quality information |
| Partition ID, when applicable | Whole-fragment A/B/C assignment |

The first joint-state digit refers to site i. Check:

\[
r_k=\sum_{a,b\in\{0,1\}}n_{ab,k},\qquad
0\le r_k\le\min(N_{ik},N_{jk}).
\]

Also require valid methylated totals and nonshared methylated counts:

\[
0\le M_{ik}-(n_{10,k}+n_{11,k})\le N_{ik}-r_k
\]

and the corresponding condition for site j. Compare reconstructed site counts
and coordinates with independently generated site summaries where available.
Investigate disagreements rather than silently accepting a permissive match.

### Fragment identity and duplicate policy

Count overlapping mates once at each CpG, retaining all usable CpGs on their
physical fragment. Record how discordant mate calls are resolved. Geometrically
spanning a site is insufficient when the methylation call is missing.

The development cohort is RRBS. Do not impose generic coordinate deduplication:
[Bismark's guidance](https://felixkrueger.github.io/Bismark/usage/deduplication/)
advises against that operation for RRBS and related enrichment libraries.
Choose molecule handling from the actual library design and use verified UMIs
where available. Shared coordinates alone do not identify PCR copies in RRBS.
When molecular independence cannot be verified, record the limitation and
include dependent-fragment sensitivity analyses.

Audit PAT generation and multiplicities. Confirm whether each counted unit is
a fragment, a read end, or another summary. Do not interpret one compressed
pattern row as one molecule. A PAT-based split is acceptable only if its
provenance establishes how counted fragments and all their calls are retained;
otherwise obtain appropriate alignment data or mark fragment validation unmet.

### Initial eligibility

- Require N_ik >= 2 and N_jk >= 2.
- Accept r_k = 0 or r_k >= 2. Mark r_k = 1 unsupported for this estimator.
- Retain a pair with at least 20 eligible donors. This threshold is a design
  choice, not evidence of stability or adequate power.
- Use identical donors and pairs for matched comparisons. Report every
  exclusion and the resulting coverage of sites, donors, and chromosomes.
- For the three-way benchmark, enforce eligibility jointly across A/B/C as
  specified below, using count/QC criteria rather than the observed score.

## 5. Implement the moment estimators

For each donor set b_ik = M_ik/N_ik and b_jk = M_jk/N_jk. For r_k >= 2 compute
the shared-subset proportions:

\[
\widehat p_{11,k}=\frac{n_{11,k}}{r_k},\quad
\widehat p^R_{i,k}=\frac{n_{10,k}+n_{11,k}}{r_k},\quad
\widehat p^R_{j,k}=\frac{n_{01,k}+n_{11,k}}{r_k}.
\]

Estimate sampling variance and covariance:

\[
\widehat v_{ik}=\frac{b_{ik}(1-b_{ik})}{N_{ik}-1},\qquad
\widehat c_k=\frac{r_k^2}{N_{ik}N_{jk}(r_k-1)}
\left(\widehat p_{11,k}-\widehat p^R_{i,k}\widehat p^R_{j,k}\right).
\]

Equivalently, if s^R_ij,k is the ordinary unbiased sample covariance of binary
calls on the shared fragments,

\[
\widehat c_k=\frac{r_k}{N_{ik}N_{jk}}s^R_{ij,k}.
\]

This identifies the correction as an overlap-scaled measurement-error
covariance estimate. Do not substitute all-read betas for the shared-subset
marginals. Set chat to zero for r = 0 under the model; do not evaluate the
shared-state formula at r = 0 or 1.

Let s_ij, s_i^2, and s_j^2 be observed sample moments over the same eligible
donors, with denominator n-1. Return:

\[
\widehat C_{ij}=s_{ij}-\frac1n\sum_k\widehat c_k,\qquad
\widehat V_i=s_i^2-\frac1n\sum_k\widehat v_{ik},\qquad
\widehat V_j=s_j^2-\frac1n\sum_k\widehat v_{jk}.
\]

Under the stated model, their conditional expectations are T_ij, T_ii, and
T_jj. For disagreement and agreement return:

\[
\widehat D^2_{ij}=\frac1n\sum_k
\left[(b_{ik}-b_{jk})^2-\widehat v_{ik}-\widehat v_{jk}
+2\widehat c_k\right],\qquad
\widehat S_{ij}=1-\widehat D^2_{ij}.
\]

Their conditional expectations are Q_ij and A_ij. This is a correction for
squared disagreement, not for Manhattan absolute disagreement.

The secondary diagnostic ratio is

\[
\widehat\rho_{ij}=\frac{\widehat C_{ij}}
{\sqrt{\widehat V_i\widehat V_j}}.
\]

Unbiased moment ingredients do not make this ratio unbiased or well behaved.
Keep raw moments and diagnostic status codes. If either corrected variance is
nonpositive, do not report a valid correlation. If both are positive but the
ratio lies outside [-1,1], retain the diagnostic value and mark it invalid.
Never silently clip it, replace it with zero, or omit its failure rate.

Signed moment estimates may violate parameter bounds; they are not a validated
bounded score. A constrained or regularized estimator is a separate candidate
requiring its own specification and comparison. Pairwise outputs using
different donor sets also do not automatically form a valid positive
semidefinite covariance matrix.

## 6. Compare estimators with matching targets

| Evaluation | Primary comparisons | Required reporting |
|---|---|---|
| Covariance recovery | Raw s_ij versus corrected Chat | Bias, MSE, signed correction, uncertainty |
| Correlation recovery | Pearson; variance-only s_ij/sqrt(Vihat Vjhat); full correction | Error on defined targets and failures across all eligible pairs |
| Squared disagreement | Raw mean squared difference versus corrected Dhat^2 | Bias, MSE, negative/out-of-range estimates |
| Pair ranking | Manhattan similarity, Spearman, dCor, Chatterjee's xi, frozen scores | Common held-out endpoint and selection rule |
| Existing error models | Adapted correlated measurement-error estimator where justified | Assumptions, inputs, convergence, parameter tuning |

Variance-only correction is not a separate covariance baseline because it
leaves s_ij unchanged. Correlation and dependence measures such as dCor or xi
do not all estimate the same parameter; do not compare their numerical values
as if each were a covariance estimate.

For the prospective agreement-ranking comparison, rank pairs using A only and
evaluate their average B/C reference squared disagreement. Proposed selection
fractions are 1%, 5%, 10%, and 25%, with fractional weighting at score ties.
Use the same pair universe, report valid-score coverage, and label this a
secondary utility comparison rather than a test of identical estimands.
Resolve undefined-score handling and the exact xi implementation in development
before locking this comparison.

Adaptations of existing error models and any regularization choices must be
documented and fixed before external evaluation. Mathematical equivalence is
an algebra/model comparison; empirical superiority is an estimation or utility
comparison. Neither establishes the other.

## 7. Development and simulation program

Use the already examined GTEx colon data to establish count feasibility,
runtime, eligibility, and correction magnitudes. Do not tune on the external
cohort. Report correction distributions by depth, overlap fraction, distance,
and observed methylation spread. These diagnostics do not estimate latent
truth in real data.

### Controlled simulation grid

The expanded simulator must vary both between-donor structure and
within-fragment coupling, which are separate quantities.

| Axis | Proposed coverage |
|---|---|
| Donors | 20, 29, 50, 100 |
| Depth | 2, 3, 5, 10, 30; asymmetric 5/30; heterogeneous donor depths |
| Shared fraction | 0, 0.25, 0.5, 1; report unsupported r = 1 cases |
| Latent structure | Constant, independent, positively related, inversely related, low spread, extreme means, systematic mean offset |
| Within-fragment association | Negative, zero, positive, and feasible-boundary cases |
| Target | Conditional sample moments and population moments, scored separately |

For given p_i and p_j, use feasible joint probabilities. One explicit
construction interpolates between p_i p_j and either Frechet bound:

\[
q_{11}=(1-|\eta|)p_i p_j+|\eta|q_{bound},\quad
q_{bound}=\begin{cases}
\min(p_i,p_j),&\eta\ge0,\\
\max(0,p_i+p_j-1),&\eta<0.
\end{cases}
\]

Set q_10 = p_i-q_11, q_01 = p_j-q_11, q_00 = 1-p_i-p_j+q_11. Include eta =
-1, -0.5, 0, 0.5, 1 where relevant. Generate shared calls jointly and nonshared
calls separately. Verify valid probabilities and all count identities.

Use a proposed 2,000 Monte Carlo cohorts per selected scenario for the expanded
screen, with fixed seeds and Monte Carlo standard errors. Preserve the existing
20,000-replicate pilot as a separate result. A machine-readable scenario list
must specify latent distributions and every evaluated combination; the axis
table alone is not a frozen simulation protocol. Select a declared subset for
the more expensive nested-bootstrap calibration, before seeing its coverage.

Assess:

- Bias and MSE of covariance and disagreement relative to known targets.
- Correlation errors only when the true correlation is defined, together with
  all estimator failures and the fraction of outputs deemed usable.
- No-overlap behavior: the covariance correction is exactly zero. Marginal
  variance and disagreement corrections need not vanish.
- Signed within-fragment covariance, heterogeneous depths, and finite-cohort
  versus population risk differences.
- Interval coverage and false-positive rates, with Monte Carlo uncertainty.

### Assumption violations

Add separate experiments for PCR-family dependence, nonrepresentative shared
fragments, incomplete and inappropriate conversion, state-dependent call
loss, informative coverage, and variable cell mixtures. Treat conversion rates
as declared sensitivity inputs rather than measured characteristics of GTEx.
Report failures without redefining the target after seeing them. If a coverage
restriction or model revision is needed, select it on development data and
record a new protocol version before external testing.

### Inference calibration

Start with 2,000 whole-donor percentile-bootstrap draws for population-moment
intervals under the independent, representative-donor model. Evaluate coverage
against known population parameters in repeated-donor simulations. These are
not automatically conditional intervals for the fixed donors' T_ij or Q_ij.

If conditional moment intervals are needed, specify and calibrate a separate
fragment-level resampling or model-based procedure that preserves the overlap
structure and declared coverage conditioning. The initial release need not
claim such intervals. Do not create correlation intervals by dropping invalid
bootstrap ratios. A calibrated ratio interval remains a separate requirement.

## 8. External cohort selection and protocol lock

Select an independent cohort from metadata before inspecting its CpG-pair
results. Initial eligibility is at least 30 donors from one tissue, with raw
fragment information or independently prepared replicate libraries and enough
coverage to support the planned split. Verify independence from development
donors, including duplicate submissions or reused libraries.

The accession has not been selected. Record candidates, access restrictions,
assay, tissue, available depth metadata, replicate design, provenance, and the
reason for choosing one. If none qualifies, report external validation unmet.
Do not substitute a convenient development subset and call it external.

Before reading validation outcomes, record a versioned protocol containing:

- Accession and sample list; input checksums; reference genome and coordinates.
- Library/molecule handling, quality filters, eligibility, and pair definition.
- Estimator code, baseline implementations, any development-selected settings.
- Scenario list and calibration findings; interval construction and scope.
- Fragment-assignment algorithm, seed, primary split, and sensitivity splits.
- Pair/block assignments, diagnostic bins, score-tie and missing-output rules.
- Primary endpoint, pass criterion, secondary outcomes, and reporting template.
- Software versions and code hashes.

Use a single declared primary split. Additional seeds or block sizes are
sensitivity analyses with complete reporting. A post-outcome implementation
fix requires a transparent deviation record and retention of affected outputs.

## 9. Three-way physical-fragment benchmark

Randomly assign each physical fragment within a donor to A, B, or C with
probability 1/3 using a reproducible algorithm and fixed seed. All CpG calls and
both mates stay in that partition. The assignment is global within each
library, not redrawn separately for each pair. Where PCR families are known,
handle them as specified by the audited molecule policy to prevent copies of
one original molecule supporting multiple partitions.

Compute all candidate estimates and pair rankings using A only. Use B/C only
for reference construction and count-based eligibility. Require N >= 2 at both
sites in each partition and r_A = 0 or r_A >= 2, with at least 20 common
eligible donors per pair. B and C need not each satisfy an r >= 2 condition
because their cross-partition reference does not use chat.

Construct the noisy references:

\[
R_{ij}=\frac12\{s(b_i^B,b_j^C)+s(b_i^C,b_j^B)\},
\qquad
Q^{ref}_{ij}=\frac1n\sum_k
(b_{ik}^B-b_{jk}^B)(b_{ik}^C-b_{jk}^C).
\]

Under the generative fragment model and appropriate conditioning,
E[R | P,design] = T and E[Qref | P,design] = Q. Here P denotes the eligible
donors' latent proportions. The reference is not error-free ground truth.

Let U = Chat^A and W = s_ij^A. Conditional independence of A from B/C and an
unbiased reference give:

\[
E[(U-R)^2-(W-R)^2\mid P,design]
=E[(U-T)^2-(W-T)^2\mid P,design].
\]

This is the justification for the benchmark. It concerns the covariance of
the sampled donors. If theta is a population covariance, replacing T by theta
in the unconditional risk identity introduces the term

\[
\Delta_{split}-\Delta_{population}
=-2E[(U-W)(R-\theta)]
=2E[h\{T-\theta\}],
\]

where h is the true average shared-sampling covariance in A. Under iid donor
sampling with E[T] = theta this becomes 2 Cov(h,T), which need not be zero.
The [audit](Results/Item5_SharedReadNoise_DraftAudit.md) gives a checked example.

Disjoint read assignments do not eliminate PCR dependence, library-wide
artifacts, selection bias, or donor-sampling fluctuations. Splitting a fixed
observed library is not equivalent to obtaining independent libraries.
Independently prepared libraries can strengthen the technical validation,
with their target and error assumptions specified separately.

## 10. Primary endpoint, uncertainty, and decision rules

Average equally over the eligible pair universe:

\[
\Delta_{MSE}=\frac1m\sum_{(i,j)}
\left[(\widehat C_{ij}^A-R_{ij})^2-(s_{ij}^A-R_{ij})^2\right].
\]

Negative values favor correction. The expected statistic compares conditional
finite-cohort covariance-estimation risks averaged over pairs. Its units are
squared covariance units. Also report raw losses, relative changes where
denominators are meaningful, and the distribution of pairwise differences.

Use a proposed 2,000-replicate joint bootstrap of whole donors and genomic
blocks for uncertainty in this aggregate endpoint. Start with nonoverlapping
1-Mb blocks within chromosomes. Assign each pair once using the lower CpG's
coordinate: block = floor((position_1based-1)/1,000,000). Resample blocks within
chromosomes and donors consistently across both methods and all partitions;
recompute both estimators and the reference in each draw. Keep the primary
eligibility mask fixed and report any bootstrap-specific estimation failures.

The target and calibration of this joint bootstrap must be checked in
development simulations with donor and spatial dependence. It is a proposed
procedure, not guaranteed valid by its block size. Millions of adjacent pairs
are not millions of independent replicates. Define the spatial and donor
sampling interpretation in the lock record.

**Primary improvement criterion:** the upper endpoint of the prespecified
95% interval for Delta_MSE is below zero, using the calibrated procedure on
the untouched external cohort.

| Outcome | Permitted conclusion and next action |
|---|---|
| Primary criterion met, assumptions and calibration supported | Evidence of improved finite-cohort covariance estimation in this external setting |
| Interval crosses zero | Improvement not established; report precision and effect size |
| Correction has higher error | Report loss of performance and its coverage/spread pattern |
| Calibration or provenance fails | Endpoint cannot support the planned inferential claim; return to development |
| No qualifying external cohort | External validation remains unmet |

Agreement, ranking, and correlation results are secondary and cannot rescue a
failed primary endpoint. A significant primary result does not prove novelty,
population-level risk improvement, or universally reliable correlation.

Report distance bins (0,10], (10,20], (20,40], (40,60], (60,100], (100,150],
and (150,200] bp. Also report depth, overlap, methylation spread, annotation,
and chromosome strata with denominators. Freeze exact additional cutpoints
using development data. Avoid selecting only favorable strata after testing.

## 11. Prior art and the potential contribution

| Work | Established contribution relevant here | Comparison still required |
|---|---|---|
| [MHL / Guo et al.](https://www.nature.com/articles/ng.3805) | Read-level methylation coupling and haplotype summaries | Distinguish its target from across-donor error correction |
| [MPCI](https://pmc.ncbi.nlm.nih.gov/articles/PMC13035127/) | Weighted Manhattan similarities across reads and CpGs | The broad Manhattan-based similarity idea is already established |
| [DNAmBERT](https://academic.oup.com/bib/article/27/5/bbag455/8780345) | Sequence and methylation-haplotype token modeling | Do not describe it as an estimator based solely on four pairwise counts |
| [MethylPCA](https://link.springer.com/article/10.1186/1471-2105-14-74) | Distinguishes shared-fragment and biological correlation in an MBD-seq workflow | Assess conceptual overlap despite assay and algorithm differences |
| [Buonaccorsi et al.](https://pubmed.ncbi.nlm.nih.gov/27264206/) | Bisulfite binomial measurement-error correction | Review full estimator details before asserting non-equivalence |
| [epiG](https://link.springer.com/article/10.1186/s13059-017-1168-4) | Error-aware methylation inference and independent-assay benchmarking | Differentiate validation of states from validation of covariance |
| [dSOMNiBUS](https://onlinelibrary.wiley.com/doi/10.1002/sim.10149) | Mismeasured multivariate methylation outcomes and dependence modeling | Compare assumptions and intended estimands |
| [Saccenti et al.](https://www.nature.com/articles/s41598-019-57247-4) and [MeasurementError.cor](https://bioconductor.org/packages/release/bioc/manuals/MeasurementError.cor/man/MeasurementError.cor.pdf) | Correlated-error treatment and latent-correlation estimation | Include relevant correlated-error baselines; independence is not a universal limitation |

The exact question is whether shared-fragment joint counts have already been
used to estimate this sampling covariance and correct across-donor CpG moments,
including mathematically equivalent formulations. Three-way partitioning also
needs its own prior-art comparison. Failure to locate a matching paper is not
proof that an estimator or benchmark is new.

Current wording for a proposal or manuscript:

> We investigate a count-based specialization of correlated measurement-error
> correction for adjacent CpGs. Per-donor shared-fragment joint states estimate
> a sampling-covariance term in observed across-donor beta covariance. Under
> explicit assumptions, a fragment-partitioned benchmark evaluates conditional
> finite-cohort estimation risk. Equivalence to prior methods, calibration,
> practical benefit, and independent-cohort generalization remain to be
> established.

## 12. Execution sequence and deliverables

The following are planned deliverables, not files claimed to exist.

| Phase | Work | Required output / gate |
|---|---|---|
| A. Resolve provenance | Audit library/PAT generation, calls, fragments, reference coordinates | Input manifest, duplicate policy, and pass/fail provenance report |
| B. Reconstruct counts | Stream per donor/chromosome; checkpoint counts; verify identities and positions | Per-donor joint-count cache and eligibility tables |
| C. Implement estimators | Keep count extraction separate from deterministic moment calculations | Pairwise estimates, raw ingredients, diagnostic codes, meaningful algebra checks |
| D. Simulate and compare | Expand model cases, existing-method comparisons, misspecification and calibration | Scenario manifest, bias/MSE/coverage tables, failure rates, runtime estimates |
| E. Develop on GTEx | Assess feasibility and actual correction magnitudes | Development report with complete diagnostics and selected settings |
| F. Select and lock | Choose independent cohort from metadata and finalize all pending specifications | Dated protocol, sample list, hashes, seeds, primary endpoint and analysis code |
| G. Evaluate externally | Build A/B/C once, compute references and paired endpoint, run calibrated uncertainty | Primary result, secondary tables, sensitivity results, deviation log |
| H. Assess contribution | Compare algebra, assumptions, software behavior, and practical performance | Bounded novelty statement and reproducible final report |

Meaningful checks include exhaustive small-count identities; positive and
negative shared-state covariance; zero/partial/full overlap; r = 1 and other
unsupported inputs; fragment assignments shared by mates and CpGs; coordinate
and genome-build consistency; constant-site failures; matched donor masks; and
the conditional versus population benchmark identity. Do not create tests
that merely reproduce the implementation without checking an independent
identity, reference calculation, or failure condition.

Store run manifests, input and code hashes, seeds, package versions, exclusion
counts, logs, and intermediate summaries. Stream large inputs rather than
duplicating whole cohorts in memory. Use new output locations for this project
stage and preserve frozen Task 47 scores and validation outputs as historical
records. Record deviations when fixes follow inspection of outcomes.

**Immediate next step:** audit the GTEx PAT/alignment provenance and determine
whether independent fragment units and per-donor joint counts can be recovered.
Then implement the count rebuild and deterministic estimators. External data
evaluation begins only after the development, calibration, and protocol-lock
gates above are satisfied.
