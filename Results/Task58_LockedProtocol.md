# Task 58: locked protocol for the external shared-read noise test

Version 2, locked on 2026-10-01 before any CpG-pair result from the external
samples was computed. It follows plan sections 8 to 10 of
`plan_shared_read_noise.md`. SHA-256 hashes of this file, the code and the
input files are in `Results/Task58_LockedProtocol.json`. The test script
(`Scripts/Task59.SharedReadNoise.ExternalTest.Oct01.2026.R`) refuses to run
if any of them has changed.

## 1. Data

- **Cohort.** The 57 white-blood-cell RRBS samples of GEO GSE233417, chosen
  from metadata (`Results/Task56_CohortChoice.md`). All 57 passed the
  independence screen (`Results/Task56_ExternalSamples.csv`):
  - same-person control pairs scored 0.975 to 0.997;
  - all pairs of different people scored 0.32 to 0.49;
  - no blood sample matched another blood sample or a colon donor.
- **Files.** The 57 `.pat.gz` files, hashed in the JSON. Donors are ordered
  by GEO accession.
- **Build and coordinates.** hg19. The wgbstools genome CpG index is mapped
  to positions from the UCSC hg19 sequence
  (`Data/hg19_seq/cpg_positions_hg19.rds`, hashed).
- **Molecules.** UMI-deduplicated per GEO. One `.pat` count is one fragment
  with its mates merged (results.md section 16.2).

## 2. Pairs and eligibility

- **Pairs.** Adjacent genome CpGs (index i and i + 1) on autosomes, at most
  200 bp apart.
- **Pair universe.** This is fixed by counts only. A pair is included if at
  least 20 donors have at least 6 reads at both sites, which is necessary
  for the split rule.
- **Donor eligibility.** A donor counts for a pair if both sites have at
  least 2 reads in each of A, B and C, and A has either 0 or at least 2
  shared fragments.
- **Pair eligibility.** A pair enters the endpoint if it has at least 20
  eligible donors.

## 3. The three-way split

Every `.pat` row's count is split between A, B and C by a multinomial draw
with probability 1/3 each. All calls of a fragment stay together. The seed
is 20261102 + 1000 x donor rank + chromosome number. This seed series is
new and was never used in development.

## 4. Estimators and reference

- **Observed covariance** s_A: the sample covariance of the A betas across
  eligible donors.
- **Corrected covariance** Chat_A: s_A minus the mean shared-read noise term
  from A's joint counts (`Scripts/SharedReadNoise.Core.R` formulas, as
  implemented in `Scripts/SharedReadNoise.Benchmark.R`).
- **Reference** R: (s(b_i^B, b_j^C) + s(b_i^C, b_j^B)) / 2. B and C share no
  fragments.

## 5. Primary endpoint and decision

**Endpoint:**

    Delta_MSE = mean over eligible pairs of [(Chat_A - R)^2 - (s_A - R)^2]

Negative values favour the correction.

**Interval (version 2).** There are 2,000 joint bootstrap draws (seed
20261103):
- donors are resampled whole;
- 1-Mb blocks are resampled within chromosomes, with block = (position - 1)
  div 1,000,000 of the lower CpG;
- the draws use a fixed subsample of the pairs whose CpG index is divisible
  by 25.

The interval is the **bias-corrected percentile** interval: the percentile
interval shifted by the mean of the draws minus the subsample estimate.

**Why version 2.** The plan proposed the plain percentile interval. In
development simulations (Task 57) it covered the true endpoint only 38% to
82% of the time, because resampling donors adds variance to the squared
errors. The bias-corrected version covered it 96% to 100%, and never
falsely showed a benefit in a null scenario. The plan says such a change
must be chosen on development data and recorded before external testing,
which this does.

**What the interval targets.** It targets the average endpoint over donors
and genomic blocks, as in Task 57. It is not an interval conditional on
these 57 donors.

| Result | Conclusion |
|---|---|
| Upper end below 0 | Primary criterion met: the corrected covariance is closer to an independent reference in this external setting |
| Interval crosses 0 | Improvement not established |
| Lower end above 0 | The correction has higher error |

## 6. Secondary outcomes

These cannot rescue a failed primary result.

- The same endpoint for squared disagreement (corrected against observed,
  with the reference Qref).
- The relative change in squared error, and the share of pairs where the
  correction is closer.
- Results by distance ((0,10] to (150,200] bp, as in the plan), by mean
  depth in A, and by overlap in A.

The plan's pair-ranking comparison is not part of this version. It was not
developed.

## 7. Deviations from the plan, stated before the test

1. **Interval.** The bias-corrected percentile replaces the plain
   percentile (section 5).
2. **Bootstrap subsample.** The bootstrap runs on a fixed 1-in-25 subsample
   of pairs, to keep it computable. The point estimate on all pairs is
   reported alongside.
3. **Pair universe.** It is built from counts only, before the split.
4. **Same lab.** The external donors are independent people, but the samples
   come from the same study, lab and pipeline as the development data.
5. **Limits from earlier sections still apply.** These are the provenance
   and assumption limits in results.md sections 16.2 and 16.4. In
   simulations, PCR copies, nonrepresentative shared fragments and
   conversion errors bias the correction.

## 8. Software

R 4.6.0 with data.table. Code hashes are in the JSON.
