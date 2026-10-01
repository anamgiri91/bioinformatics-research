# Task 64: locked protocol for the second external test

Locked on 2026-10-01, before any CpG-pair result from the plasma samples
was computed. It repeats the Task 58 test (version 2) on a cohort from a
different lab, assay and pipeline. It also adds the plan's section 6
comparisons as descriptive secondaries.

SHA-256 hashes of this file, the code and every input file are in
`Results/Task64_LockedProtocol.json`. The test script
(`Scripts/Task65.SharedReadNoise.ExternalTest2.Oct01.2026.R`) refuses to
run if any of them has changed. Its dry run on GTEx colon chromosomes 21
and 22 (`Scripts/Task65.DryRun.Rout`) tested the code before the lock.

## 1. Data

- **Cohort.** The 46 healthy-control plasma cfDNA samples of GSE149438,
  chosen from metadata (`Results/Task61_SecondCohortChoice.md`). All 46
  passed the independence screen (`Results/Task62_PlasmaSamples.csv`):
  - halves of one sample scored 0.953 to 1.000;
  - same-donor GTEx controls scored 0.920 to 0.993;
  - every pair of different people scored at most 0.593;
  - no sample matched another, or any development sample.
- **Files.** mHapBrowser `.mhap` files, hashed in the JSON. Donors are
  ordered by GEO accession.
- **Conversion.** Each line becomes a `.pat`-like row: the global hg19 CpG
  index of its first CpG, the 0/1 string as T/C, and its count. Both
  strands are kept. The script stops if any autosomal line is not on the
  hg19 CpG list or has the wrong length.
- **Molecules.** Duplicates were removed upstream (sambamba). Mates are
  merged when their CpG ranges overlap. Otherwise they are two lines (see
  section 7).

## 2. Pairs and eligibility

As in Task 58:
- adjacent genome CpGs on autosomes, at most 200 bp apart;
- the pair universe uses counts only: at least 20 donors with at least 6
  reads at both sites;
- a donor counts for a pair if both sites have at least 2 reads in each of
  A, B and C, and A has either 0 or at least 2 shared molecules;
- a pair needs at least 20 eligible donors.

## 3. The three-way split

Each line's count is split between A, B and C (multinomial, 1/3 each).
The seed is **20261302** + 1000 x donor rank + chromosome number. This
seed series is new.

## 4. Estimators, reference and primary endpoint

These are unchanged from Task 58:

    Delta_MSE = mean over eligible pairs of [(Chat_A - R)^2 - (s_A - R)^2]

with R = (s(b_i^B, b_j^C) + s(b_i^C, b_j^B)) / 2.

**Interval (version 2):**
- 2,000 joint bootstrap draws (seed **20261303**), resampling donors and
  1-Mb blocks within chromosomes;
- the bias-corrected percentile interval;
- the bootstrap pair set is every pair if the universe has at most
  300,000 pairs, otherwise every pair whose CpG index is divisible by 25.

| Result | Conclusion |
|---|---|
| Upper end below 0 | Primary criterion met in a second external setting, from another lab and pipeline |
| Interval crosses 0 | Improvement not established in this setting |
| Lower end above 0 | The correction has higher error in this setting |

## 5. Secondary outcomes

These cannot rescue a failed primary result.

- The squared-disagreement endpoint, the relative change in squared error,
  and the share of pairs where the correction is closer.
- Strata by distance, by mean depth in A and by overlap in A, as in Task
  58.
- **Sensitivity analysis on pairs at most 40 bp apart,** with its own
  bias-corrected interval from the same draws. This is declared because
  of the unmerged-mates limit (section 7). Pairs this close are rarely
  split across two lines.

## 6. Section 6 comparisons (descriptive)

These use the Task 63 definitions in `Scripts/SharedReadNoise.Section6.R`.
They are reported as found and have no pass rule.

1. **Correlation recovery.**
   - The estimators are Pearson, variance-only and the full correction.
   - The reference is R / sqrt(Ri Rj), with Ri = s(b_i^B, b_i^C).
   - An estimate is valid when its variances are positive and it lies in
     [-1, 1].
   - Reported: the valid shares, and the MSE against the reference on
     common pairs and on each estimator's own valid pairs.
2. **Pair ranking,** on the bootstrap pair set.
   - Ten scores from A only (the Section6 file lists them).
   - The top 1, 5, 10 and 25% of each score, with fractional weights at
     ties and undefined scores ranked last.
   - Each top set is judged by the mean reference squared disagreement
     Qref. This is done on all pairs, and on pairs where both sites vary
     in A.
   - Chatterjee's xi is the larger of its two directions, with ties broken
     by seed **20261305**.
   - **The main comparison, stated now:** the corrected positive part
     plus 1.645 SE, the observed squared disagreement, and the frozen
     F+. The noise-aware variants were added in development after the
     plain corrected score ranked noisy pairs first (Task 63).
3. **Existing correlated-error method.** Ding & Gentleman
   (MeasurementError.cor 1.84.0) on 2,000 random pairs with a defined
   reference (seed **20261304**), with se = sqrt(vhat) per donor. It is
   compared with the three estimators on common pairs.

**What development showed (Task 63, GTEx colon).** This was known before
the lock:
- corrected correlations were worse than Pearson;
- the observed squared disagreement and Manhattan agreement ranked best;
- the Ding & Gentleman baseline had the highest error.

## 7. Limits and deviations, stated before the test

1. **Unmerged mates.** When the mates' CpG ranges do not overlap, they are
   two lines. That fragment's shared calls are then hidden, and the split
   can send the two lines to different parts. For those pairs, the
   reference can carry some shared noise and the correction misses some.
   This mostly affects pairs farther apart than the mates' overlap. Hence
   the 40 bp sensitivity analysis.
2. **Deduplication by position, upstream.** cfDNA fragments often share
   ends, so some distinct molecules may have been merged. That lowers depth
   but creates no shared reads. Any PCR copies left would bias the
   correction (results.md 16.4).
3. **Library tail.** Swift Methyl-Seq adds a short tail to read 2. Whether
   it was trimmed is not stated. Untrimmed bases act like call errors.
4. **Biology close to blood.** cfDNA from healthy people comes mostly from
   blood cells, so this tests the lab, assay and pipeline more than the
   tissue.
5. **The bootstrap pair-set rule** (section 4) is new. It keeps every pair
   when the capture panel gives few enough.
6. **The 40 bp sensitivity analysis and the section 6 comparisons** are new
   relative to Task 58.

## 8. Software

R 4.6.0 with data.table, matrixStats, MeasurementError.cor 1.84.0 and
digest. Code hashes are in the JSON.
