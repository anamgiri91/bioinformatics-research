# Task 67: draft endpoint for the co-methylation decay analysis

**Superseded on 2026-10-02.** The implemented development analysis uses an
absolute contrast and count/sequence-based matching, replacing the ratio
and observed methylation/spread matching below. See the
[current report](Task67_DistanceDecay_Report.md). This historical draft
was never locked and must not be used as the current external protocol.

Draft written 2026-10-01. **It is not locked.** This is step 2 of the next
steps in `plan_shared_read_noise.md`: pre-specify a standard analysis that
the correction could change, before evaluating it on new data. It must
be frozen, with hashes, before it is applied to any cohort not yet
examined.

## Why

Tasks 53, 59 and 65 showed descriptive distance profiles. In plasma, for
example, the observed covariance fell 84.8% from 0 to 10 bp to 150 to
200 bp, the corrected one 42.4% and the B/C reference 47.3%.

Codex's audit pointed out that depth, overlap, CpG density and
methylation level all differ between distance bins. So these falls mix
the shared-read noise with pair composition. A fair comparison has to
hold the composition fixed.

## Question

Once pairs are matched on composition, how much does subtracting the
shared-read term change the estimated fall of co-methylation with
distance? Is the corrected fall closer to the fall in a reference that
shares no fragments?

## Estimands, per distance bin

- **Raw:** mean across-donor covariance of the A betas.
- **Corrected:** the same, after subtracting the Task 53 noise term
  (A only).
- **Reference:** mean B-with-C cross covariance. This requires partitions
  assigned per whole fragment (Task 66 rules), so that B and C share no
  physical fragment.

## Matching on composition

The bins are (0,10], (10,20], (20,40], (40,60], (60,100], (100,150]
and (150,200] bp. Pairs in each bin are reweighted to the joint
distribution, in the (0,10] bin, of four covariates:
- mean depth in A, in quartiles;
- local CpG density (CpGs within 100 bp, in quartiles);
- CpG island status (inside, shore, other);
- mean methylation level (5 bins) by spread (3 bins).

Cells with fewer than 20 pairs in either bin are dropped from both, and
the dropped share is reported.

## Primary endpoint

For each estimate, the decay ratio is the weighted mean covariance in
(150,200] over that in (0,10]. The endpoint has two parts:

    Delta_decay = (corrected ratio - reference ratio)^2
                - (raw ratio - reference ratio)^2

Negative means the corrected decay is closer to the reference decay.
Also reported: the three ratios, and the full matched curves.

## Uncertainty and decision

- **Interval.** The joint donor and 1-Mb-block bootstrap, re-centred by
  the bootstrap mean shift, as in Task 58. It must first be calibrated
  in simulation for this ratio endpoint, using the same three-scenario
  design as Task 57 plus a scenario with real decay.
- **Primary criterion.** The upper end of the 95% interval for
  Delta_decay is below 0.
- **Secondary results.** Each bin's matched covariance under the three
  estimates. Neither the dropped-cell share nor any change to the weights
  can rescue a failed primary result.

## Data

- **Main use:** the first locked test of a new cohort (Task 68
  recommends GSE165915). This endpoint would be a pre-declared secondary
  result there.
- **Plasma:** on the 46 plasma donors rebuilt with linked mates (Task 66),
  it can only be a post-result sensitivity analysis, because those
  donors have already been examined.

## Still to settle before the lock

- The CpG island annotation source and build: UCSC cpgIslandExt for hg19
  or hg38, chosen to match the cohort.
- Whether to use a standardised measure, such as covariance divided by
  the product of the reference SDs. Correlations failed in Task 63, so
  covariance stays primary.
- The simulation scenarios for the calibration, written down before it
  is run.
