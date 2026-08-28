# Code & Methods Review — Outlier Methylation Pipeline

Reviewed 2026-08-27 against the full tree pulled from
`leap2.txstate.edu:/mmfs1/home/wln26/Experiments.Outlier.July31.2026`.
Findings are ordered by how much they change a conclusion.

---

## What the pipeline does

BRCA TCGA 450k methylation, 380,355 CpGs across 53 matched normal/tumor pairs
(plus a 32-sample "dead" cohort). Each CpG carries a precomputed methylation
state (`L`, `LM`, `M`, `HM`, `H`, `R`). The `OutlierMeth` package (Downs & Cope,
JHMI) flags each CpG x sample cell as `-1` (hypo), `0`, or `+1` (hyper) by
comparing its beta value against per-CpG thresholds.

Thresholds come from one of two places, and the distinction is the spine of the
whole project:

- **Self-reference** — `referenceMeth()` derives thresholds from the same 53
  samples being tested. Experiments 7, Tasks 10, 11, 13, 15, 17.
- **External reference** — the prebuilt `tcga.rda` panel (~2,000 independent
  samples). Tasks 12, 14, and the `_extref_` half of Task 13.

The scientific question is whether flags track real biology or just track the
methylation state a CpG already sits in — i.e. whether calling a low-methylation
site "hypomethylated" is a finding or a tautology.

**That question is well posed, and the external-reference arm answers it. The
self-reference arm cannot, for the reason in Finding 1.**

---

## Finding 1 — CRITICAL: the self-reference arm is a mathematical identity

**Every self-referential result in this project is a fixed constant determined
by the sample size. None of it carries biological signal.**

`referenceMeth()` sets the hyper threshold at `quantile(x, 0.99)` over the 53
values at a CpG. R's default type-7 quantile places that at position
`(n-1)p + 1 = 52 x 0.99 + 1 = 52.48` — strictly between the 52nd and 53rd order
statistics. Exactly one sample can exceed it: the maximum. The 0.01 quantile
lands at position 1.52 and exactly one sample falls below it: the minimum.

So **every CpG receives exactly one `-1` and one `+1`, always**, regardless of
its methylation state, its variance, or anything biological.

Simulated at n = 53 over 2,000 random CpGs per significance level:

| p | hyper flags/CpG | hypo flags/CpG | total |
|---|---|---|---|
| 0.01 | 1 | 1 | 2 |
| 0.001 | 1 | 1 | 2 |
| 0.0001 | 1 | 1 | 2 |
| 0.00001 | 1 | 1 | 2 |

Note the `p` column does nothing. Every threshold from 0.01 to 0.00001 lands in
the same inter-order-statistic gap, so the significance level is inert at
n = 53. Any result reported "at p = 0.01" is identical to one at p = 0.00001.

The predicted constants are `1/53 = 1.89%`, `51/53 = 96.23%`, and `2/53 = 0.0377`
mean flags per CpG. The recorded outputs match exactly:

- `task13_stateFlagPct_normal_noExtRef.csv` — every state reads
  `1.89 / 96.23 / 1.89`, to the digit.
- `Task15_Normal_flagbias_by_state.csv` — `mean.neg1 = 1`, `mean.0 = 51`,
  `mean.1 = 1` for all six states.
- Task 12 log — self-referential mean per-CpG count `1.999982` (Normal),
  `1.999997` (Tumor). Both are 2.

### What this invalidates

- **Task 13 / 15 Part A ("no flag bias across methylation states")** is not a
  negative result. The tables were incapable of showing a difference; identical
  rows were guaranteed before the data was read.
- **Task 11 ("does removing samples 14-17 change anything?")** reports
  `1.999982 -> 1.999971` for Normal and `1.999997 -> 1.999997` for Tumor. The
  statistic is pinned at 2 by construction, so removing samples could not have
  moved it. This experiment cannot answer its own question — and it burned
  ~28 minutes of compute across four `referenceMeth()` builds to do so.
- **Task 17's biological rule** re-implements the same `quantile(values, 0.99)`
  over the same 53 values, so `biological.flag.count` is 0 or 1 for every CpG
  (71 ones, 29 zeros — the 29 being `M`/`R` sites with no rule assigned).

### What survives

**Per-sample analysis is still informative.** Each CpG contributes exactly one
min and one max, but *which* sample supplies them varies, so per-sample totals
carry real signal — `Task15_Normal_flagbias_by_sample.csv` ranges from 4,423
(N1) to 9,328 (N3) hypo flags. Task 13/15 **Part B** is sound. It is only the
**by-state** analysis that is degenerate.

### Fix

**Correction to an earlier draft of this review: increasing the cohort does not
fix this.** The number of values above an empirical quantile is fixed by `n` and
`p` alone — it is a property of order statistics, and the data never enters into
it. Simulated over 5,000 CpGs spanning flat, low-skewed, high-skewed and bimodal
distributions:

| n | distinct per-CpG flag counts across 5,000 heterogeneous CpGs |
|---|---|
| 53 | one value: `2` |
| 201 | one value: `4` |
| 1000 | one value: `20` |

A larger cohort raises the constant; it never creates variation. **No
self-referential empirical-quantile design can answer a per-CpG or per-state
question at any sample size.** The by-state result is unobtainable this way, full
stop — which is precisely why the external-reference arm exists and is the one
that works.

If a self-referential arm is still wanted, the threshold has to stop being an
order statistic of the same sample. Measured at n = 53 over 3,000 heterogeneous
CpGs:

| Method | Per-CpG count range | Distinct values | Variance |
|---|---|---|---|
| Current (self, empirical quantile) | 2-2 | 1 | 0.000 |
| Leave-one-out | 2-4 | 3 | 0.471 |
| Parametric (mean +/- z*SD per CpG) | 0-5 | 6 | 1.692 |

Leave-one-out restores some variation by removing the guaranteed self-flag.
Parametric thresholds restore the most, and are the only option that lets a CpG
have **zero** outliers — which is what you actually want, since most CpGs should
have none. Recommended order:

1. **Report the external-reference results.** They are unaffected and already
   computed. The self-vs-external contrast is still worth publishing as a
   methodological caution, with this artifact as the explanation.
2. **Parametric per-CpG thresholds** if a reference-free arm is required.
3. **Leave-one-out** as a minimal-change fallback.

---

## Finding 2 — HIGH: a typo silently deletes an entire methylation state

`Task10.FlagStateTable.Aug6.2026.R:42`,
`Task13.FlagBias.StateAndSample.NoExtRef.Aug2026.R:24`,
`Task14.FlagBias.StateAndSample.ExtRef.Aug2026.R:24`:

```r
state_order <- c("L", "LM", "M", "HM", "H", "Rc")   # data has "R", not "Rc"
```

The subsequent `tab[ord, ]` and `summ[match(ord, summ$state), ]` keep only
matched states, so every `R`-state CpG is dropped without warning.

Row counts confirm it. Task 13/14 report `104,540 + 56,636 + 36,785 + 97,094 +
61,845 = 356,900` CpGs; Task 15 (which does not filter) reports those five plus
`R = 23,455`, totalling exactly 380,355.

| Dataset | R-state CpGs dropped | Share of data |
|---|---|---|
| Normal | 23,455 | 6.2% |
| **Tumor** | **153,479** | **40.3%** |

**This matters most for the result you'd actually publish.** The Task 14 tumor
external-reference table — currently the strongest finding in the project — was
computed with the single largest category of tumor CpGs missing.

`Task13.FullFlagMatrix.Aug20.2026.R` and `Task15` already document and avoid
this. Tasks 10, 13, and 14 have not been rerun since.

**Fix:** change `"Rc"` to `"R"` and rerun Tasks 10, 13, 14. Better, derive the
order from the data (`sort(unique(...))`) and `stopifnot()` that the row count
still totals 380,355 — a guard that would have caught this immediately.

---

## Finding 3 — MEDIUM: Task 12's headline comparison mixes denominators

`Task12.ExternalReference.Aug6.2026.R:19-31`. The self-referential mean is taken
over `rownames(ref.self.normal)` (NA-omitted CpGs only), but the external mean
is taken over all 380,355 rows. `tcga.rda` covers 370,201 of them; the ~10,154
unmatched rows come back all-`NA`, and `rowSums(abs(...), na.rm = TRUE)` turns
each into a `0` that is then averaged in as though it were a real measurement.

The external means are therefore diluted by `370,201 / 380,355 = 0.973`:

| | reported | over matched CpGs only |
|---|---|---|
| Normal | 1.351 | 1.388 |
| Tumor | 6.505 | 6.683 |

The conclusion — tumors carry ~5x the outlier burden of normals against an
independent panel — is unaffected and looks real. The numbers are ~2.7% low.
Restrict to `intersect(rownames(beta), rownames(tcga))` before averaging.

---

## Finding 4 — MEDIUM: Task 17's agreement metric counts shared zeros

`Task17.Chr22.BiologicalValidation.Aug21.2026.R`:

```r
sum(!is.na(a) & !is.na(b) & a == b)
```

`a == b` is satisfied when both flags are `0`, which is the overwhelmingly
common case: `bio_flag` is 0 for every `M`/`R` CpG outright and for 52 of 53
samples elsewhere. Agreement is therefore near-total by construction — the
observed distribution runs 44-53 out of 53, with 18 CpGs at a perfect 53 and 32
at 52. It measures the shared sparsity of two sparse matrices, not concordance.

Use Cohen's kappa, or condition on cells where at least one method flagged:
`sum(a == b & (a != 0 | b != 0)) / sum(a != 0 | b != 0)`.

---

## Finding 5 — LOW: positional `cbind` with no key check

`Experiment6.Combine.July31.2026.R:26-36` joins the 53-sample alive cohort to
the 32-sample dead cohort with `cbind`, assuming both files carry identical CpG
row order. Nothing verifies it. Correct today; a silent, invisible mis-join the
first time either file is re-sorted or re-exported. Add:

```r
stopifnot(identical(alive.normal$Composite.Element.REF,
                    dead.normal$Composite.Element.REF))
```

---

## Finding 6 — LOW: housekeeping

- `Scripts/fill_Task13_excel.R` is byte-identical to
  `Scripts/Task13.CreateSampleSummary.R`. Delete one.
- Six batch scripts omit `--no-save`, producing the `.RDataTmp` permission
  errors that end the Task 10/11/12 and Experiment 6/7/8 logs: `Task9.sbatch`,
  `Task10.sbatch`, `Task11.sbatch`, `Task12.sbatch`,
  `Experiment7.July31.2026.slurm`, `FastTasks.345.July31.2026.slurm`,
  `Task6.July31.2026.slurm`, `Task8.July31.2026.slurm`. Cosmetic — results were
  written before the failure — and already fixed from Task 13 onward.
- Scripts hardcode `/mmfs1/home/wln26/...` and `/home/s_s355/...` throughout.
  One `DATA_DIR` / `OUT_DIR` pair at the top of each would make the repo
  runnable by anyone and stop leaking a colleague's account name.

---

## Note on the sync

The 2026-08-27 rsync overwrote the local Task 13 and Task 14 scripts with the
leap2 copies. The leap2 versions are authoritative for *what actually ran*
(they carry the absolute output paths), but they had lost their documentation
headers — including the Aug 6, 2026 meeting notes that recorded why each
analysis exists. No executable logic was lost. The headers are preserved in
commit `b9b1281` and can be reattached with:

```
git show b9b1281:Scripts/Task14.FlagBias.StateAndSample.ExtRef.Aug2026.R
```

Worth restoring — that rationale is the clearest statement of the project's
motivation anywhere in the tree.

---

## Suggested order of work

1. Fix `"Rc"` -> `"R"`; rerun Tasks 10, 13, 14. Cheap, and it changes a headline
   number for tumor.
2. Decide the self-reference question (Finding 1). It determines whether Tasks
   11, 13A, 15A, and 17 get rerun under a valid design or get retired.
3. Fix the Task 12 denominator and the Task 17 agreement metric.
4. Reattach the lost headers; parameterise the paths.
