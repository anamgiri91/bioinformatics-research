# Plan

Ordered by how much each step changes a conclusion, not by effort.
Status as of 2026-09-03, after Tasks 29-32.

Context for every item below: [REVIEW.md](REVIEW.md) has the code-level
findings, [README.md](README.md) the results so far, [LITERATURE.md](LITERATURE.md)
the review of `OutlierMeth`'s reference panels and parameters, Borealis and
epimutacions.

---

## Where the project stands

The original question — *do outlier flags track biology, or the methylation
state the CpG already sits in?* — now has an answer, and it is not the answer
the by-state count tables were built to find.

The **rate** of flags is roughly flat across states. The **magnitude** is not:
98% of external-reference flags at `H` sites and 82–87% at `L` sites
correspond to a beta shift under 0.10, against a median shift of 0.21–0.27 at
`R` sites.

Task 22 established the mechanism and closed off the alternatives:

1. **Not panel miscalibration.** If the BRCA cohort simply sat outside the
   747-sample `tcga` panel's range at these CpGs, many samples per CpG would flag.
   They do not — at most 11 of 53, zero CpGs with ≥50% flagged, 0% of Normal
   flags in such CpGs. The flags really are per-sample calls.
2. **Not a join or register error.** Reconstructing the thresholds from the
   flag matrix found zero CpGs where a flagged sample's beta sat inside the
   unflagged range, across 13,618 CpG × tissue combinations.
3. **It is threshold geometry.** At `H` sites the hyper-outlier zone is 0.037
   wide and the hypo zone 0.886, yet 87.8% of flags are `+1` — a 171× higher
   flag density per unit of available beta space in the constrained direction.
4. **The boundary is inside the noise.** At `H` sites the two samples
   straddling the threshold are 0.003 apart. At a perturbation of sd = 0.01,
   inside the array's own technical error, 21% of flags vanish and more calls
   change than there were flags to begin with.
5. **It is not an `OutlierMeth` bug.** Tukey's 3×IQR rule shows a state bias
   too — a different one, driven by the IQR itself ranging from 0.008 at `L`
   to 0.157 at `R` — and overlaps the percentile rule at Jaccard 0.05–0.29.

Tasks 29–32 close off the remaining escape routes:

6. **It cannot be tuned away.** The `p` argument has two usable settings on the
   747-sample `tcga` panel, not four — at n = 747 the three strictest levels all
   interpolate between the two most extreme reference samples. And because
   thresholds are monotone in `p` the flag set is nested, so the whole reachable
   family can be enumerated: at the strictest setting that still flags anything,
   chr22 Normal `H`-site `+1` flags are *still* 100% below |Δβ| = 0.10, with 0
   of 536 surviving a 0.10 floor. `p` changes how many flags you get, never how
   big they are.
7. **The noise level is now measured, not assumed.** Adjacent-probe pairs bound
   technical noise at 0.0029 (Normal) / 0.0036 (Tumour) at the tenth percentile
   and 0.014 / 0.017 at the median. The simulated sd = 0.01 sits inside that
   range, so the reproducibility argument no longer rests on a literature value.
8. **It is not a few bad probes, and it is not a few bad samples.** The contrary
   flags concentrate on sites (half the sites of a state carry all of them) but
   spread across samples (31 of 53 in Normal). The concentrating sites are the
   ones where the external threshold happens to land closest to this cohort's
   median, at 1.0–1.9 cohort MADs — geometry again, not biology, which is why
   probe masking will move the pattern rather than remove it.

That reframes the deliverable. The paper is not "flags are biased toward state
X". It is **"scale-free thresholding on bounded, heteroscedastic beta values
manufactures irreproducible outliers at the methylation extremes, and here is
what it costs you"** — with a concrete, cheap fix that applies to the whole
class of methods, not one package.

Two things are still needed before it can be written: the probe-artifact
confound has to be excluded (step 2), and the geometry and stability results
have to be shown to hold beyond chr22 (step 3, 7c.2). The `"Rc"` bug is fixed
and the magnitude result is already genome-wide.

---

## 1. ~~Rerun Tasks 10, 13, 14 with `"Rc"` → `"R"`~~ — DONE 2026-08-28

All three scripts are fixed and now carry a guard that stops on any
`methy.state` present in the data but absent from `state_order`. Corrected
tables were regenerated in Task 23 directly from the Task 13 flag matrices — no
`referenceMeth()` rebuild, no `tcga.rda` needed, and the five surviving rows
reproduce the originals to the digit.

The bug was **worse than estimated here**: framed as "40.3% of tumour CpGs", it
was 65.0% of tumour *flag mass* in the external-reference arm. The deleted `R`
rows have the highest flag rate of any state (21.6%) and the largest effect
sizes. See `results.md` §8.

The original text is kept below for the record.

**Blocking. Cheap. Changes a headline number.**

`state_order <- c("L","LM","M","HM","H","Rc")` in
`Task10.FlagStateTable.Aug6.2026.R:42`,
`Task13.FlagBias.StateAndSample.NoExtRef.Aug2026.R:24` and
`Task14.FlagBias.StateAndSample.ExtRef.Aug2026.R:24`. The data has `"R"`. The
downstream `tab[ord, ]` keeps only matched states, so every `R` CpG is dropped
silently — 6.2% of Normal, **40.3% of Tumour**.

Task 14's tumour external-reference table is the strongest existing result in
the project and was computed with the single largest category of tumour CpGs
missing. Task 21 now shows those `R` sites carry 19,170 of the chr22 tumour
external flags at a median |Δβ| of 0.214 — the highest-magnitude flags in the
dataset are exactly the ones being discarded.

Do not just fix the string. Derive the order from the data and assert the
total, so this class of bug cannot recur:

```r
state_order <- intersect(c("L","LM","M","HM","H","R"), sort(unique(df$methy.state)))
stopifnot(sum(tab) == 380355)
```

Tasks 15, 18, 20 and 21 already use `"R"` and need no change.

## 2. Mask cross-reactive and polymorphic probes

**Blocking for any biological claim.**

Nothing in the pipeline has excluded probes that are known to be unreliable:
cross-hybridising probes, probes overlapping common SNPs at the CpG or single-
base-extension site, and probes on repetitive sequence. Apply the standard
450k masks (Chen et al. 2013; Zhou, Laird & Shen 2017 / `sesameData`).

This is not housekeeping. Two current results depend on it:

- N37's focal event sits at chr22:16,601,097–16,602,837 — the pericentromeric
  region of chr22's acrocentric short arm, repetitive and poorly mapped. A
  seven-probe hypermethylation run there is exactly what a mapping artifact
  looks like. Until masked, it cannot be reported as an epimutation.
- The samples carrying the most burden (N15, N17, N14) may be carrying probe
  failures rather than methylation.

The Task 20 annotation already checks CpG-island and SNV overlap for contrary
flags; in the tight window both came back zero. Extend it chromosome-wide and
add the masking lists.

## 3. Confirm the magnitude result outside chr22

**The main result rests on one chromosome.**

Task 21's magnitude table is chr22 only (6,809 CpGs). The Task 13 flag
matrices already cover all 380,355 sites for both references and both tissues,
so this is a straight extension of Task 21 Part E — no new flagging, just a
wider denominator. Report median |Δβ| and the fraction under 0.05 / 0.10 /
0.15 per state, genome-wide, both tissues.

If the pattern holds at 380k sites it is a general property of percentile
thresholding on beta values, which is a much stronger claim than a chr22
observation.

## 3b. ~~Get a real noise estimate instead of a simulated one~~ — DONE 2026-09-03

Option 2 is done (Task 32 Part A). For probe pairs within 100 bp the systematic
difference between the two positions is constant across samples and drops out of
`sd(beta_i − beta_j)`, leaving `sd/sqrt(2)` as an upper bound on technical noise.
chr22, tenth percentile over pairs: **0.0029 (Normal), 0.0036 (Tumour)**; median
over pairs 0.014 / 0.017. Distant pairs (>1 Mb) give 0.0081 / 0.0107 at the same
percentile, 2.8–3.0× looser, so the bound is measuring probe-level noise rather
than biology.

**The assumed sd = 0.01 sits inside the measured range and is conservative
relative to the median, so every stability result in §7 and Tasks 24/26/32
stands as published.** `Results/Task32_MeasuredNoise.csv`

Option 1 — technical replicates — would still be tighter and is the only way to
separate technical noise from genuine divergence between neighbouring probes.
It needs the TCGA barcodes upstream of de-identification. Low priority now that
the assumption is validated.

**A negative result worth keeping:** scaling the magnitude floor by this
measured sigma makes things *worse*, not better. `ext.floor.noise` (floor =
k × sigma per state) scores 0.584 / 0.705 stability against the constant 0.10
floor's 0.793 / 0.907, because sigma is smallest exactly where the state is most
compressed, so the `L`-site floor lands at 0.023 and keeps the trivial flags.
Noise scaling and interpretability are different requirements — see
`results.md` §13.

## 4. Re-flag with an absolute-difference floor and quantify what changes

**Mostly done.** Task 26 swept `offset ∈ {0, 0.05, 0.10, 0.15, 0.20}` and
Task 32 ran a head-to-head of five floor variants at one flag rate. Summary:

| variant | Normal stability | Tumour | contrary-flag share, Normal |
|---|---|---|---|
| no floor (as published) | 0.329 | 0.648 | **52.9%** |
| `\|beta − median\| ≥ 0.05` | 0.748 | 0.871 | 28.0% |
| **`\|beta − median\| ≥ 0.10`** | **0.793** | **0.907** | 21.0% |
| per-state floor, calibrated | 0.781 | 0.751 | 14.9% |
| `k ×` measured noise, per state | 0.584 | 0.705 | 15.7% |
| `deltMeth ≥ d` (distance to threshold) | 0.649 | 0.880 | 17.2% |

The flat 0.10 floor on `|beta − cohort median|` wins in both tissues.
`epimutacions` uses 0.15; the 0.10–0.15 gap is small (`REVIEW.md` Q9).

**Still to do: the per-sample question.** Whether N15/N17/N14 stay on top once
trivial flags are removed has not been run. If they do, their burden is real;
if they collapse, it was heteroscedasticity all along. Task 29 gives a partial
hint — over the 100-CpG window's contrary flags the leaders are N15, N48 and
N17, so N48 needs adding to the list of samples to explain — but the
chromosome-wide burden ranking under a floor has not been recomputed.

**Do not spend time on M-values.** Task 24 tested them: a median ± k·MAD rule
on M-values is the least stable of seven methods benchmarked at matched flag
rates (Jaccard 0.308 Normal / 0.474 Tumour, against 0.558 / 0.628 for the same
rule on raw beta, and 0.899 / 0.953 for the external reference plus a 0.10
floor). The logit stabilises variance but amplifies noise at the extremes,
which is where the defect lives. Recorded as a negative result in
`results.md` §9; an earlier draft of this plan recommended the opposite.

## 5. Report event burden, not flag burden

The per-sample summaries count correlated neighbours as independent outliers.
N37's eleven flags are three events; its largest is seven CpGs across 1,740 bp.

Task 21 already computes `events.1kb`, `epimutations` (runs of ≥3, the
`epimutacions` definition) and `max.run.cpgs` at all three scales. Make the
event count the headline number in every per-sample table from here on, with
the raw flag count kept alongside. Where the two diverge, the raw count is
double counting.

Consider building co-methylated regions properly (CoMeBack) rather than using
a fixed 1 kb gap — the gap is defensible but arbitrary, and roughly 10% of
450k probes fall in genuinely co-methylated blocks.

## 6. Test what drives the per-sample burden

N15, N17 and N14 lead at every scale in both tissues (Spearman 0.761 Normal /
0.843 Tumour between chr22 external burden and genome-wide self burden). That
is a global sample property and it is currently unexplained. The candidates,
in the order they should be excluded:

1. **Technical** — array batch, plate, scan date, detection-p failure rate,
   bisulfite conversion controls. McCartney et al. found outlier burden
   "strongly associated with cell proportions and other technical factors";
   this is the first thing a reviewer will ask.
2. **Composition** — for tumour, purity/ploidy (ABSOLUTE or ESTIMATE); for
   normal, immune and adipose fraction by reference-based deconvolution.
3. **Biological** — global hypomethylation, age, subtype.

Until (1) and (2) are excluded, per-sample burden cannot be interpreted.

## 7. Decide the fate of the self-reference arm

Currently: `bio` and `self` are retained in Tasks 18/20/21 as a negative
control, clearly labelled. That is honest and it is enough for the paper —
the self-vs-external contrast is itself a methodological caution worth
publishing, with the order-statistic argument as the explanation.

If a genuinely reference-free arm is wanted, it must stop being an order
statistic of the same sample. Measured at n = 53 over 3,000 heterogeneous CpGs
([REVIEW.md](REVIEW.md) Finding 1):

| method | per-CpG count range | distinct values | variance |
|---|---|---|---|
| current (self, empirical quantile) | 2–2 | 1 | 0.000 |
| leave-one-out | 2–4 | 3 | 0.471 |
| parametric (per-CpG mean ± z·SD) | 0–5 | 6 | 1.692 |

Parametric is the only option that lets a CpG have **zero** outliers, which is
what you want, since most CpGs should have none. `bio.loo` is already
implemented in Tasks 18–21 as the minimal-change fallback. Prefer a robust
parametric form — median ± z·MAD, or a beta-distribution fit, which
`epimutacions` also offers and which respects the [0,1] bound.

**Retire Task 11.** It asks whether removing samples 14–17 changes the mean
per-CpG flag count, a statistic pinned at exactly 2 by construction. It cannot
answer its own question and cost ~28 minutes of compute doing so. Delete or
mark superseded; do not rerun under the current design.

## 7b. Items raised by the 2026-08-28 supervisor review

Answered empirically in Task 26 and written up in `REVIEW.md` Part 1:
denominator asymmetry, floor sensitivity, rate-matched stability, genome-wide
magnitude, the `outliers.coef2/3` definition, and the matched-pairs test.
Four items remain open.

1. **Confirm the normal/tumour pairing from the TCGA barcodes.** A
   methylation-based identity test puts `cor(N_i, T_i)` as the row maximum in
   only 10 of 53 cases on the most variable CpGs — above chance, but far from
   what genuine pairs give. The filtered 380k set appears to drop the `rs`
   identity probes. Nothing currently uses pairing, but the phrase "53 matched
   pairs" is inherited, not verified.

2. **Establish that the external panel is technically comparable.** Array
   version, preprocessing, normalisation, batch structure, and the definition
   of "tumour-adjacent" are all unverified. Gross miscalibration is ruled out
   (no CpG flags more than 11 of 53 samples), but that is not comparability.
   This is the largest un-addressed threat to the external arm.

3. **Add direction-specific agreement and precision/recall to §6.** Cohen's
   kappa is unreliable under >95% zeros. Task 20 already writes the full
   contingency counts, so this is an aggregation over existing output: do the
   methods agree on `+1` and `−1` separately, and what is precision/recall
   against a designated reference method? Task 32 adds one direction-aware
   metric (`pct.contrary`, the share of a method's flags running against the
   site's state) but not the pairwise direction-split agreement this asks for.

4. **Uniform provenance blocks.** Every table names its output CSV and every
   CSV has one producing script with a `.Rout` transcript, but the review asked
   for a consistent block — script, inputs, output, denominator, reference
   mode, parameters — under each of ~20 tables and 5 figures. Formatting, not
   analysis.

Also extend the §7 threshold-geometry and noise-stability experiments beyond
chr22. The magnitude result is now genome-wide; these two are not.

## 7c. Items opened by the 2026-09-03 email and Tasks 29-32

1. **The `bio.stat.dev` arm needs a decision.** Task 32 makes the state-pooled
   deviation rule the best method tested that needs no external panel at all
   (0.779 / 0.833 stability, state-rate ratio 1.0, zero contrary flags by
   construction). That is a genuinely useful result — it means the paper's
   conclusion does not depend on trusting a pan-tissue reference panel. It
   should either become the second recommended rule or be explicitly retired;
   leaving it as a third unlabelled column repeats the `bio`/`self` confusion.

2. **Extend Tasks 29–31 beyond chr22.** Task 31's ceiling, `p`-envelope and
   `deltMeth`/`relMeth` results are chr22-only. The Task 13 matrices cover all
   380,355 CpGs, so this is a wider denominator on existing code, as
   step 3 was.

3. **Mask the sites Task 29 identified first.** The contrary flags concentrate
   on a handful of probes — `cg15668074` alone carries 19 of the 29 tumour `LM`
   `−1` flags, and the Normal `LM` block sits at chr22:16.60–16.61 Mb, the same
   pericentromeric region as N37's candidate event. Probe masking (step 2) will
   remove some of these. It will not remove the pattern: Task 29 Part B shows
   the concentrating sites are the ones where the external threshold happens to
   sit closest to this cohort's median, so once they are masked the same
   concentration will reappear at the next-closest sites. Report both.

4. **N48.** Task 29 puts N48 level with N15 as the top carrier of contrary
   flags in the normal window (9 each). N48 is not in the N15/N17/N14 trio that
   §5/Task 21 identified from chromosome-wide burden. Either it is a
   window-local artefact or the burden list is incomplete — one query against
   `Task21_Chr22_SampleBurden_Normal.csv` settles it.

5. **Say which panel, everywhere.** The reference is `tcga` (747 samples,
   21 tissue types), not `all` (2,015 / 25). Corrected in the four root docs on
   2026-09-03. Any figure caption, slide or draft carried over from before that
   date needs the same fix.

## 8. Housekeeping

- Fix Task 12's denominator — the self-referential mean is taken over
  NA-omitted CpGs, the external mean over all 380,355 rows, and the ~10,154
  rows absent from `tcga.rda` enter as zeros. Restrict to
  `intersect(rownames(beta), rownames(tcga))`. The reported means are ~2.7%
  low; the 5× tumour/normal conclusion is unaffected.
- Task 17's agreement metric counts shared zeros and is near-total by
  construction. Superseded by Task 20's kappa and non-zero conditional
  agreement — mark it so rather than citing the old number.
- Delete `Scripts/fill_Task13_excel.R`; it is byte-identical to
  `Scripts/Task13.CreateSampleSummary.R`.
- Add `stopifnot(identical(...))` on the CpG key before the positional `cbind`
  in `Experiment6.Combine.July31.2026.R:26-36`.
- Add `--no-save` to the eight batch scripts that omit it (the `.RDataTmp`
  errors ending the Task 10/11/12 and Experiment 6/7/8 logs). Fixed from
  Task 13 onward.
- Parameterise `DATA_DIR` / `OUT_DIR` in the pre-Task-18 scripts. They embed a
  colleague's home directory (`/home/s_s355/...`), which matters if the repo
  is ever public. Tasks 18–21 already resolve paths at runtime.
- `git rm -r --cached Packages/OutlierMeth` and record the dependency at
  commit `408b110` instead — see [CLASSIFICATION.md](CLASSIFICATION.md).
- Reattach the Aug 6, 2026 meeting-note headers lost to the rsync:
  `git show b9b1281:Scripts/Task14.FlagBias.StateAndSample.ExtRef.Aug2026.R`.

---

## Two questions for the next meeting

1. **Is the deliverable the methodological result?** The state-magnitude bias
   plus the absolute-difference fix is a cleaner, more defensible paper than
   the by-state flag-rate comparison the tables were originally built for, and
   it stands on the external-reference arm alone. Confirm the target before
   step 3 is run genome-wide.

2. **Is the `bio` rule meant to be a detector or a contrast?** As specified it
   cannot produce a contrary flag — it only ever assigns `+1` at `L`/`LM` and
   `-1` at `H`/`HM` — so a contrary-flag count for it is structurally zero and
   its comparison against the external reference measures the definitional
   difference, not a disagreement about the data. If it is meant as a
   detector, it needs the parametric threshold of step 7 and a magnitude
   floor. If it is meant as a contrast that shows what state-aware flagging
   would look like, it is already doing its job and should be labelled that
   way in the paper.

3. **Which of the two fixes is the recommendation?** `ext` + a 0.10 floor is
   the most stable rule tested and keeps the external panel. `bio.stat.dev` is
   only slightly behind, needs no panel at all, and cannot produce a contrary
   flag by construction. They are answers to different questions — "how do I
   make this published method usable" versus "what should replace it" — and the
   paper can carry both, but the abstract can only lead with one.
