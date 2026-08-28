# Review

Two reviews in one file.

- **Part 1** answers the supervisor's leading questions from the 2026-08-28
  critical review. Where a question was answerable from data, it was run rather
  than argued; where it was not, that is said plainly.
- **Part 2** is a readability audit of `results.md` against the standard "can a
  student with basic statistics follow this?"
- **Part 3** is the original code and methods review of 2026-08-27, updated
  with the current status of each finding.

The supervisor's overall assessment is accepted: the analysis is stronger than
the original Task 17 summary, and several statements were stated as conclusions
when the evidence supported hypotheses. Those have been revised in `results.md`
and the specific changes are listed at the end of Part 1.

---

# Part 1 — The leading questions

## Q1. What is the primary scientific question?

**Conceded — the report moved among three.** It is now committed to one:

> **Does `OutlierMeth`-style thresholding produce outlier calls that are
> reproducible and biologically interpretable, and if not, why not?**

That is a *methods-evaluation* question. Detecting real epimutations is a
downstream goal that this project is not yet positioned to claim, because probe
masking and technical QC have not been done. `results.md` §4's N37 language has
been softened accordingly.

## Q2. Should the `bio` rule be labelled a heuristic rather than a biological method?

**Yes, and it now is.** At n = 53, `quantile(x, 0.99)` under R's type-7 rule
falls at position `52 × 0.99 + 1 = 52.48`, strictly between the 52nd and 53rd
order statistics, so exactly one sample can exceed it. The rule assigns a fixed
flag count to every eligible CpG regardless of the data, at every significance
level from 0.01 to 0.00001.

It is labelled a **negative control** throughout: it shows what "no signal"
looks like. It is not a biological comparator, and no conclusion in `results.md`
rests on it.

## Q3. Are the 53 normal/tumour samples truly matched pairs?

**Cannot be confirmed from the data available.** Tested directly rather than
assumed. Correlating each `N_i` against every `T_j` over the most
between-individual-variable CpGs:

| CpGs used | times `cor(N_i, T_i)` is the row maximum | expected if unmatched |
|---|---|---|
| all 380,355 | 1 / 53 | ~1 |
| top 10,000 by variance | 4 / 53 | ~1 |
| top 2,000 by variance | **10 / 53** | ~1 |

10 of 53 is well above chance, so the column ordering probably does encode
pairing — but it is nothing like the near-perfect matching that genuine paired
samples give on identity probes. The filtered 380k set appears to exclude the
`rs` SNP probes that would settle it.

**Action:** confirm from the TCGA barcodes upstream of de-identification before
any paired analysis. **No result in `results.md` currently uses pairing**, so
nothing is invalidated — but the phrase "53 matched normal/tumour pairs" is an
inherited assumption, not something this project has verified.

## Q4. Is the external TCGA-GEO reference technically comparable to this cohort?

**Not established, and this is the largest un-addressed threat to the external
arm.** The panel is 2,015 normal and tumour-adjacent samples across 25 tissue
types (Downs, Thursby & Cope 2023). Unverified here: preprocessing and
normalisation pipeline, array version, batch structure, and the definition of
"tumour-adjacent". The paper itself concedes that "tissue-specific methylation
patterns will drive the threshold values at some CpG sites" and that a
multi-tissue reference is a deliberately conservative compromise.

**Partial reassurance, not a substitute:** if the panel were badly miscalibrated
for BRCA, whole cohorts would flag together at affected CpGs. They do not — at
most 11 of 53 samples flag at any chr22 CpG, no CpG reaches 50%, and 0% of
Normal flags sit in CpGs where half the cohort is flagged. That rules out gross
miscalibration; it does not establish comparability.

## Q5. Are all external-reference percentages using the same denominator?

**No — and the supervisor is right that this needed stating.** CpGs absent from
`tcga.rda` return all-NA, so `ext` is scored on fewer cells than `bio`/`self` at
the same sites. Coverage is **systematically uneven across states**, which is
the part that matters:

| state | Normal: CpGs absent | Tumour |
|---|---|---|
| `H` | 1.05% | 0.46% |
| `L` | 2.94% | 2.34% |
| `HM` | 3.30% | 3.38% |
| `R` | 8.91% | 8.27% |
| `LM` | **14.39%** | **11.47%** |
| `M` | **16.29%** | **24.82%** |

`Results/Task26_DenominatorAudit_*.csv`

Every reported percentage is computed over evaluable cells and so is internally
correct. But a raw **count** comparison across methods at `LM` or `M` sites
compares different denominators, and the uneven coverage is itself a mild
confound for any between-state comparison. Now stated in a preamble to
`results.md` and carried in the tables.

## Q6. Does "contrary flag" mean biologically implausible, or merely opposite to a coarse label?

**Merely opposite to a coarse label** — the stronger reading was not intended
and the wording invited it. `methy.state` is a six-level summary of a
continuous quantity; within-state heterogeneity is real and expected, and a
genuinely hypomethylated sample at an `L`-state CpG is biologically possible.

What makes the contrary flags *suspicious* is not their direction but their
**magnitude**: 97% of `H`-site flags move beta by less than 0.10. A contrary
flag with a large beta difference would be interesting; almost none are.

## Q7. What analysis excludes biological heterogeneity as the cause?

**None — and the claim has been revised.** "The cause is threshold geometry,
not biology" overstated the evidence. `results.md` now reads: *strongly
consistent with threshold geometry; biological and technical contributors have
not been excluded.*

What the evidence does establish:

1. Flag density per unit of available beta space is 44–171× higher in the
   state-constrained direction — a geometric property of the threshold that
   holds whatever the biology.
2. The decision boundary at `H` sites sits 0.003 apart in beta, inside array
   technical error, and 21% of flags do not survive noise of sd = 0.01. **A
   biological effect would not be destroyed by measurement noise this small.**
   This is the strongest single argument, and it is an argument about
   reproducibility, not about cause.
3. An independent statistic (Tukey 3×IQR) shows state dependence too, in its
   own direction, because the IQR is itself state-dependent.

A design that *would* separate the two: technical replicates, or an independent
platform (WGBS/EPIC) on the same samples. Neither is available here.

## Q8. Is N37's cluster a focal epimutation before probe QC?

**No, and the headline now matches the caveat.** Revised to *"a candidate focal
epimutation pattern, pending probe masking and technical-quality checks."*

The concern is specific: chr22:16.6 Mb (hg19) is pericentromeric on the
acrocentric short arm, repetitive and poorly mapped. A seven-probe run of
elevated beta there is also exactly what cross-hybridisation looks like.
Standard masking (Chen 2013; Zhou, Laird & Shen 2017) has not been applied.

What survives QC regardless: the per-sample summary counted one contiguous run
seven times. That is an arithmetic statement about the summary, not a
biological claim.

## Q9. Why 0.10 for the floor? Is it justified independently?

**It was chosen after seeing these data — a fair hit.** A sweep now replaces
the single point (`Results/Task26_FloorSensitivity_Normal.csv`):

| floor | flag rate | flags retained | stability | mag H÷R |
|---|---|---|---|---|
| 0 (current) | 2.67% | 100% | 0.328 | 0.117 |
| 0.05 | 1.69% | 63.3% | 0.748 | 0.243 |
| 0.10 | 1.03% | 38.6% | 0.794 | 0.492 |
| 0.15 | 0.54% | 20.2% | 0.820 | 0.526 |
| 0.20 | 0.31% | 11.6% | 0.868 | — |

Stability rises **monotonically**, so the curve cannot select a value — pushed
far enough the floor flags almost nothing and scores nearly perfectly. At 0.20
no `H`-state flag survives at all.

**Conclusion: the cut-off must be justified biologically**, as the smallest beta
difference worth calling real on the 450k platform, not read off this table.
`epimutacions` uses 0.15 and the 0.10-to-0.15 gap is small (0.794 vs 0.820), so
0.10 is not privileged. `results.md` §9 says this explicitly.

## Q10. Should ext+floor be recalibrated so its advantage is not a lower-rate artifact?

**Yes — done, and the conclusion survives.** Every competing method was re-tuned
*down* to the floor's exact rate (`Results/Task26_RateMatchedStability_*.csv`):

| method | Normal @ 1.031% | Tumour @ 6.272% |
|---|---|---|
| **ext + 0.10 floor** | **0.794** | **0.907** |
| tukey.iqr | 0.474 | 0.654 |
| mad.beta | 0.424 | 0.680 |
| beta.fit | 0.221 | 0.537 |
| mad.mvalue | 0.191 | 0.429 |

At identical sensitivity the floor still leads by a wide margin in both tissues.
The advantage is a property of the method, not of its flag rate. This was a real
weakness in the first draft and is now tested rather than caveated.

## Q11. Is Cohen's κ informative under extreme class imbalance?

**Only partly, and it should not stand alone.** With >95% zeros, κ is sensitive
to the marginal flag rates and can move for reasons unrelated to agreement.

`results.md` §6 already reports `agreement | flagged` (of cells flagged by at
least one method, the fraction where both agreed) alongside κ, and Task 20 also
writes Jaccard and the raw contingency counts to CSV. Both metrics point the
same way — κ 0.21 / agreement 12.2% in Normal, κ 0.11 / 6.8% in Tumour — which
is why the conclusion stands. The wording now says the methods *operationalise
outliers differently* rather than that they are "close to independent".

**Not yet done:** direction-specific agreement (do they agree on `+1` and `−1`
separately?) and precision/recall against a designated reference method. Both
are cheap; added to `plan.md`.

## Q12. Are chr22 conclusions representative genome-wide?

**Partly answered — and the magnitude result now is genome-wide.** Repeated over
all 380,355 CpGs (`Results/Task26_GenomeWideMagnitude.csv`):

| state | tissue | flags | median \|Δβ\| | under 0.10 | chr22 said |
|---|---|---|---|---|---|
| `H` | Normal | 78,259 | 0.031 | 96.97% | 97.7% |
| `L` | Normal | 174,567 | 0.042 | 81.53% | 82.2% |
| `R` | Normal | 38,131 | 0.269 | 3.68% | 3.6% |
| `H` | Tumour | 68,114 | 0.028 | 96.16% | 97.4% |
| `R` | Tumour | 1,608,650 | 0.213 | 19.92% | 17.6% |

Reproduces to within about a percentage point.

**Still chr22-only:** the threshold-geometry table and the noise-stability
experiment (§7). This limitation now appears beside those results, not only in
Next Steps.

## Q13. Were normal and tumour methylation states assigned independently?

**Yes, and this is now stated where the state tables are, not buried.**
`methy.state` is annotated per tissue, so the same cgID can sit in different
states in the two. On chr22, 2,856 of 6,809 sites change state, and 1,835 that
are not `R` in Normal become `R` in Tumour (`R` goes 308 → 2,095).

Consequence a reader must hold: **the Normal and Tumour rows of a state table
are not the same CpGs.** The site *set* is held constant; the state *labels*
are not. Comparing "the `H` row" across tissues compares different probes.

## Q14. What are `outliers.coef2` and `outliers.coef3`?

**Settled empirically.** They are Tukey outlier counts at two coefficients —
per CpG, the number of the 53 samples falling outside
`[Q1 − k·IQR, Q3 + k·IQR]`:

| column | k | agreement with recomputation |
|---|---|---|
| `outliers.coef2` | **2** | 100.0% identical |
| `outliers.coef3` | **3** | 100.0% identical |

Checked on 20,000 randomly sampled CpGs in the Normal cohort; k = 1.5 matches
neither. This is the same family as the Tukey comparison in §7, which uses
k = 3 — so `outliers.coef3` *is* that comparison, already computed and shipped
in the source files. Now defined once in `results.md`.

## Q15. Can every figure/table carry script, inputs, output, denominator, parameters?

**Partly, and honestly: not yet uniformly.** Every table in `results.md` names
its output CSV, and every CSV is produced by exactly one named script whose
`.Rout` transcript records inputs and parameters. Tissue and reference mode are
in the table captions; window and denominators are stated in §1's preamble and
in the per-table `n.evaluated` columns.

What is missing is a *uniform* provenance block under each object. That is a
formatting pass over ~20 tables and 5 figures, and is the one review item not
completed. Added to `plan.md` as a discrete task.

---

## Statements revised in `results.md`

| was | now |
|---|---|
| "the cause is threshold geometry, not biology" | "strongly consistent with threshold geometry; biological and technical contributors have not been excluded" |
| "the flags are statistically real and biologically empty" | "many flags carry small absolute beta differences, raising the concern that statistical flagging is not tracking biologically meaningful effect sizes" |
| "That is a focal epimutation." | "This is a candidate focal epimutation pattern, pending probe masking and technical-quality checks" |
| "The two methods are close to independent." | "Agreement is low once the shared-zero class is accounted for (κ ≈ 0.21 / 0.11), indicating the methods operationalise outliers differently" |
| "The absolute-difference floor is the fix." | "An absolute-difference floor is the strongest candidate tested so far, but the cut-off is not yet justified" — with the sensitivity sweep and rate-matched comparison |

---

# Part 2 — Readability audit

Standard applied: **a student who has taken one statistics course** — comfortable
with median, quantile, SD, correlation, a p-value — but who has not met Cohen's
κ, Jaccard, order statistics, or heteroscedasticity.

## What already works

- The `bio`/`self`/`ext` table appears before any result, so the three methods
  are defined once.
- "Answers in brief" gives the conclusion before the evidence.
- The degeneracy is shown rather than asserted — every `self` row reads
  `1.89 / 96.23 / 1.89`, and a reader can verify `1/53 = 1.89%` themselves.
- Concrete numbers throughout instead of significance language.

## What will lose that reader, in order of severity

**1. Four undefined technical terms.** `κ`, `Jaccard`, "order statistic",
"heteroscedastic". Each is load-bearing.
*Fix applied:* κ and Jaccard are now glossed at first use in §6; a plain-language
gloss of the order-statistic argument sits in the §-preamble. "Heteroscedastic"
should be replaced by "the spread of beta values is much smaller near 0 and 1
than in the middle" — **not yet done**, flagged below.

**2. The unit of analysis shifts without warning** — CpG, cell, sample burden,
event. §2 says "39 flags", §3 says "23 distinct CpGs across 20 samples", §4
says "11 flags = 3 events". A reader can easily think these contradict.
*Fix applied:* a unit-of-analysis table now sits before §1.

**3. `1.89%` looks like a result.** It is `1/53`. Without that, a reader sees a
suspiciously stable number and wonders what is special about 1.89.
*Fix applied:* stated inline as `1/53` and `51/53`.

**4. "Contrary flag" sounds like "wrong flag".** Addressed under Q6; §2 now
distinguishes the coarse label from biological implausibility.

**5. The tumour `H` row has n = 1 CpG.** "26.42%" is 14 cells at a single
probe. A reader will over-read it.
*Fixed:* the row is marked ⚠ with a note that a percentage over one CpG
estimates nothing, and readers are pointed to the genome-wide table instead.

**6. Zone-asymmetry arithmetic is compressed.** "30,363 flags per unit of hyper
zone against 177" asked the reader to reconstruct a division they were never
shown.
*Fixed:* the four-line calculation is now printed in a code block, from 1,283
flags through to the 171× ratio, with the note that an unbiased rule would give
a ratio near 1.

**7. "Heteroscedastic" carried the argument while going undefined.**
*Fixed:* replaced with "beta is trapped between 0 and 1, so the spread of values
is much smaller near 0 and near 1 than in the middle", with the technical term
kept in parentheses for anyone who wants to look it up. The companion idea —
that a percentile rule asks "who is in the tail?" and never "by how much?" — is
now stated in the same paragraph.

**8. Figure 3 needed a reading instruction.** Its coloured bands are threshold
*zones*, not data, which is not guessable.
*Fixed:* a three-line caption now says so and tells the reader which two bands
to compare.

## Verdict

All eight issues are now fixed in `results.md`. The remaining barrier is
structural rather than lexical: §7 and §9 assume the reader has accepted the
argument of §3, so they do not stand alone. A reader who starts at §9 will not
know why a "floor" is desirable. That is addressed by the "Answers in brief"
table and the §8-first pointer, both of which route a hurried reader through
the right order, but anyone extracting a single section should take §3 with it.

One item from the review is **not** done: a uniform provenance block (script,
inputs, output, denominator, parameters) under every table and figure. Every
object names its source CSV and every CSV has exactly one producing script with
a `.Rout` transcript, so the information is recoverable — but it is not uniform.
Tracked in `plan.md`.

---

# Part 3 — Code and methods review (2026-08-27, status-updated)

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

> **STATUS: FIXED 2026-08-28.** All three scripts corrected and guarded;
> corrected tables regenerated in Task 23. The impact was larger than estimated
> below — 65.0% of tumour *flag mass* in the external arm, not 40.3%. See
> `results.md` §8.

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

> **STATUS: SUPERSEDED 2026-08-28.** Task 20 reports Cohen's kappa and
> non-zero conditional agreement. See Part 1 Q11 for the caveat that kappa
> alone is unreliable under this much class imbalance.

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
