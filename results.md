# Results

Findings for the week of 2026-08-24, against the four questions set in the
pre-meeting email. Every number here is reproducible from a script in
`Scripts/` and a CSV in `Results/`; the source table for each is named beneath
it.

Cohort: TCGA-BRCA, 53 matched normal/tumour pairs, 380,355 CpGs on the
Illumina 450k array. Chromosome 22 carries 6,809 of them.

**Three flagging methods are compared throughout:**

| | how the threshold is set | |
|---|---|---|
| `bio` | state-aware rule: `L`/`LM` flag `+1` above the 99th percentile of the 53 values; `H`/`HM` flag `−1` below the 1st; `M`/`R` no rule | the rule proposed in the meetings |
| `self` | `referenceMeth()` on the same 53 samples | how the package was run in Tasks 10–15 |
| `ext` | the packaged TCGA-GEO panel: 2,015 independent normal and tumour-adjacent samples | how the package was **designed** to be run |

> **Read the `ext` column.** At n = 53, `quantile(x, 0.99)` sits between the
> 52nd and 53rd order statistic, so exactly one sample can exceed it. `bio` and
> `self` therefore assign a fixed number of flags to every eligible CpG no
> matter what the data says, at every significance level from 0.01 to 0.00001.
> They are kept in every table below as a **negative control** — they show what
> "no signal" looks like — not as results.

---

## Answers in brief

| # | Question | Answer |
|---|---|---|
| 1 | Summarise chr22 by state | Done, both tissues, all 6 states — §1 |
| 2 | Any `+1` at `H`/`HM`? Any `−1` at `L`/`LM`? | **Yes, routinely.** Up to 6.3% of cells at tumour `H` sites — §2 |
| 3 | If non-zero, which sample, which site, why? | Three samples drive a third of them; the cause is threshold geometry, not biology — §3 |
| 4 | Debug the sample summary, N37 on top | N37 is **41st of 53** chromosome-wide. Its window rank was an artifact — but it does carry one real focal event — §4 |
| 5 | Explore genuinely adjacent CpGs | The original 100 sites span 5.17 Mb. A true 73 kb window changes the picture — §5 |
| 6 | Compare the state-aware rule with the external reference | They barely agree: κ = 0.21 (normal), 0.11 (tumour) — §6 |

Beyond the four questions, §7–§9 report why the flags behave this way, a
bug fix that changes a headline number, and a benchmark of seven alternatives.

---

## 1. Chromosome 22 summarised by state

The email's example — *"H, 11 CG sites"* — identifies the window exactly: the
first 100 chr22 CpGs by position, where `H` has 11 sites. Reproduced here in
full.

### Normal — 100 CpGs, chr22:11,915,060–17,085,407

| state | CpGs | method | −1 | 0 | +1 | % −1 | % 0 | % +1 |
|---|---|---|---|---|---|---|---|---|
| **L** | 8 | bio | 0 | 416 | 8 | 0.00 | 98.11 | 1.89 |
| | | self | 8 | 408 | 8 | 1.89 | 96.23 | 1.89 |
| | | **ext** | **8** | 359 | 4 | **2.16** | 96.77 | 1.08 |
| **LM** | 14 | bio | 0 | 728 | 14 | 0.00 | 98.11 | 1.89 |
| | | self | 14 | 714 | 14 | 1.89 | 96.23 | 1.89 |
| | | **ext** | **21** | 718 | 3 | **2.83** | 96.77 | 0.40 |
| **M** | 24 | bio | 0 | 1272 | 0 | 0.00 | 100.00 | 0.00 |
| | | self | 24 | 1224 | 24 | 1.89 | 96.23 | 1.89 |
| | | ext | 12 | 1241 | 19 | 0.94 | 97.56 | 1.49 |
| **HM** | 38 | bio | 38 | 1976 | 0 | 1.89 | 98.11 | 0.00 |
| | | self | 38 | 1938 | 38 | 1.89 | 96.23 | 1.89 |
| | | **ext** | 16 | 1952 | **46** | 0.79 | 96.92 | **2.28** |
| **H** | 11 | bio | 11 | 572 | 0 | 1.89 | 98.11 | 0.00 |
| | | self | 11 | 561 | 11 | 1.89 | 96.23 | 1.89 |
| | | **ext** | 2 | 571 | **10** | 0.34 | 97.94 | **1.72** |
| **R** | 5 | bio | 0 | 265 | 0 | 0.00 | 100.00 | 0.00 |
| | | self | 5 | 255 | 5 | 1.89 | 96.23 | 1.89 |
| | | ext | 5 | 256 | 4 | 1.89 | 96.60 | 1.51 |

`Results/Task20_Chr22_index100_StateSummary_Normal.csv`

Note every `self` row reads `1.89 / 96.23 / 1.89` — that is `1/53` and `51/53`,
to the digit, in all six states. This is the degeneracy, visible.

### Tumour — same 100 sites

`methy.state` is annotated per tissue, so these sites sit in different states in
tumour: `R` grows from 5 to 52 of the 100.

| state | CpGs | method | −1 | 0 | +1 | % −1 | % 0 | % +1 |
|---|---|---|---|---|---|---|---|---|
| **L** | 7 | ext | 4 | 307 | 7 | 1.26 | 96.54 | 2.20 |
| **LM** | 8 | ext | **29** | 362 | 33 | **6.84** | 85.38 | 7.78 |
| **M** | 13 | ext | 178 | 473 | 38 | 25.84 | 68.65 | 5.52 |
| **HM** | 19 | ext | 204 | 760 | **43** | 20.26 | 75.47 | **4.27** |
| **H** | 1 | ext | 14 | 39 | 0 | 26.42 | 73.59 | 0.00 |
| **R** | 52 | ext | 551 | 2010 | 195 | 19.99 | 72.93 | 7.08 |

`Results/Task20_Chr22_index100_StateSummary_Tumor.csv` (bio/self rows omitted —
identical constants)

### Whole chromosome, and the whole genome

The 100-site window is small. The same summary over all 380,355 CpGs
(`Results/Task23_StateFlagTable_extref_*.csv`) is in §8, because producing it
required fixing a bug first.

![Flag rate by state](Results/Fig1_FlagRateByState.png)

---

## 2. Contrary flags: `+1` at `H`/`HM`, `−1` at `L`/`LM`

A `+1` at a site that is already high, or a `−1` at a site that is already low,
is a flag pointing where there is no room. The email asked whether any exist.

**Yes — in every state, in both tissues, at every scale tested.**

| tissue | state | contrary flag | count | of cells | rate |
|---|---|---|---|---|---|
| Normal | `L` | `−1` | 39 | 2,120 | 1.84% |
| Normal | `LM` | `−1` | 5 | 530 | 0.94% |
| Normal | `HM` | `+1` | 9 | 1,113 | 0.81% |
| Normal | `H` | `+1` | 16 | 848 | 1.89% |
| Tumour | `L` | `−1` | 52 | 2,014 | 2.58% |
| Tumour | `LM` | `−1` | 33 | 583 | **5.66%** |
| Tumour | `HM` | `+1` | 47 | 1,060 | **4.43%** |
| Tumour | `H` | `+1` | 40 | 636 | **6.29%** |

External reference, 100-CpG tight window (§5).
`Results/Task20_Chr22_tight100_ContraryCounts_*.csv`

The `bio` rule returns **zero** in every row — it can only assign `+1` at
`L`/`LM` and `−1` at `H`/`HM`, so a contrary flag is impossible by
construction. That is worth stating explicitly, because it means the `bio`
column in this table is not evidence of anything.

---

## 3. Which sample, which site, and why

### Which samples

Three samples supply about a third of all state-contrary external flags on
chr22, and they are the same three in both tissues.

| sample | contrary `L` | `LM` | `HM` | `H` | **total** | all ext flags | % contrary |
|---|---|---|---|---|---|---|---|
| N17 | 101 | 18 | 195 | 96 | **410** | 1,250 | 32.8% |
| N15 | 84 | 13 | 192 | 92 | **381** | 1,297 | 29.4% |
| N14 | 118 | 18 | 149 | 84 | **369** | 1,039 | 35.5% |
| N48 | 153 | 93 | 72 | 38 | 356 | 496 | 71.8% |
| N31 | 48 | 37 | 147 | 69 | 301 | 405 | 74.3% |
| N27 | 98 | 15 | 60 | 91 | 264 | 278 | **95.0%** |

`Results/Task21_Chr22_ContraryDrivers_Normal.csv`

N15, N17 and N14 also lead the burden ranking at every scale in both tissues,
and their chr22 external burden tracks their genome-wide burden at Spearman
0.761 (normal) / 0.843 (tumour). **This is a whole-sample property, not a
chr22 story** — global hypomethylation, tumour purity, or a technical batch
effect. None has been tested yet; that is the first item in `plan.md`.

Note the last column. For the lower-burden samples, *almost every* external
flag is state-contrary — 95% for N27. That is the opposite of what a
biological reading would predict.

### Which sites

The 39 contrary `L`-state flags in the tight window fall on 23 distinct CpGs
across 20 samples. Of the 69 contrary external flags there:

- **0** lie in a CpG island
- **0** have an SNV within 10 bp
- **2** are also Tukey 3×IQR outliers within this cohort
- **48** are isolated — no flagged neighbour within 1 kb

`Results/Task20_Chr22_tight100_ContraryDetail_Normal.csv` (private — carries
per-sample beta at named positions)

So the usual technical explanations do not apply. Which leaves the question:

### Why

Because the difference being flagged is negligible. Measuring
`|β − cohort median|` for every chr22 cell the external reference flagged:

| state | tissue | flags | median \|Δβ\| | under 0.05 | under 0.10 |
|---|---|---|---|---|---|
| `H` | Normal | 1,283 | **0.031** | 78.3% | **97.7%** |
| `L` | Normal | 3,774 | **0.035** | 52.7% | **82.2%** |
| `LM` | Normal | 1,495 | 0.117 | 27.0% | 44.2% |
| `HM` | Normal | 2,019 | 0.117 | 7.1% | 39.3% |
| `M` | Normal | 513 | 0.144 | 0.0% | 19.3% |
| `R` | Normal | 564 | **0.267** | 0.0% | **3.6%** |
| `H` | Tumour | 1,473 | 0.031 | 82.6% | 97.4% |
| `L` | Tumour | 3,822 | 0.026 | 59.5% | 87.3% |
| `R` | Tumour | 19,170 | 0.214 | 6.7% | 17.6% |

`Results/Task21_Chr22_ExtFlagMagnitude.csv`

**98% of the flags at `H` sites move beta by less than 0.10.** At `R` sites the
median flag moves it by 0.21–0.27. The flags are statistically real and
biologically empty, and exactly where the state predicts. §7 gives the
mechanism.

![Magnitude by state](Results/Fig2_FlagMagnitudeByState.png)

---

## 4. Debugging the sample summary — N37

N37 topped the Task 17 summary with 11 of 71 biological flags, ~8× the uniform
expectation of 1.34. Widening the denominator settles it:

| scale | CpGs | method | flags | expected | **rank of 53** |
|---|---|---|---|---|---|
| index100 | 100 | `bio` | 11 | 1.34 | **1** |
| index100 | 100 | `ext` | 1 | 2.83 | 30 |
| tight100 | 100 | `bio` | 1 | 1.66 | 17 |
| chr22 | 6,809 | `bio` | 24 | 108.3 | **41** |
| chr22 | 6,809 | `ext` | 114 | 182.0 | 25 |

`Results/Task21_Chr22_SampleBurden_Normal.csv`

Three independent reasons N37 came out on top, none of them "N37 is an outlier
sample":

1. **The window is not representative.** Chromosome-wide, N37 carries less than
   a quarter of its expected burden. Spearman between index100 rank and chr22
   rank is only 0.344.
2. **The external reference does not see it** — 1 flag, rank 30.
3. **Ten flags are one event.** They collapse to 3 runs at a 1 kb gap; the
   largest holds 7 CpGs.

### But the event itself is real

This corrects the Task 19 write-up, which implied N37 had merely won a series
of narrow maxima. It had not:

| cgID | position | state | N37 β | cohort median | 2nd highest | **margin** | gap |
|---|---|---|---|---|---|---|---|
| cg01836687 | 16,601,097 | LM | 0.553 | 0.052 | 0.141 | **0.413** | — |
| cg26822097 | 16,601,691 | LM | 0.402 | 0.039 | 0.341 | 0.061 | 594 |
| cg12431879 | 16,601,897 | LM | 0.442 | 0.089 | 0.289 | 0.152 | 206 |
| cg00332021 | 16,602,522 | LM | 0.435 | 0.048 | 0.213 | 0.222 | 625 |
| cg00816224 | 16,602,592 | LM | 0.438 | 0.046 | 0.243 | 0.195 | 70 |
| cg15554678 | 16,602,673 | LM | 0.356 | 0.015 | 0.111 | 0.245 | 81 |
| cg23572163 | 16,602,837 | LM | 0.276 | 0.065 | 0.211 | 0.064 | 164 |

`Results/Task21_Chr22_N37_Cluster_Normal.csv` (private)

Seven consecutive `LM` CpGs across 1,740 bp, N37 at β 0.28–0.55 against a
cohort median of 0.015–0.089, median margin 0.195 over the next-highest sample.
**That is a focal epimutation.** The per-sample summary was right for the wrong
reason: it counted one event seven times.

The correct unit is the event, not the flag — following the `epimutacions`
definition of ≥3 outlier CpGs within 1 kb. By that unit N37 has one.

![N37 cluster](Results/Fig5_N37_LocalCluster.png)

> **Caveat before this is called biology.** chr22:16.6 Mb (hg19) is
> pericentromeric on the acrocentric short arm — repetitive and poorly mapped.
> A seven-probe run there is also what a cross-hybridisation artifact looks
> like. Standard 450k probe masking (Chen 2013 / Zhou 2017) has **not** been
> applied to this project yet. That is `plan.md` step 2.

---

## 5. Genuinely adjacent CpG sites

The email asked for sites *close to each other*. The 100 CpGs used in Tasks
17–19 are consecutive **by index**, not by distance:

| window | span | median gap | largest gap | pairs within 1 kb |
|---|---|---|---|---|
| `index100` (Tasks 17–19) | **5,170,347 bp** | 512 bp | **3,410,179 bp** | 61 / 99 |
| `tight100` (new) | **73,399 bp** | 187 bp | 13,823 bp | 80 / 99 |

`Results/Task20_Chr22_WindowComparison.csv`

`tight100` is the 100 consecutive chr22 CpGs with the smallest genomic span,
found by sliding a 100-probe window across all 6,809: chr22:50,459,312–50,532,711.
It also has a better state mix (`L` 40, `HM` 21, `H` 16, `LM` 11, `M` 9, `R` 3).

Both windows are analysed side by side in every Task 20 output, so the legacy
result stays reproducible and the difference is visible. The window choice
matters: N37 drops from rank 1 to rank 17 between them, and the Spearman
correlation of sample ranks against the chromosome-wide truth improves from
0.344 to 0.553.

---

## 6. The state-aware rule vs the external reference

Raw percent agreement is meaningless here — all four flag matrices are >95%
zeros, so any two of them "agree" ~96% of the time by sharing blanks. (This is
why the Task 17 agreement figure of 44–53 out of 53 was uninformative.) Two
honest metrics instead:

- **agreement | flagged** — of cells flagged by *at least one* method, the
  fraction where both gave the same value
- **Cohen's κ** — chance-corrected

| tissue | pair | flagged by either | agreed | agreement \| flagged | κ |
|---|---|---|---|---|---|
| Normal | `bio` vs `ext` | 172 | 21 | **12.2%** | **0.206** |
| Normal | `self` vs `ext` | 235 | 63 | 26.8% | 0.411 |
| Normal | `bio` vs `self` | 200 | 88 | 44.0% | 0.604 |
| Tumour | `bio` vs `ext` | 415 | 28 | **6.8%** | **0.109** |
| Tumour | `self` vs `ext` | 455 | 99 | 21.8% | 0.333 |

`Results/Task20_Chr22_tight100_Concordance_*.csv`

By state (Normal, `bio` vs `ext`): κ = 0.26 at `H`, 0.17 at `HM`, 0.17 at `L`,
0.35 at `LM`, and 0.00 at `R` — where `bio` has no rule at all.

**The two methods are close to independent.** They are not two views of one
outlier set; they are two different definitions of the word.

---

## 7. Why the flags behave this way

Task 22 reconstructs the external thresholds from the flag matrix itself:
`flagMeth` sets `+1` iff `β > P`, so `max(β | unflagged) ≤ P < min(β | flagged)`
brackets `P` per CpG.

> This doubles as a correctness check on the whole pipeline. Across 13,618
> CpG × tissue combinations there were **zero** cases where a flagged sample's
> beta sat inside the unflagged range — so the Task 13 flag/beta join is in
> register.

### The outlier zones are wildly asymmetric

| state | `P` | hyper zone (1−P) | hypo zone (N−0) | gap between the two samples straddling `P` |
|---|---|---|---|---|
| `L` | 0.100 | 0.900 | **0.021** | 0.024 |
| `LM` | 0.226 | 0.774 | 0.049 | 0.057 |
| `M` | 0.682 | 0.318 | 0.360 | 0.028 |
| `HM` | 0.910 | **0.090** | 0.677 | 0.013 |
| `H` | 0.963 | **0.037** | 0.886 | **0.003** |
| `R` | 0.762 | 0.238 | 0.269 | 0.052 |

Normal, chr22 medians. `Results/Task22_Chr22_ThresholdGeometry_Normal.csv`

At an `H` site the hyper-outlier zone is 0.037 of the beta scale and the hypo
zone is 0.886 — a 24-fold asymmetry — yet **87.8% of observed flags are `+1`**.
Per unit of room available that is 30,363 flags in the hyper direction against
177 in the hypo: a **171× higher flag density in the direction the state
already constrains**. At `L` sites, 44×.

![Threshold geometry](Results/Fig3_ThresholdGeometry.png)

### The decision boundary sits inside the array's noise

At `H` sites the two samples straddling the threshold are **0.003 apart**. 450k
technical replicate SD is roughly 0.01–0.03. Re-flagging after adding N(0, sd)
noise, 20 replicates per level (Normal, 9,648 baseline flags):

| noise sd | flags lost | flags gained | % baseline lost | changed calls ÷ baseline |
|---|---|---|---|---|
| 0.005 | 1,487 | 6,833 | 15.4% | 0.86 |
| **0.010** | **2,062** | **13,467** | **21.4%** | **1.61** |
| 0.020 | 2,693 | 23,187 | 27.9% | 2.68 |
| 0.050 | 3,485 | 43,726 | 36.1% | 4.91 |

`Results/Task22_Chr22_FlagStability_Normal.csv`

At sd = 0.01 — inside the array's own error — **more calls change than there
were flags to begin with.** Instability tracks the compressed states exactly:
`H` 34.2%, `L` 25.4%, `HM` 15.5%, against `R` 7.5%.

![Flag stability](Results/Fig4_FlagStability.png)

### It is not an `OutlierMeth` bug

Tukey's 3×IQR — a different statistic, and the one behind the `outliers.coef2/3`
columns already sitting unused in the source files — has a state bias too, a
*different* one:

| state | median cohort IQR | percentile rate | Tukey rate | Jaccard |
|---|---|---|---|---|
| `L` | **0.008** | 3.30% | 3.39% | 0.294 |
| `LM` | 0.052 | 3.38% | 2.61% | 0.242 |
| `M` | 0.081 | 1.41% | 0.29% | 0.065 |
| `HM` | 0.068 | 2.37% | 0.40% | 0.076 |
| `H` | 0.023 | 2.69% | 0.77% | 0.050 |
| `R` | 0.157 | 3.91% | 0.71% | 0.101 |

`Results/Task22_Chr22_PercentileVsTukey_Normal.csv`

The reason is the second column: the IQR that Tukey thresholds on is *itself* a
function of the state. **The defect follows from applying any scale-free
dispersion threshold to bounded, heteroscedastic beta values** — not from one
package.

### What was ruled out

- **Panel miscalibration.** If the BRCA cohort simply sat outside the
  2,015-sample panel's range, many samples per CpG would flag. At most 11 of 53
  do; zero CpGs reach 50%; 0% of normal flags sit in such CpGs. These are
  genuine per-sample calls.
- **A join or register error.** Zero bracket violations, as above.

---

## 8. Bug fix: the `R` state was being deleted

`Task10`, `Task13.FlagBias` and `Task14.FlagBias` all declared:

```r
state_order <- c("L", "LM", "M", "HM", "H", "Rc")   # the data has "R"
```

The subsequent `tab[ord, ]` kept only matched states, so **every `R`-state CpG
was silently dropped**. All three scripts are now fixed and carry a guard that
stops on any state present in the data but missing from `state_order`.

The bug was worse than the earlier review estimated — it was framed as "40% of
tumour CpGs", but by flag mass:

| reference | tissue | `R` CpGs | % of CpGs | `R` flags | % of **all flags** |
|---|---|---|---|---|---|
| self | Normal | 23,455 | 6.2% | 46,910 | 6.2% |
| self | Tumour | 153,479 | 40.4% | 306,958 | 40.4% |
| ext | Normal | 23,455 | 6.2% | 38,131 | 7.4% |
| **ext** | **Tumour** | **153,479** | **40.4%** | **1,608,650** | **65.0%** |

`Results/Task23_RStateImpact.csv`

Task 14's tumour external-reference table — the strongest result in the project
— was built on **35% of the flag mass**. Corrected, genome-wide:

| state | CpGs | % −1 | % 0 | % +1 | **% flagged** |
|---|---|---|---|---|---|
| `L` | 76,106 | 2.08 | 95.84 | 2.08 | 4.16 |
| `LM` | 40,824 | 3.41 | 89.66 | 6.93 | 10.34 |
| `M` | 5,436 | 8.61 | 82.63 | 8.76 | 17.37 |
| `HM` | 74,430 | 4.55 | 89.56 | 5.89 | 10.44 |
| `H` | 30,080 | 0.50 | 95.71 | 3.79 | 4.29 |
| **`R`** | **153,479** | **10.14** | **78.44** | **11.42** | **21.56** |

Tumour, external reference, all 380,355 CpGs.
`Results/Task23_StateFlagTable_extref_Tumor.csv`

The five original rows reproduce to the digit, which validates the
regeneration. The missing row turns out to have **the highest flag rate of any
state** — and, from §3, the largest effect sizes. The category that was deleted
is the one carrying the signal.

---

## 9. What to use instead — benchmark of seven methods

Seven outlier definitions on the same chr22 cells. **Every method is calibrated
by bisection to the external reference's own flag rate**, so stability and
state-dependence are compared at equal sensitivity — an uncalibrated comparison
of outlier methods is close to meaningless.

- `stability` — Jaccard between the flag set on clean data and on data
  perturbed by N(0, 0.01), the array's technical error. Self-derived methods
  are **fully recomputed** on the perturbed data, so threshold re-estimation is
  included. Higher is better.
- `state rate ratio` — max ÷ min flag rate across the six states; 1 = perfect
  state-independence.
- `mag H÷R` — median |Δβ| of flags at `H` sites ÷ the same at `R` sites; 1 =
  equally large effects regardless of state.
- `% trivial` — share of flags with |Δβ| < 0.10.

### Normal

| method | rate | **stability** | state rate ratio | mag H÷R | % trivial |
|---|---|---|---|---|---|
| **`ext` + 0.10 floor** | 1.09 | **0.899** | 62.0 | **0.492** | **0.0** |
| `tukey.iqr` | 2.82 | 0.583 | 4.53 | 0.135 | 59.5 |
| `mad.beta` | 2.82 | 0.558 | 5.54 | 0.128 | 60.0 |
| `loo.quantile` | 3.77 | 0.397 | 1.00 † | 0.145 | 58.5 |
| `beta.fit` | 2.82 | 0.390 | 1.75 | 0.151 | 45.4 |
| `ext.percentile` (current) | 2.82 | 0.328 | 2.78 | 0.117 | 61.5 |
| `mad.mvalue` | 2.82 | **0.308** | 3.22 | 0.140 | 52.6 |

### Tumour

| method | rate | **stability** | state rate ratio | mag H÷R | % trivial |
|---|---|---|---|---|---|
| **`ext` + 0.10 floor** | 6.61 | **0.953** | 128.3 | **0.463** | **0.0** |
| `ext.percentile` (current) | 10.78 | 0.648 | 4.65 | 0.144 | 38.7 |
| `mad.beta` | 10.78 | 0.628 | 1.96 | 0.117 | 48.6 |
| `tukey.iqr` | 10.78 | 0.600 | 1.37 | 0.115 | 50.5 |
| `beta.fit` | 10.78 | 0.578 | 1.71 | 0.117 | 40.9 |
| `loo.quantile` | 10.78 | 0.522 | 1.01 † | 0.116 | 52.0 |
| `mad.mvalue` | 10.78 | **0.474** | 1.24 | 0.114 | 49.1 |

`Results/Task24_MethodBenchmark_Scores.csv`

### Three conclusions, one of which reverses an earlier recommendation

**1. The absolute-difference floor is the fix.** Adding `|β − cohort median| ≥
0.10` on top of the existing external reference roughly triples stability
(0.33 → 0.90 normal, 0.65 → 0.95 tumour), quadruples the magnitude ratio
(0.12 → 0.49), and removes trivial flags by construction. It is one line of
code on top of what already runs. `epimutacions` ships the same idea as
`offset_abs = 0.15`.

Its state-rate ratio *rises* to 62–128, and that is the point rather than a
defect: it removes nearly all `H`/`L` flags because those flags had no
magnitude. Whether that is correct depends on accepting that sites with no room
genuinely have fewer real outliers — which §7 argues they do.

*Caveat:* the floor is a filter, not a calibrated method, so its rate (1.09%)
is below the others and some of its stability advantage is a lower-rate effect.

**2. M-values do not help — they make it worse.** This reverses the
recommendation in the earlier README and plan. `mad.mvalue` is the **least**
stable method tested in both tissues (0.308 / 0.474), below the same rule on
raw beta (0.558 / 0.628). The logit transform stabilises *variance* but
amplifies *noise* at the extremes: a 0.01 shift near β = 0.97 is a large move
in M-space. Since the extremes are exactly where the problem lives, the
standard textbook correction is the wrong tool here. Recorded as a negative
result.

**3. † `loo.quantile`'s perfect state-rate ratio is degeneracy, not virtue.**
It flags a fixed count per CpG by construction, so its rate cannot vary by
state. Read it alongside its stability (0.40 / 0.52), which is poor.

---

## Files

| | |
|---|---|
| `Scripts/Task20.*` | tight window, concordance (κ), annotated contrary flags |
| `Scripts/Task21.*` | burden across three scales, contrary drivers, magnitude |
| `Scripts/Task22.*` | threshold geometry, stability, percentile vs Tukey |
| `Scripts/Task23.*` | corrected state tables, `R` restored |
| `Scripts/Task24.*` | seven-method benchmark |
| `Scripts/Task25.*` | figures |
| `Results/Fig1–5*.png` | the five figures above |
| `Slurm/Task2*.sbatch` | batch equivalents |

Each script has a matching `.Rout` console transcript. Tables naming a sample
*and* a genomic position are individual-level and stay private — see
`CLASSIFICATION.md`.

`README.md` has the project overview and method; `plan.md` has the ordered next
steps; `REVIEW.md` has the code-level findings.

---

## Next, in order

1. **Probe masking** (Chen 2013 / Zhou 2017). Blocks any biological claim,
   including N37's event.
2. **Confirm §7 genome-wide.** The mechanism rests on chr22; the Task 13
   matrices already cover all 380,355 sites.
3. **Measure the noise instead of assuming it.** §7's stability result uses a
   literature range. Technical replicates — or adjacent co-methylated probe
   pairs — would replace it with this cohort's own value.
4. **Test what drives N15/N17/N14.** Batch, detection-p, cell composition,
   purity. Until excluded, per-sample burden cannot be interpreted.
