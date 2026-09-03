# Literature and methods review

Written for the 2026-09-03 task list. Three questions were asked:

1. What data is the `ext` reference actually made of, and how was it used here?
2. Does `OutlierMeth` identify outlying **samples**, or only outlying cells?
3. Can the parameters be adjusted to make it less sensitive?

Then two packages to read against ours: **Borealis** and **epimutacions**.

Everything below that is marked *verified* was read in the cited source.
Anything not marked was not checked and should not be quoted.

---

## 1. The `ext` reference

### What is in it

`OutlierMeth` ships four reference panels. Their composition is given in
Downs, Thursby & Cope 2023 (*Epigenetics* 18(1):2213874, `PMC10208159`) —
*verified*:

| panel | samples | tissue types | what they are |
|---|---|---|---|
| `tcga` | **747** | 21 | TCGA normal tissue, centralised QC and processing |
| `geo` | 1,268 | 15 | 29 GEO series, normal / tumour-adjacent; **preprocessing varies by study** |
| `all` | 2,015 | 25 | `tcga` + `geo` combined |
| `blood` | 656 | 1 | healthy whole blood, GSE40279, ages 19–101 |

All are Illumina HumanMethylation450. `blood` is kept separate by the authors
because whole blood has "very distinctive patterns of DNA methylation".

The thresholds are built by `referenceMeth()`, which is nothing but eight calls
to R's `quantile()` per probe — the 0.99 / 0.999 / 0.9999 / 0.99999 points and
their lower mirrors — stored as columns `P0.01 … P0.00001` and `N0.01 …
N0.00001`. No model, no shrinkage, no covariates.

### Correction to this repository's earlier statements

**Every doc in this repo previously described the external reference as
"2,015 samples across 25 tissue types". That is the `all` panel. This project
uses `tcga`, which is 747 samples across 21 tissue types.**

The code is unambiguous — `Task12`, `Task13.FullFlagMatrix` and `Task14` all do

```r
load(".../OutlierMeth/data/tcga.rda")
flag <- flagMeth(beta, reference = tcga, p = 0.01)
```

so the reference is `tcga`, not `all`. Corrected in `README.md`, `results.md`,
`plan.md` and `REVIEW.md` on 2026-09-03. The consequences are worked through in
Task 31 Part D and they are not cosmetic — at n = 747 three of the four `p`
levels collapse onto the reference maximum.

### How it was used here, exactly

- `flagMeth(beta, reference = tcga, p = 0.01)`, on all 380,355 CpGs, for
  Normal and Tumour separately (`Task13.FullFlagMatrix.Aug20.2026.R`).
- `tcga` covers 370,201 of our 380,355 CpGs. `flagMeth` initialises its output
  to `NA` and only fills `intersect(rownames(beta), rownames(reference))`, so
  the ~10,154 unmatched CpGs stay `NA` — which is why every denominator in this
  project is "evaluable cells", not "all cells". Task 12 got this wrong once
  (`plan.md` §8) and Task 26 Part A audits it.
- The result is a CpG × sample matrix of `-1 / 0 / +1`. Nothing else.

### Is the panel comparable to this cohort?

Still open, and it is the largest un-addressed threat to the external arm
(`REVIEW.md` Q4). Two things are now known:

- Gross miscalibration is ruled out — no CpG flags more than 11 of 53 samples,
  and no CpG flags ≥ 50% (Task 22, `plan.md` §1).
- But the panel pools 21 tissue types, and 29 GEO series with heterogeneous
  preprocessing enter the `geo` and `all` panels. A per-probe percentile over
  21 tissues is dominated by **between-tissue** variation, so it is not a
  breast-normal reference in any meaningful sense. That inflates the spread
  and should make the panel *less* sensitive, not more — yet the observed
  chr22 Normal flag rate is 2.67%, about 2.7× the 1% the `p = 0.01` setting
  nominally implies.

---

## 2. Does `OutlierMeth` identify outlying samples?

**No.** The package exports exactly four functions
(`Packages/OutlierMeth/NAMESPACE`):

| function | returns |
|---|---|
| `referenceMeth(beta)` | the 8-column per-probe threshold table |
| `flagMeth(beta, reference, p)` | CpG × sample matrix of `-1 / 0 / +1` |
| `deltMeth(beta, reference, p)` | the same cells as `beta − threshold` |
| `relMeth(beta, reference, p)` | the same as `(beta − P)/(1 − P)`, or `(N − beta)/N` |

There is no sample-level call, no burden statistic, no QC function and no
multiple-testing control. The unit is the **cell**, not the sample.

The paper does compute per-sample burden — "for each sample, the proportion of
CpG sites falling beyond the thresholds" — but as an analysis step, not as
package output (*verified*). Every per-sample table in this project (Tasks 13,
14, 15, 19, 21) is a column sum we compute ourselves, which is why
`plan.md` §5 can change the unit from flags to events without touching the
package.

That matters for the supervisor's question. If "outlying sample" is a wanted
output, it has to be defined here, with its own null — and `plan.md` §6 is the
list of confounders that must be excluded first.

### The two magnitude functions nobody is using

`deltMeth` returns `beta − threshold` for every flagged cell. **A magnitude
floor is one comparison away from something the package already computes**, and
this is the cheapest possible fix.

`relMeth` divides that by the remaining head-room, `1 − P`. At an `H` site the
head-room is 0.037, so relMeth *multiplies* a trivial shift by ~27. Measured on
chr22 Normal (Task 31 Part C): at `H` sites median `deltMeth` = 0.0022 while
median `relMeth` = 0.0746, an inflation of **34×**; at `HM` sites, 11×. A user
reaching for `relMeth` to get an interpretable effect size would get the
opposite — the flags whose direction the state made inevitable come back
looking like the largest effects in the dataset.

---

## 3. Can the parameters be adjusted?

Two knobs: `reference` (4 panels) and `p` (4 levels). Task 31 answers this in
full; the short version:

**The `p` argument has two usable settings on this panel, not four.**
`referenceMeth()` uses R's default type-7 quantile, which places level `1 − p`
at position `h = (n − 1)(1 − p) + 1` in the sorted reference. Level `p` is a
real tail estimate only when `n ≳ 1/p + 2`. For `tcga` at n = 747:

| `p` | position in the sorted 747 | rank from the top | expected n above |
|---|---|---|---|
| 0.01 | 739.54 | 8th | 7.5 |
| 0.001 | 746.25 | 1st–2nd | 0.75 |
| 0.0001 | 746.93 | 1st–2nd | 0.075 |
| 0.00001 | 746.99 | 1st–2nd | 0.007 |

The three strictest levels all land between the two most extreme reference
samples. `p = 0.0001` and `p = 0.00001` differ by 6% of one inter-sample gap.
No packaged panel is large enough to fix this — `p = 0.0001` needs n ≥ 10,002
and the largest panel is 2,015.

**This is the same order-statistic degeneracy the README documents for the
n = 53 self-reference arm, and it reaches the external arm too.**

**And tightening `p` would not help even if it worked.** Because thresholds are
monotone in `p`, the flag set is nested: any stricter setting keeps a prefix of
each CpG's current flags, ordered by beta. Task 31 Part B enumerates that whole
family from the `p = 0.01` matrix. At the strictest reachable setting — one
flag per CpG per direction — chr22 Normal `H`-site `+1` flags are still
**100% below |Δβ| = 0.10, with 0 of 536 surviving a 0.10 floor**. At
`p = 0.01` the same cell is 1,127 flags, also 100% trivial, also 0 surviving.

The ceiling is structural: at `H` sites the median cohort beta is 0.934, so
there is 0.066 of beta space above it, and only **0.24%** of all `H` cells sit
≥ 0.10 from their site median. No threshold can find what is not there.

So the answer to "adjust the parameters to get results" is: **no. The `p`
argument changes how many flags you get, never how big they are.** What is
needed is a term the package does not have — a magnitude requirement. Task 32
builds and scores five.

---

## 4. Borealis

*Bioconductor. Vignette read 2026-09-03 — verified. The methods paper was not
read; nothing here is attributed to it.*

**What it is.** Outlier methylation detection for **bisulfite sequencing
(BS-seq) read counts**, aimed at rare-disease diagnosis. Fits a
**beta-binomial** per CpG: `x` methylated reads out of `n` total, with a cohort
mean `mu` and a dispersion `theta`. Each sample's counts are then tested
against that fitted distribution; `runBorealis()` is the pipeline entry point,
`plotCpGsite()` the per-site view. Output carries `pVal`, an `effSize`, and
`isHypo` (yes / no / `NA` when the change is too small to call). Multiple
testing is left to the user via `p.adjust()`.

**Why it matters to us.** Borealis is the model-based version of what this
project keeps concluding must happen:

1. **It models dispersion explicitly.** The beta-binomial has a variance term
   per CpG, so a fixed absolute shift is not equally significant everywhere.
   This is the principled answer to the heteroscedasticity that makes a
   percentile threshold meaningless at `H` and `L` sites.
2. **Read depth carries its own uncertainty.** `n` is in the model. A 450k beta
   value has no `n` — its uncertainty has to be supplied from outside, which is
   exactly what our measured-noise estimate (Task 32 Part A) is doing by hand.
3. **It has an explicit "too small to call" state.** `isHypo = NA`. That is a
   magnitude floor, built into the output rather than bolted on.

**Why it cannot be used directly.** It needs Bismark-style read counts. We have
450k beta values, which are intensity ratios with no underlying count. A
beta-binomial cannot be fitted to them without inventing an `n`. The
transferable idea is the *shape* — per-CpG dispersion plus an explicit
undecidable class — not the code.

The closest thing already benchmarked here is `beta.fit` in Task 24: a per-CpG
Beta(a, b) fitted by moments with a two-sided tail probability. It scored
0.390 / 0.578 stability — better than M-values, worse than the external
reference with a floor. A proper beta-binomial with shrunk dispersion would be
a fairer test of the idea than that moment fit.

---

## 5. epimutacions

*Bioconductor 3.23. Vignette read 2026-09-03 — verified against the release
vignette; the reference manual PDF could not be text-extracted on this machine.*

**What it is.** Six outlier methods for 450k/EPIC beta values, run either as
**case vs a reference panel** (`epimutations()`) or **leave-one-out within a
cohort** (`epimutations_one_leave_out()`).

| method | parameters and defaults |
|---|---|
| `manova` | `pvalue_cutoff = 0.05` |
| `mlm` | `pvalue_cutoff = 0.05` |
| `iForest` | `outlier_score_cutoff = 0.7`, `ntrees = 100` |
| `mahdist` | `nsamp = "deterministic"` |
| `quantile` | `window_sz = 1000`, `qsup = 0.995`, `qinf = 0.005`, **`offset_abs = 0.15`** |
| `beta` | `pvalue_cutoff = 1e-06`, **`diff_threshold = 0.1`** |

An **epimutation** is "a consecutive window of a minimum of 3 outlier CpGs with
a maximum distance of 1 kb between them".

**The three things worth taking.**

1. **Both threshold-based methods ship a mandatory magnitude floor.** The
   quantile method requires `offset_abs = 0.15` on top of the 0.005/0.995
   quantiles; the beta method requires `diff_threshold = 0.1` on top of
   `p < 1e-6`. Neither will report a statistically extreme but biologically
   trivial cell. This is independent precedent for the fix this project
   arrived at from the data, and it is why `plan.md` §4 uses 0.15 as the
   reference offset.

   Note the structure: `quantile` is *the same rule as `flagMeth`* — an
   empirical quantile of a reference panel — with one extra term. The
   difference between the two packages at that method is precisely
   `offset_abs`.

2. **The region definition is the unit.** ≥3 outlier CpGs within 1 kb, not a
   single cell. `plan.md` §5 already adopted this (Task 21 computes
   `epimutations` on that definition) after N37's eleven flags turned out to be
   three events. Our `events.1kb` and `max.run.cpgs` are directly comparable.

3. **Leave-one-out is a first-class mode, not a workaround.** `bio.loo` in
   Tasks 18–21 is the same idea. Worth noting that `epimutations()` proper uses
   an external *case-control* design, so the honest comparison to our `ext` arm
   is `epimutations()`, and to our `self` arm is `epimutations_one_leave_out()`.

**Where our contribution sits.** epimutacions supplies a floor but fixes it at
one constant across all CpGs. Task 32 asks whether a **state-specific** floor
does better, and the calibrated per-state floors (Normal chr22, 96th percentile
of each state's own deviation distribution) come out as:

| state | `L` | `LM` | `M` | `HM` | `H` | `R` |
|---|---|---|---|---|---|---|
| floor | 0.040 | 0.154 | 0.162 | 0.148 | 0.049 | 0.327 |

A single constant of 0.15 is close to right for `LM`, `M` and `HM`, three times
too strict for `L` and `H`, and half of what `R` needs. So epimutacions' 0.15
is a reasonable global compromise, and a per-state floor is a refinement of it
rather than a replacement — which is what Task 32's scores show (0.781 vs 0.793
stability in Normal; the constant floor still wins on the headline metric).

---

## 6. What is claimed, and on what evidence

| claim | source | status |
|---|---|---|
| `tcga` panel = 747 samples, 21 tissue types | Downs et al. 2023 | verified in the paper |
| `all` = 2,015 / 25; `geo` = 1,268 / 15; `blood` = 656 | Downs et al. 2023 | verified in the paper |
| `OutlierMeth` has no sample-level output | `NAMESPACE`, `R/*.R` | verified in the source |
| Paper cautions against outlier flagging for gradual/continuous processes | Downs et al. 2023 | verified in the paper |
| No delta-beta filter in `flagMeth` | `R/flagMeth.R` | verified in the source |
| `p` levels collapse at n = 747 | arithmetic on `quantile` type 7 | computed, Task 31 Part D |
| relMeth inflates `H`-site flags 34× | Task 31 Part C | computed on our data |
| Borealis: beta-binomial on BS-seq counts, `isHypo = NA` class | vignette | verified in the vignette |
| Borealis dispersion shrinkage details | — | **not checked** — the vignette does not describe it |
| epimutacions defaults incl. `offset_abs = 0.15`, `diff_threshold = 0.1` | vignette | verified in the vignette |
| epimutation = ≥3 CpGs within 1 kb | vignette | verified in the vignette |
| Where epimutacions applies `delta_beta` internally | — | **not checked** |

---

## References

- Downs BM, Thursby S-J, Cope L. *Detecting aberrant DNA methylation in
  Illumina DNA methylation arrays: a toolbox and recommendations for its use.*
  Epigenetics 18(1):2213874, 2023.
  [PMC10208159](https://pmc.ncbi.nlm.nih.gov/articles/PMC10208159/) ·
  [publisher](https://www.tandfonline.com/doi/full/10.1080/15592294.2023.2213874)
- Borealis (Bioconductor).
  [package page](https://bioconductor.org/packages/release/bioc/html/borealis.html) ·
  [vignette](https://bioconductor.org/packages/release/bioc/vignettes/borealis/inst/doc/borealis.html)
- epimutacions (Bioconductor).
  [vignette](https://www.bioconductor.org/packages/release/bioc/vignettes/epimutacions/inst/doc/epimutacions.html) ·
  [reference manual](https://bioconductor.posit.co/packages/3.23/bioc/manuals/epimutacions/man/epimutacions.pdf)
- Du P, Zhang X, Huang C-C, et al. *Comparison of Beta-value and M-value
  methods for quantifying methylation levels by microarray analysis.*
  BMC Bioinformatics 11:587, 2010.
- McCartney DL, et al. *DNA methylation outlier burden, health, and ageing in
  Generation Scotland and the Lothian Birth Cohorts.*
  [PMC7098133](https://pmc.ncbi.nlm.nih.gov/articles/PMC7098133/)
