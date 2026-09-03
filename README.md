# Outlier DNA Methylation in TCGA-BRCA

Does an "outlier methylation" flag report real biology, or does it mostly
report the methylation state the CpG already sits in?

That is the whole project. A CpG that is unmethylated in every normal tissue
cannot become much less methylated; calling one of 53 samples "hypomethylated"
there may be a tautology dressed as a finding. This repository tests that,
using the `OutlierMeth` toolbox on 53 matched normal/tumour breast pairs.

**Current answer: the bias is real and large, it is a bias of *magnitude*
rather than *rate*, and it is not a bug in any one package — it follows from
thresholding bounded, heteroscedastic beta values on a scale-free statistic.
At the methylation extremes the decision boundary falls inside the array's own
noise, and the resulting flags do not survive a perturbation smaller than the
array's technical error.**

**Over half the external reference's normal-tissue flags on chr22 run *against*
the state the CpG sits in, and no setting of the package's parameters removes
them: the ceiling on effect size at an `H` site is 0.066 of beta space, so even
the strictest reachable threshold leaves 100% of `H`-site hyper flags below
|Δβ| = 0.10. The fix has to be a term the package does not have.** See
[Findings](#findings).

---

## Data

| | |
|---|---|
| Cohort | TCGA-BRCA, 53 patients with matched normal + tumour (the "53 Alive" set); a separate 32-patient "Dead" cohort is used in Experiment 6 |
| Platform | Illumina HumanMethylation450 (beta values) |
| Sites | 380,355 CpGs passing upstream filtering |
| Annotation | each CpG carries a precomputed `methy.state`: `L`, `LM`, `M`, `HM`, `H`, `R` |
| chr22 | 6,809 CpGs — `L` 2,211, `LM` 981, `M` 762, `HM` 1,648, `H` 899, `R` 308 (Normal annotation) |

Source matrices live outside this repository, on `leap2.txstate.edu` under
`/home/s_s355/research3.data/TCGA/filter.BRCA.UCEC.Aug17.2022/`. They are
individual-level genomic data and are deliberately not committed; see
[CLASSIFICATION.md](CLASSIFICATION.md) for the public/private split and the
reasoning behind each call.

`methy.state` is annotated **per tissue**, so the same cgID can sit in
different states in Normal and Tumour. On chr22, 2,856 of 6,809 sites change
state between the two, and 1,835 that are not `R` in Normal become `R` in
Tumour (`R` goes 308 → 2,095). Every table below states which annotation it
used.

---

## Method

`OutlierMeth` flags each CpG × sample cell as `-1` (hypo), `0`, or `+1`
(hyper) by comparing the beta value against per-CpG thresholds. Thresholds
come from one of two places, and the distinction is the spine of the project:

- **External reference** — the packaged `tcga.rda` panel: **747 TCGA normal
  samples across 21 tissue types**, independent of this cohort. This is the
  design the package was published for.

  *(Corrected 2026-09-03. Earlier drafts described this as "2,015 samples
  across 25 tissue types" — that is `all`, the combined TCGA+GEO panel. Every
  script here loads `tcga.rda` and calls `flagMeth(..., reference = tcga)`, so
  the reference is the 747-sample one. The distinction matters: see Finding 6.)*
- **Self reference** — `referenceMeth()` run on the same 53 samples being
  tested. This is *not* how the package is meant to be used, and at n = 53 it
  is degenerate (see below).

A third, project-specific method is compared against both:

- **`bio`** — a state-aware rule proposed in the Aug 2026 meetings:
  `L`/`LM` sites flag `+1` only for beta above the 99th percentile;
  `H`/`HM` sites flag `-1` only for beta below the 1st percentile;
  `M`/`R` sites get no rule. `bio.loo` is the same rule with the threshold
  for sample *i* computed from the other 52.

### The self-reference degeneracy

`quantile(x, 0.99)` over n = 53 values is an order statistic of those same
values. R's type-7 rule places it at position `52 × 0.99 + 1 = 52.48`, strictly
between the 52nd and 53rd order statistics, so **exactly one sample can exceed
it — the maximum**. The 0.01 quantile admits exactly one below — the minimum.

Every CpG therefore receives exactly one `-1` and one `+1`, always, regardless
of its state, its variance, or anything biological. Every significance level
from 0.01 to 0.00001 lands in the same inter-order-statistic gap, so `p` is
inert. Raising *n* raises the constant (n = 201 → 4 flags/CpG; n = 1000 → 20)
but never creates variation.

Consequences, worked through in [REVIEW.md](REVIEW.md) Finding 1: the
self-referential by-state tables were incapable of showing a difference before
the data was read, and Task 11 ("does removing samples 14–17 change the mean
flag count?") could not have moved a statistic pinned at 2.

**The `bio` and `self` columns are retained in every table as a negative
control — they show what "no signal" looks like. Read the `ext` column.**

---

## Findings

### 1. Contrary flags exist, and they are common

The Aug 2026 email asked: at `H`/`HM` sites, is there ever a `+1`? At `L`/`LM`
sites, is there ever a `-1`? Against the external reference, yes — routinely.

chr22:50,459,312–50,532,711 (the 100-CpG window of minimum span, Task 20):

| state | tissue | contrary flag | count | % of evaluated cells |
|---|---|---|---|---|
| `L` | Normal | `-1` | 39 / 2,120 | 1.84% |
| `H` | Normal | `+1` | 16 / 848 | 1.89% |
| `L` | Tumour | `-1` | 52 / 2,014 | 2.58% |
| `HM` | Tumour | `+1` | 47 / 1,060 | 4.43% |
| `H` | Tumour | `+1` | 40 / 636 | 6.29% |

The `bio` rule produces zero by construction — it only ever assigns `+1` at
`L`/`LM` and `-1` at `H`/`HM`. That is the point of comparing them.

### 2. The bias is in magnitude, not rate

The contrary-flag *rate* is unremarkable. The contrary-flag *size* is not.
Measuring `|beta − cohort median|` for every chr22 cell the external reference
flagged (Task 21, Part E):

| state | tissue | flags | median \|Δβ\| | % under 0.05 | % under 0.10 |
|---|---|---|---|---|---|
| `H`  | Normal | 1,283 | **0.031** | 78.3% | **97.7%** |
| `L`  | Normal | 3,774 | **0.035** | 52.7% | **82.2%** |
| `HM` | Normal | 2,019 | 0.117 | 7.1% | 39.3% |
| `LM` | Normal | 1,495 | 0.117 | 27.0% | 44.2% |
| `M`  | Normal | 513 | 0.144 | 0.0% | 19.3% |
| `R`  | Normal | 564 | **0.267** | 0.0% | **3.6%** |
| `H`  | Tumour | 1,473 | 0.031 | 82.6% | 97.4% |
| `L`  | Tumour | 3,822 | 0.026 | 59.5% | 87.3% |
| `R`  | Tumour | 19,170 | **0.214** | 6.7% | **17.6%** |

**This is the project's central result.** Nearly every outlier flag at an `H`
or `L` site corresponds to a beta shift too small to mean anything — 98% of
`H`-site flags move less than 0.10. At `R` sites the same method produces
flags with a median shift of 0.21–0.27.

The mechanism is known and not specific to `OutlierMeth`: beta values are
severely heteroscedastic, their standard deviation compressed below 0.2 and
above 0.8 (Du et al. 2010). A percentile threshold is scale-free and does not
care. `L` and `H` sites live exactly in the compressed zone, so a statistically
extreme value there is a biologically trivial one.

An absolute-difference floor fixes it. `epimutacions` already ships one
(`offset_abs = 0.15`). At that floor, 99% of the `H`-site and `L`-site flags in
this dataset disappear, while 68% (Tumour) to 83% (Normal) of `R`-site flags
survive. The floor removes almost exactly the flags whose direction the state
made inevitable, and Task 24 measures what that buys: stability rises from 0.33
to 0.90 (Normal) and 0.65 to 0.95 (Tumour).

**M-values are not the fix, despite being the textbook correction for beta
heteroscedasticity.** Benchmarked against six alternatives at matched flag
rates, a median ± k·MAD rule on M-values is the *least* stable method tested in
both tissues (0.308 / 0.474), below the same rule on raw beta (0.558 / 0.628).
The logit stabilises variance but amplifies noise at the extremes — a 0.01
shift near β = 0.97 is a large move in M-space — and the extremes are exactly
where the problem lives. See `results.md` §9.

### 2b. The mechanism: the threshold sits inside the array's noise band

Task 22 reconstructs the external thresholds from the flag matrix itself —
`flagMeth` sets `+1` iff `beta > P`, so `max(beta | unflagged) ≤ P < min(beta |
flagged)` brackets `P` per CpG. Across 13,618 CpG × tissue combinations there
were **zero** cases where a flagged sample's beta sat inside the unflagged
range, which also independently confirms the Task 13 flag/beta join is in
register.

chr22, Normal, median per state:

| state | `P` threshold | hyper zone (1−P) | hypo zone (N−0) | asymmetry | gap between the two samples straddling P |
|---|---|---|---|---|---|
| `L`  | 0.100 | 0.900 | **0.021** | 44× | 0.024 |
| `LM` | 0.226 | 0.774 | 0.049 | 16× | 0.057 |
| `M`  | 0.682 | 0.318 | 0.360 | 0.9× | 0.028 |
| `HM` | 0.910 | **0.090** | 0.677 | 0.13× | 0.013 |
| `H`  | 0.963 | **0.037** | 0.886 | 0.04× | **0.003** |
| `R`  | 0.762 | 0.238 | 0.269 | 0.9× | 0.052 |

At an `H` site the hyper-outlier zone is 0.037 of the beta scale and the
hypo-outlier zone is 0.886 — a 24-fold asymmetry — yet **87.8% of the observed
flags are `+1`**. Normalising by the room available: 30,363 flags per unit of
hyper zone against 177 per unit of hypo zone, a **171× higher flag density in
the direction the state already constrains**. At `L` sites the same ratio is
44×.

And the decision boundary is not resolvable. At `H` sites the two samples
straddling `P` are **0.003 apart**; 450k technical replicate SD is roughly
0.01–0.03.

**The flag is therefore not reproducible.** Re-flagging chr22 after adding
N(0, sd) noise and clamping to [0,1], 20 replicates per level (Task 22 Part C,
Normal, 9,648 baseline flags):

| noise sd | flags lost | flags gained | % baseline lost | churn / baseline |
|---|---|---|---|---|
| 0.005 | 1,487 | 6,833 | 15.4% | 0.86 |
| **0.010** | **2,062** | **13,467** | **21.4%** | **1.61** |
| 0.020 | 2,693 | 23,187 | 27.9% | 2.68 |
| 0.050 | 3,485 | 43,726 | 36.1% | 4.91 |

At sd = 0.01 — inside the array's own technical error — more calls change than
there were flags to begin with. Instability tracks the compressed states
exactly: `H` 34.2%, `L` 25.4%, `HM` 15.5%, against `R` 7.5%.

### 2c. It is not an `OutlierMeth` bug

Tukey's 3×IQR rule, computed within this cohort, shares nothing with the
external panel but the data. It is also the rule behind the `outliers.coef2/3`
columns already sitting unused in the source files. It has a state bias too —
a *different* one (Task 22 Part E, Normal):

| state | median cohort IQR | percentile flag rate | Tukey flag rate | Jaccard |
|---|---|---|---|---|
| `L`  | **0.008** | 3.30% | 3.39% | 0.294 |
| `LM` | 0.052 | 3.38% | 2.61% | 0.242 |
| `M`  | 0.081 | 1.41% | 0.29% | 0.065 |
| `HM` | 0.068 | 2.37% | 0.40% | 0.076 |
| `H`  | 0.023 | 2.69% | 0.77% | 0.050 |
| `R`  | 0.157 | 3.91% | 0.71% | 0.101 |

Two rules, two different state profiles, almost no overlap. The reason is in
the second column: the IQR that Tukey thresholds on is *itself* a function of
the state (0.008 at `L`, 0.157 at `R`), so Tukey inherits the same pathology in
its own direction.

**The defect is not specific to `OutlierMeth`.** It follows from applying any
scale-free dispersion threshold to bounded, heteroscedastic beta values. Which
is why the fix — an absolute-difference floor, or M-values — applies to the
whole class.

### 2d. Over half the flags run against the state

The sharpest single statement of the problem. Scoring every method on the share
of its flags that are *contrary* — `−1` at an `L`/`LM` site, `+1` at an `H`/`HM`
site — on chr22 (Task 32):

| method | contrary share, Normal | Tumour |
|---|---|---|
| external reference, as published | **52.9%** | 27.0% |
| + a 0.10 magnitude floor | 21.0% | 6.5% |
| + a per-state floor | 14.9% | 21.8% |
| state-pooled deviation rule | **0%** by construction | **0%** |

### 6. Neither `OutlierMeth` parameter can fix it

The package has two knobs, `reference` (4 panels) and `p` (4 levels), and
Task 31 closes off both.

`referenceMeth()` builds every threshold with R's default type-7 `quantile()`,
which places level `1 − p` at position `(n − 1)(1 − p) + 1`. That is a real tail
estimate only when `n ≳ 1/p + 2`. On the 747-sample `tcga` panel:

| `p` | rank from the top of 747 | expected n above |
|---|---|---|
| 0.01 | 8th | 7.5 |
| 0.001 | 1st–2nd | 0.75 |
| 0.0001 | 1st–2nd | 0.075 |
| 0.00001 | 1st–2nd | 0.007 |

The three strictest levels all land between the two most extreme reference
samples. `p = 0.0001` would need n ≥ 10,002 and the largest packaged panel is
2,015. **The same order-statistic degeneracy that invalidates the n = 53
self-reference arm also caps the external arm at two usable settings.**

And tightening `p` would not help. Thresholds are monotone in `p`, so the flag
set is nested and the whole reachable family can be enumerated from the
`p = 0.01` matrix. At the strictest reachable setting, chr22 Normal `H`-site
`+1` flags are still **100% below |Δβ| = 0.10, with 0 of 536 surviving a 0.10
floor** — identical to `p = 0.01`. The ceiling is structural: only 0.24% of all
`H` cells sit ≥ 0.10 from their site median.

**`p` changes how many flags you get. It never changes how big they are.**

The package's own `relMeth`, which divides the shift by the remaining head-room
`1 − P`, points the wrong way: at `H` sites it turns a median `deltMeth` of
0.0022 into a median `relMeth` of 0.0746, an inflation of **34×**. `deltMeth`,
by contrast, is one comparison away from the magnitude floor this project
recommends.

### 7. The contrary flags concentrate on sites, not on samples

The 2026-09-03 email asked whether the contrary flags in the 100-CpG window
come from a couple of CpGs or from all of them, and from a couple of samples or
from any of them. Both halves have answers, and they differ (Task 29).

**Sites: concentrated.** Half the sites of a state carry every contrary flag it
has, and the top two hold 57–86% of them. One probe, `cg15668074`, carries 19
of the 29 tumour `LM` `−1` flags on its own.

**Samples: spread.** 31 of 53 normal samples and 37 of 53 tumour samples carry
at least one, with the top three holding 27% against a uniform expectation of
6%. A mild gradient, not a two-sample story.

And the concentrating sites are not biologically special — they are the sites
where the external threshold happens to sit closest to this cohort's own
distribution, at a median of **1.0–1.9 cohort MADs** from the median. Spearman
ρ between a site's contrary-flag count and that distance is negative in all six
testable state × tissue cells (−0.95 to −0.12).

### 8. A state-pooled threshold fixes the degeneracy but not the magnitude

The email's "bio-stat" proposal — one threshold per *state* from the pooled
`n_sites × 53` cells, rather than one per CpG from 53 — was tested two ways
(Task 30). Pooling the raw beta flags whole sites; pooling the **deviation**
from each site's own median asks the intended question.

`bio.stat.dev` is the best rule tested that needs no external panel at all:
stability 0.779 / 0.833 against `mad.beta`'s 0.429 / 0.678 at the same flag
rate, a state-rate ratio of 1.0, and zero contrary flags by construction.

It does not fix the magnitude problem, because the 96th percentile of the `L`
and `H` states' own deviation distributions is only 0.037 and 0.046. Read as
data-calibrated floors, those percentiles answer the "is 0.05 the right number"
question directly (chr22 Normal):

| state | `L` | `LM` | `M` | `HM` | `H` | `R` |
|---|---|---|---|---|---|---|
| calibrated floor | 0.037 | 0.142 | 0.144 / 0.179 | 0.125 | 0.046 | 0.306 / 0.341 |

**0.05 is about right at `L` and `H` and three to seven times too lenient
everywhere else.** Benchmarked head to head, a flat 0.10 still beats both 0.05
and the per-state version on stability (`results.md` §13).

### 8b. Genome-wide: 134 contrary flags out of 154,810 clear a 0.10 floor

Everything in Findings 2b, 2c and 6 was measured on chr22. Task 33 repeats it
on all 380,355 CpGs and it replicates closely.

| tissue | cell | flags | median \|Δβ\| | under 0.10 | **survive a 0.10 floor** |
|---|---|---|---|---|---|
| Normal | `H` `+1` | 68,878 | 0.029 | 99.95% | **36** |
| Normal | `L` `−1` | 85,932 | 0.015 | 99.89% | **98** |
| Normal | `R` `+1` | 22,221 | **0.294** | 1.27% | **21,938** |
| Tumour | `H` `+1` | 60,218 | 0.026 | 99.97% | **17** |
| Tumour | `L` `−1` | 81,948 | 0.013 | 99.95% | **39** |

**0.087% of the external reference's contrary flags at `H` and `L` sites carry
a beta shift of 0.10 or more, against 98.7% at `R` sites.** The threshold
geometry replicates too — 169.8× flag density in the constrained direction at
`H` sites genome-wide against 171× on chr22, and 88.0% of `H` flags `+1`
against 87.8%.

And the plainest version of the whole argument: **at `L` sites, where every
sample is unmethylated, the external reference calls almost exactly as many
samples abnormally low as abnormally high — 85,932 against 88,635 — inside a
hypo zone 0.022 of the beta scale wide.**

### 8c. The per-sample burden is not heteroscedasticity

Applying a 0.10 floor genome-wide and re-ranking (Task 33 Part E), normal
tissue keeps the **same top six in the same order** — N15, N17, N14, N48, N31,
N34 — at Spearman 0.910, retaining 41% of flags overall. One sample separates
out: **N27 falls from rank 7 to rank 28**, keeping only 10.7% of its flags.

So the leaders' burden survives a magnitude requirement. That does not make it
biological — batch, detection-p, purity and cell composition are all still
unexcluded ([plan.md](plan.md) §6) — but it does rule out the cheapest
explanation, and it raises the stakes on excluding the rest. Tumour reshuffles
more (T17 2 → 8, T16 4 → 12), so quote the floored ranking there.

### 9. The state-aware rule and the external reference do not agree

Cohen's kappa and non-zero conditional agreement over the tight window
(Task 20, Part D) — raw percent agreement is meaningless on matrices that are
>95% zeros, which is why [REVIEW.md](REVIEW.md) Finding 4 rejected Task 17's
original metric:

| tissue | pair | cells flagged by either | agreed | agreement | kappa |
|---|---|---|---|---|---|
| Normal | `bio` vs `ext` | 172 | 21 | 12.2% | 0.206 |
| Normal | `self` vs `ext` | 235 | 63 | 26.8% | 0.411 |
| Tumour | `bio` vs `ext` | 415 | 28 | 6.8% | 0.109 |
| Tumour | `self` vs `ext` | 455 | 99 | 21.8% | 0.333 |

The two methods are close to independent. They are not two views of one
outlier set; they are two different definitions of outlier.

### 10. N37 tops the Task 17 summary for three separable reasons

N37 carried 11 of 71 biological flags in the original 100-CpG window, ~8× the
uniform expectation of 1.34. Tasks 19 and 21 decompose that:

| scale | method | flags | expected | rank of 53 |
|---|---|---|---|---|
| index100 | `bio` | 11 | 1.34 | **1** |
| index100 | `ext` | 1 | 2.83 | 30 |
| tight100 | `bio` | 1 | 1.66 | 17 |
| chr22 (6,809 CpGs) | `bio` | 24 | 108.3 | **41** |
| chr22 | `ext` | 114 | 182.0 | 25 |

- **The window is not representative.** Chromosome-wide, N37 carries less than
  a quarter of its expected burden and ranks 41st. Spearman correlation between
  index100 rank and chr22 rank is 0.344.
- **The external reference does not see it.** One flag, rank 30.
- **Ten flags are one event.** N37's flags collapse to 3 runs at a 1 kb gap;
  the largest holds 7 CpGs.

But the run itself is real, and this is a correction to the Task 19 write-up:
across `cg01836687`–`cg23572163` (chr22:16,601,097–16,602,837, 1,740 bp, 7
consecutive `LM` sites) N37 sits at beta 0.276–0.553 while the cohort median
runs 0.015–0.089. Its median margin **over the second-highest sample** is 0.195
within that run, 0.152 across all eleven flags. That is a focal
hypermethylation event of substantial magnitude, not a sample that narrowly won
seven maxima.

So the per-sample summary is right for the wrong reason: it counts one
epimutation seven times. The correct unit is the event, and by that unit N37
has one.

*Caveat:* chr22:16.6 Mb (hg19) lies in the pericentromeric region of the
acrocentric short arm, which is repetitive and poorly mapped. Cross-reactive
and polymorphic probe masking has not yet been applied — see
[plan.md](plan.md) step 2.

### 11. Which samples actually carry burden

N15, N17, N14 and N48 lead at every scale and in both tissues — N48 was
originally missed and is confirmed 4th of 53 both on chr22 and genome-wide
(Finding 8c). Their chr22
external-reference burden tracks their genome-wide self-referential burden at
Spearman 0.761 (Normal) and 0.843 (Tumour), so this is a global sample
property, not a chr22 story — plausibly global hypomethylation, tumour purity,
or a technical batch effect. None of those have been tested yet.

---

## Task inventory

| | What it does | Status |
|---|---|---|
| Exp 1–3 | Load BRCA matrices; chrX vs chr22; inspect 5 high-outlier CpGs | done |
| Exp 5–8 / Tasks 6, 8 | Dead cohort; combine cohorts; filter low-variance H/L sites; remove outlier samples | done |
| Task 9 | Boxplots + histograms of `outliers.coef2/3` by state | done |
| Task 10 | Flag × state contingency table | **fixed + guarded**; corrected output in Task 23 |
| Task 11 | Effect of removing samples 14–17 | **retired** — measures a constant |
| Task 12 | Self vs external reference, genome-wide means | **rerun needed** — denominator mismatch, ~2.7% low |
| Task 13 | Full 4-matrix flag pipeline (self + ext × Normal + Tumour) | done — the data spine |
| Task 14 | External-reference flag bias by state and sample | **fixed + guarded**; the `R` rows carried 65% of tumour flag mass — Task 23 |
| Task 15 | Self-referential bias by state and sample, + bar plots | Part A degenerate; Part B (per-sample) sound |
| Task 17 | chr22 pilot: 100 CpGs, state-aware rule vs external | superseded by 18/20; agreement metric invalid |
| Task 18 | chr22 state summary, contrary flags, drill-down | done |
| Task 19 | Per-sample debug, N37, event collapsing | done |
| **Task 20** | Tight window, method concordance (kappa), annotated contrary flags | **new** |
| **Task 21** | Burden across three scales, contrary-flag drivers, magnitude analysis | **new** |
| **Task 22** | Threshold geometry, flag stability under noise, percentile vs Tukey | **new** |
| **Task 23** | Tasks 10/13A/14 regenerated with the `R` state restored | **new** |
| **Task 24** | Seven-method benchmark at matched flag rates | **new** |
| **Task 25** | Figures | done |
| **Task 26** | Supervisor-review response: denominators, floor sensitivity, rate-matched stability, genome-wide magnitude | done |
| **Task 27** | The complete 100-CpG sheet, site by site | done |
| **Task 28** | Window selection: minimise span or maximum gap | done |
| **Task 29** | Contrary-flag decomposition over sites and samples; per-site threshold geometry | **new** |
| **Task 30** | State-pooled ("bio-stat") thresholds, absolute and deviation variants | **new** |
| **Task 31** | Parameter headroom: effect-size ceiling, `p` envelope, `deltMeth` vs `relMeth`, `p`-level resolution | **new** |
| **Task 32** | Measured noise from adjacent probes; five candidate fixes scored | **new** |
| **Task 33** | Ceiling, `p` envelope, geometry, stability and burden-under-floor, genome-wide | **new** |
| **Task 34** | Figures 8–12 | **new** |
| **Task 35** | Direction-split agreement and precision/recall | **new** |

[plan.md](plan.md) has the ordered next steps.

---

## Running it

Every script from Task 18 onward resolves its own paths — it uses the leap2
locations when they exist and a local working copy otherwise, so the same file
runs unmodified in both places. Earlier scripts hardcode absolute HPC paths.

```sh
R CMD BATCH --no-save --no-restore \
  Scripts/Task20.Chr22.TightWindow.Concordance.Aug28.2026.R \
  Scripts/Task20.Chr22.TightWindow.Concordance.Aug28.2026.Rout
```

`Slurm/*.sbatch` holds the batch equivalents. `data.table` is used when
available and falls back to base R readers. Tasks 20 and 21 run in about 7 and
30 seconds respectively on a laptop, given the Task 13 flag matrices.

Note that Tasks 20 and 21 need `Task13_*_{self,ext}ref_flags.csv` in
`Results/` — those are private (~62 MB each) and are not in the repository.
Task 20's threshold-lookup columns additionally want `tcga.rda`; it degrades
to `NA` for those columns when the panel is absent.

### Dependency

`OutlierMeth` (Downs & Cope, JHMI, GPL-2) is vendored under `Packages/` at
commit `408b110` from `github.com/bdowns4/OutlierMeth`. It should be recorded
as a dependency and removed from the tree — see [CLASSIFICATION.md](CLASSIFICATION.md).

---

## Layout

```
Data/          third-party reference tracks (CpG islands, SNVs, IMR90) — not ours to redistribute
Logs/          install logs, slurm stdout
Packages/      vendored OutlierMeth + a scratch experiment
Results/       CSV/XLSX outputs and figures; the large per-sample matrices are gitignored
Scripts/       all analysis code, each with its matching .Rout transcript
Slurm/         batch submission scripts
CLASSIFICATION.md  public/private split and the reasoning
REVIEW.md          code and methods review; six findings, ordered by impact
LITERATURE.md      the OutlierMeth reference panels, its parameters, Borealis, epimutacions
plan.md            what to do next
```

---

## References

- Downs BM, Thursby S-J, Cope L. *Detecting aberrant DNA methylation in
  Illumina DNA methylation arrays: a toolbox and recommendations for its use.*
  Epigenetics 18(1):2213874, 2023.
  [PMC10208159](https://pmc.ncbi.nlm.nih.gov/articles/PMC10208159) —
  the `OutlierMeth` paper. Note its own conclusion that outlier flagging suits
  binary tumour/normal discrimination and is *unwise* for continuous processes.
- Du P, Zhang X, Huang C-C, et al. *Comparison of Beta-value and M-value
  methods for quantifying methylation levels by microarray analysis.*
  BMC Bioinformatics 11:587, 2010.
  [link](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-587) —
  beta-value heteroscedasticity; the mechanism behind Finding 2.
- Barbosa M, et al. *epimutacions* (Bioconductor).
  [vignette](https://www.bioconductor.org/packages/release/bioc/vignettes/epimutacions/inst/doc/epimutacions.html) —
  six outlier methods (MANOVA, MLM, isolation forest, robust Mahalanobis,
  quantile, beta), an epimutation defined as ≥3 outlier CpGs within 1 kb, a
  leave-one-out mode, and an absolute-difference floor.
- McCartney DL, et al. *DNA methylation outlier burden, health, and ageing in
  Generation Scotland and the Lothian Birth Cohorts.*
  [PMC7098133](https://pmc.ncbi.nlm.nih.gov/articles/PMC7098133/) —
  3×IQR outlier definition, per-sample burden, and the finding that burden is
  strongly confounded by cell composition and technical factors.
- Borealis (Bioconductor) — beta-binomial outlier detection for bisulfite
  sequencing read counts, with an explicit "too small to call" class.
  [package](https://bioconductor.org/packages/release/bioc/html/borealis.html) —
  the model-based version of what this project argues for; needs read counts,
  which 450k beta values do not have. See [LITERATURE.md](LITERATURE.md) §4.
- CoMeBack — co-methylated region construction for 450k data.
  [Bioinformatics 36:2675](https://academic.oup.com/bioinformatics/article/36/9/2675/5716323)
