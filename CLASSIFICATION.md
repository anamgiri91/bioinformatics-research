# Public / Private Classification

Last reviewed: 2026-08-28, against the full `Experiments.Outlier.July31.2026`
tree pulled from `leap2.txstate.edu:/mmfs1/home/wln26/`.

The split is enforced by `.gitignore`, not by moving files. Everything stays
where the scripts expect it; git simply only ever sees the public half.

| | Files | Size |
|---|---|---|
| **Public** (committed to GitHub) | 115 | 1.7 MB |
| **Private** (local + leap2 only) | 9 paths | 586 MB |

The governing rule: **anything that resolves to an individual sample stays
private; anything that collapses across samples is public.**

---

## Private

### Individual-level methylation data

| Path | Size | Why |
|---|---|---|
| `Results/Task13_Normal_selfref_flags.csv` | 62 MB | 380,355 CpGs x 53 individuals |
| `Results/Task13_Normal_extref_flags.csv` | 61 MB | same |
| `Results/Task13_Tumor_selfref_flags.csv` | 62 MB | same |
| `Results/Task13_Tumor_extref_flags.csv` | 63 MB | same |
| `Results/Task13_FullFlagMatrix.xlsx` | 272 MB | all four of the above in one workbook |
| `Results/Task17_Chr22_BiologicalValidation_Normal_100CpGs.csv` | 15 KB | 100 CpGs x 53 individuals |
| `Results/Task18_Chr22_ContraryFlagDetail_*.csv` | 27 KB | names cgID + position + sample + that sample's beta |
| `Results/Task19_Chr22_TopSampleDetail_*.csv` | 41 KB | same |
| `Results/Task20_Chr22_*_ContraryDetail_*.csv` | 120 KB | same, plus cohort quartiles and MAD-z |
| `Results/Task21_Chr22_N37_Cluster_Normal.csv` | 2 KB | one named sample's betas across 11 identified sites |

Samples are de-identified (`N1..N53` / `T1..T53`) and carry no TCGA barcodes —
I grepped the whole tree and found none. That is necessary but not sufficient:
a per-individual matrix over 380k methylation sites is high-dimensional enough
to be treated as individual-level genomic data regardless of the label on the
column. These belong in a controlled repository (dbGaP / institutional store),
not a public git remote.

The Task 18-21 drill-down tables are small, but size is not the test.
Each row resolves a named sample to a named genomic position and reports its
beta value, which is the same class of disclosure as the flag matrices above
and is narrower only in count. Their aggregate siblings - `StateSummary`,
`ContraryCounts`, `Concordance`, `SampleBurden`, `ScaleComparison`,
`ContraryDrivers`, `ExtFlagMagnitude` - carry counts per sample with no
genomic coordinates and stay public, consistent with the treatment of the
Task 13/14 sample summaries below.

`Task13_FullFlagMatrix.xlsx` is independently disqualified — 272 MB exceeds
GitHub's 100 MB hard per-file limit and would reject the push outright.

### Reference threshold tables

| Path | Size | Why |
|---|---|---|
| `Results/experiment7_ref_normal.rds` | 24 MB | see below |
| `Results/experiment7_ref_tumor.rds` | 24 MB | see below |

These *look* like aggregates, which is exactly why they are worth calling out.
`referenceMeth()` stores per-CpG quantiles at 0.99 / 0.999 / 0.9999 / 0.99999.
With n = 53, every one of those upper quantiles interpolates between the 52nd
and 53rd order statistic, so the stored value is a near-exact readout of the
single most-methylated individual at that CpG — and the `N*` columns likewise
pin the least-methylated one. Publishing them publishes two individuals' beta
values at 380k sites. Private.

### Source data

`Data/` — `Mar24.r2.imr90.chr22.730883.txt` (18 MB), `SNV.all.txt` (3 MB),
`UCSC.CpGI.chr22.txt` (42 KB). Third-party reference data that is not ours to
redistribute, sourced from another lab's space on leap2
(`/home/s_s355/research3.data/`). UCSC CpG islands are publicly downloadable
anyway; the other two should be cited and fetched, not vendored.

---

## Public

- `Scripts/*.R`, `*.Rout` — all analysis code and its console transcripts.
- `Slurm/*.sbatch`, `*.slurm` — job submission scripts. Previously excluded
  wholesale by the old `Slurm/` rule; they are just paths and resource requests
  and are the most reproducibility-relevant artifact in the repo, so they are
  now included.
- `Logs/` — install logs and slurm stdout.
- **Aggregate** result CSVs: `task10_*`, `task13_state*`, `task14_state*`,
  `Task15_*_flagbias_by_state.csv`, `Task17_Chr22_Selected100CpGs.csv`,
  `Task17_Chr22_Biological_vs_External_Normal_100CpGs.csv`.
- **Per-sample count** summaries: `task13_sampleFlagSummary_*`,
  `task14_sampleFlagSummary_*`, `Task15_*_flagbias_by_sample.csv`,
  `Task17_Chr22_BiologicalFlag_SampleSummary.csv`,
  `Task13_SampleFlagSummary.xlsx`. These are 53-row tables of counts per
  de-identified sample with no genomic coordinates attached — a judgment call,
  but they carry no site-level information and are what the figures are built
  from. Move them to private if your IRB reads individual-level strictly.
- Figures: `task9_*.png`, `Task15_*_barplot.pdf`, `BRCA_chr1_*.pdf`.

---

## Two items needing your decision

### 1. `Packages/OutlierMeth/` is third-party code, already pushed

This is GPL-2 code by Bradley Downs and Leslie Cope (JHMI), cloned from
`github.com/bdowns4/OutlierMeth`. It is **already committed and pushed** to
`anamgiri91/bioinformatics-research`.

Redistribution is legally fine — GPL-2 permits it and the `LICENSE` file is
intact — so nothing is on fire. But vendoring someone else's package into your
research repo muddies authorship, and the directory carries its own `.git`
(178 MB of pack files) which git will eventually turn into a broken submodule
reference. Recommended: `git rm -r --cached Packages/OutlierMeth`, and record
the dependency plus its exact commit (`408b110`) in the README instead.
Its `.git/` is already excluded so the 178 MB can never enter your history.

### 2. Confirm the GitHub repo's visibility

`Scripts/` embeds absolute HPC paths including a colleague's home directory:
`/home/s_s355/research3.data/TCGA/filter.BRCA.UCEC.Aug17.2022/...`. That
discloses another user's account name and their data layout. Harmless if the
repo is private; worth parameterising into a single `DATA_DIR` variable before
it goes public. Nothing scanned as a credential — no keys, tokens, passwords,
or email addresses anywhere in the tree.
