# Task 68: GSE165915 cluster preparation

**Superseded execution route, 2026-10-02:** the user instructed us not to
use the cluster. Do not run the SSH or SLURM steps below. Continue via the
[local-only execution plan](Task68_LocalExecution.md). The methods and
metadata below remain a preparation record; cluster access is no longer
an active dependency.

Updated 2026-10-02. The user selected the sperm WGBS cohort and deferred the
plasma cohort redo. **No cluster job has been submitted.** An SSH check
reached `leap2.txstate.edu`, but `wln26` authentication failed. Cluster
account/access confirmation has been requested; no password or key is
stored in this repository.

## Verified cohort and input manifest

The archived GEO and refreshed ENA metadata join to exactly 52 distinct
donor labels, 52 BioSamples and 52 paired sequencing runs, with 26 donors
from each exposure tertile. This is a combined high/low DDE-exposure sperm
cohort, not 52 unexposed healthy controls. The intended primary target is
across these eligible donors with equal donor weights; exposure-related
variation remains biological variation in that target.

`Task68_SpermSamples.tsv` records the GEO/BioSample/run join, paired HTTPS
URLs, MD5s, file sizes and exposure groups. Total compressed FASTQ size is
742,806,907,554 bytes (742.8 decimal GB), about 4.70 billion read pairs.
Input metadata hashes are in `Task68_SelectionManifest.json`.

The paper and GEO methods report **roughly 4× genome coverage** after
processing. The previous screening note's ~9× was a raw-bases/genome-size
estimate, not an established usable physical-fragment depth. A three-way
split may leave too few eligible donors at many pairs. Do not lower the
20-donor eligibility rule after inspecting covariance results.

Sources: [GEO series](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165915),
[study methods](https://doi.org/10.3389/fgene.2022.929471),
[ENA project](https://www.ebi.ac.uk/ena/browser/view/PRJNA698569).

## Counts-only pilot chosen before outcomes

Choose the smallest compressed run within each exposure group, solely by
metadata. The two pilot runs are SRR13602252 (low exposure, 12.16 GB) and
SRR13602261 (high exposure, 11.19 GB). They are a feasibility check, not a
two-donor validation and not a basis for estimating across-donor effects.

The pilot pipeline:

1. Download both mates and verify ENA MD5s before processing.
2. Trim with Trim Galore: quality 20, minimum length 30, clipping 8/20 bp
   from R1/R2 starts and 8/8 bp from their ends. These are the TruSeq
   defaults in the study's [CpG_Me documentation](https://www.benlaufer.com/CpG_Me/),
   not the Swift plasma settings. The historical exact per-sample settings
   are not asserted to have been reproduced.
3. Align directionally to a complete hg38 reference with Bismark/Bowtie2,
   max insert 1,000 bp. This insert limit is a declared pilot setting;
   inspect mapping/insert-size QC before finalizing the analysis lock.
4. Deduplicate paired fragments by coordinates with Bismark. Preserve
   QNAME and both mates; no UMI independence is claimed.
5. Query-name sort the BAM. Read CpG locations through Bismark XM/XG tags
   and the CIGAR; require base quality >=20. Union the mates' locations,
   so a fragment counts once at each CpG. Unobserved CpGs between mates
   stay unobserved.
6. Hash `(seed,run,QNAME)` to A/B/C before counting. The whole fragment
   stays in one part. Record site depths and A shared-fragment counts.
   Ignore XM case: **no beta, covariance, correlation or pair ranking is
   calculated in this stage**. Shared coverage is available even when
   mates do not overlap.
7. Save private coverage arrays, packed pair-eligibility masks and an
   aggregate fragment/insert-size QC report. Save input metadata, reference,
   CpG-map manifest and counting-code hashes; a completed pilot can be
   reused only with matching commands and provenance and an existing count
   archive. Coverage is an upper bound
   on final state calls because methylation-conflict masking is deferred.

This is a new reference-aware counts module; the hg19-specific plasma
extractor is not reused with hg38 coordinates. Synthetic tests check mate
union, linkage between non-overlapping mates, shared counts, exclusion of
single-fragment overlaps, strand coordinates, clipping, quality filtering,
state-blind counting, orphans and deterministic splitting.

## Cluster execution

Required tools: Python >=3.11 with numpy and pysam, curl, Trim Galore and
cutadapt, Bismark, Bowtie2, samtools, and SLURM. The preflight script records
versions; exact cluster versions and reference/index hashes must be part
of the later external lock. The batch files request the repository's
existing `shared` partition, 8 CPUs, 64 GB and 24 hours per donor, at most
two simultaneous tasks. These requests have not been verified against
the current cluster policy.

Set the four paths to real cluster locations (no secrets):

```sh
export SRN_REPO=/path/to/bioinformatics-research
export SRN_SCRATCH=/path/to/private/sperm-wgbs
export SRN_REFERENCE_DIR=/path/to/full-hg38-bismark-reference
export SRN_CPG_MAP=/path/to/private/hg38-cpg-map
```

The reference folder must contain `hg38.fa` and Bismark CT/GA indices built
from that same reference. Build the CpG map once from the indexed FASTA:

```sh
samtools faidx "$SRN_REFERENCE_DIR/hg38.fa"
python3 Scripts/SharedReadNoise.CoverageQC.py \
  --prepare-map "$SRN_REFERENCE_DIR/hg38.fa" --map-dir "$SRN_CPG_MAP"
bash Scripts/Task68.ClusterPreflight.sh
sbatch Slurm/Task68.SpermPilot.sbatch
```

Pilot FASTQs total 23.35 GB. Allow at least 100 GB free scratch per
concurrent donor for raw reads, trimmed reads and intermediate BAMs, plus
reference/index space. A full cohort retains much more than 743 GB after
alignment; check quota and retention policy before expanding. Nothing is
deleted automatically by these scripts.

After the pilot QC supports proceeding, the full counts-only array is
`Slurm/Task68.SpermCohortCounts.sbatch`. It reuses identical completed
pilot runs. Aggregate feasibility with:

```sh
python3 Scripts/Task68.CombineCoverage.py \
  --manifest Results/Task68_SpermSamples.tsv --work "$SRN_SCRATCH" \
  --output "$SRN_SCRATCH/cohort_coverage.json"
```

The cohort cannot pass a validation criterion at this stage. Before any
pairwise outcome is evaluated, require complete provenance/QC, enough
pairs with at least 20 eligible donors, simulations using the actual
coverage/missingness design, and a finalized versioned protocol and hashes.
Use the prior covariance-error endpoint as primary and the absolute Task
67 decay contrast as secondary. Sparse coverage, near-constant sperm
methylation and exposure-group covariance must be reported, not used to
select favorable pairs after seeing outcomes. Fragment length must be
measured from alignments rather than assumed from the assay name.

## Validation performed locally

- Manifest joins: 52 unique donors/runs, two FASTQs and valid MD5s each;
  exposure counts 26/26.
- Dry-run command generation: `Task68_ClusterDryRun.json`.
- Synthetic BAM coverage tests passed.
- SLURM and shell scripts pass `bash -n`.
- No real sperm alignment or cluster resource measurement yet; SSH access
  is the immediate blocker.

## Depth check from the study's own CpG reports (2026-10-02)

Cluster access is blocked, so the pilot's main question, usable depth,
was checked here instead. GEO publishes a Bismark CpG report for each
sample, and the reports for the same two pilot donors (GSM5058020 and
GSM5058029, about 260 MB each) were used.
- **Read:** depth only, methylated plus unmethylated counts with both
  strands combined, on autosomes.
- **Not computed:** no methylation level or pair statistic.
- **Script:** `Scripts/Task68.SpermDepthCheck.Oct02.2026.R`. Its outputs
  are `Results/Task68_SpermDepthCheck.csv` and `_UsablePairs.csv`.

| | GSM5058020 (low exposure) | GSM5058029 (high exposure) |
|---|---|---|
| Mean depth per CpG | 2.91 | 2.60 |
| Median depth | 2 | 2 |
| CpGs with at least 6 reads | 15.8% | 12.8% |
| CpGs with at least 10 reads | 2.2% | 1.6% |
| Chance the donor is eligible for a pair within 200 bp | 3.9% | 3.1% |

- **Usable depth is below the published 4x:** about 3x per CpG after the
  study's own processing.
- **Projection to 52 donors.** About 608,000 of the 23.9 million pairs
  within 200 bp (2.6%) would have at least 20 eligible donors. This
  assumes the other 50 donors are covered like these two.
  - Of the pairs with a projected chance above 0.5, 78% lie outside
    repeats, at a typical depth of 8.5 reads (90th percentile 11).
  - They span 2,759 1-Mb blocks. So they are the better-covered parts of
    the genome, not repeat artefacts.
- **The projection is optimistic.** It assumes close CpGs share all
  their reads. Pairs that share fewer reads need both sites deep in
  every part, which is rarer.
- **What it means.** The sperm test would rest on a few percent of
  pairs, from the better-covered regions, not on most pairs as in Tasks
  59 and 65.
- **The alternative.** GSE173787 (53 CD19+ B-cell samples, about 15x
  raw, 316 GB) would make most pairs usable. But it mixes multiple
  sclerosis patients and controls.
- **The whole-fragment pilot is still needed**, now locally, to confirm this with
  whole-fragment counts from our own pipeline.
