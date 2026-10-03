# Task 68: local-only WGBS execution

Updated 2026-10-02 following the user's instruction: **do not use the
cluster**. No further SSH checks or SLURM submissions are part of this
plan. The earlier cluster preparation is retained as a historical record.
The plasma redo remains deferred.

**PAUSED by the user on 2026-10-02 (2026-10-03 UTC).** All 13 related
processes were stopped and the shutdown verified. No job is running or
scheduled to restart automatically. Files and logs are preserved. Read
[the session handoff](../SESSION_HANDOFF.md) for the exact checkpoint and
restart commands; resume only when the user requests continuation.

## Last local run, started 2026-10-02 and paused before alignment

Pilot row 1, **SRR13602252**, reached local setup. When restarted, its driver
waits for verified FASTQs and completed indices, then automatically runs
trimming, alignment, paired deduplication, query-name sorting and the
counts-only coverage check. Row 2 will follow QC review of row 1; it has
not been started. No methylation covariance outcomes have been computed.

- Trim Galore 0.6.10 is installed at commit
  `4edff97d22f3837d42a29e4afbfaeb6e07ffb11b`, with cutadapt 5.2. A synthetic
  paired-read check confirmed the declared 8/20-base start clips, 8/8-base
  end clips and expected filenames.
- The full initial UCSC hg38 FASTA passed the published compressed MD5
  `1c9dcaddfa41027f17cd8f7a82c7293b`. Its autosomal map contains
  27,852,739 reference CpGs. Both bisulfite indices were interrupted and
  must be rebuilt; `.bt2.tmp` files are not complete indices.
- FASTQs download in bounded chunks with four HTTPS connections. All 39
  R1 chunks and 20 R2 chunks are saved, plus four partial R2 chunks. Each
  assembled file must match the ENA MD5 before processing. Chunks and the
  unused partial file from the earlier slow download are retained.
- Stage checkpoints reject altered outputs or mismatched commands/code;
  failed stages receive no success marker. A donor-level lock prevents
  concurrent duplicate pipeline runs. Restart tests passed.
- Each processing stage records elapsed time, output sizes and the
  operating system's resource report. Maximum process RSS in that report
  is not a measurement of the summed peak memory of concurrent children.

Live status and logs are private under `Data/`:

- `gse165915_sperm_wgbs/SRR13602252/LOCAL_PIPELINE_STATUS.json`
- `gse165915_sperm_wgbs/SRR13602252/local_driver.log`
- `gse165915_sperm_wgbs/SRR13602252/download_driver.log`
- `gse165915_sperm_reference/prepare_driver.log`

On successful completion, the driver writes the aggregate report
`Results/Task68_LocalPilot_SRR13602252.json`. Until that exists and QC is
reviewed, no pilot feasibility result is claimed. This is an operational
run record, not an external statistical protocol lock.

Sources: [Trim Galore release](https://github.com/FelixKrueger/TrimGalore/releases/tag/0.6.10),
[UCSC hg38 reference directory](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/initial/).

## Local feasibility checked

- This Mac reports 24 GiB RAM and 15 CPUs; free disk space is about
  684 GiB at this check. Available memory during processing and actual
  peak memory/runtime have not been measured.
- The two metadata-selected pilots total 23.35 decimal GB of compressed
  FASTQs. Process them one at a time, allowing at least 100 GB working
  space for the active donor, plus the reference and its indices.
- All 52 donors total 742.8 decimal GB of raw FASTQs, already more than
  current free disk space before alignments or reference files. Do not
  attempt a full-cohort download. Expansion needs an explicit sequential
  storage/retention design informed by the pilot.
- Bowtie2, samtools and R are available. Bismark 0.24.2 is installed under
  `Data/shared_fragment_tools/`; the local Python environment has numpy,
  pysam and cutadapt. Trim Galore and the hg38 FASTA/CpG map are now set up;
  bisulfite indexing is paused and incomplete. The existing plasma hg19 reference is not a
  substitute for the declared hg38 sperm pipeline.

The pilot is feasible in disk capacity; its execution is **paused**
and resource feasibility is not a WGBS validation result. The separately recorded
depth check from the study's processed CpG reports remains useful, but its
cohort eligibility projection does not replace whole-fragment counting.

## Execution order

1. Finish local full hg38 indexing. Record tool
   versions, reference source and hashes. Build the CpG map from the same
   FASTA as the Bismark indices.
2. Run pilot row 1, inspect its mapping, duplicate, fragment-length and
   coverage reports, and measure runtime, peak memory and working storage.
   Then run row 2 sequentially. Keep all intermediate data private under
   `Data/`; do not automatically delete existing project data.
3. Check whether three-way splitting at the reported low usable coverage
   can support enough eligible donors/pairs. Two donors only assess
   per-donor processing and coverage, not the cohort-wide ≥20-donor
   criterion. Plan sequential cohort processing only after this gate.
4. Before methylation outcomes, finish missingness-aware calibration and
   the versioned external lock. Covariance error remains primary; the
   Task 67 absolute decay contrast remains secondary.

The library, alignment, deduplication and whole-fragment counting rules in
the [earlier preparation record](Task68_ClusterExecution.md) still apply;
its cluster commands do not. Pilot manifests and counts-only scripts are
unchanged. The local wrapper uses four requested alignment threads and
the existing two-core trimming setting, and invokes no remote commands.

Review commands without downloading or processing data:

```sh
bash Scripts/Task68.LocalPilot.sh 1 --dry-run
bash Scripts/Task68.LocalPilot.sh 2 --dry-run
```

After local prerequisites are present, remove `--dry-run` and execute one
row at a time. Each real run checks required tools, reference/map agreement,
ENA checksums and available disk space. Completed runs can only be reused
with matching commands, metadata and counting-code hashes.
