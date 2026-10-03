# Start here: shared-fragment CpG project

Saved on 2026-10-02 America/Chicago (2026-10-03 UTC).
**PAUSED at the user's request. All local Task 68 jobs have been stopped.**
No computation should restart until the user requests continuation.

## User's current scope

1. Task 67: investigate how shared-fragment noise changes the apparent
   decay of neighbouring-CpG covariance with genomic distance.
2. Task 68: test generalisation in sperm WGBS, GSE165915, locally.
3. **No cluster, SSH or SLURM. Plasma redo is deferred.**

The working roadmap is `plan_shared_read_noise.md`; root `plan.md` concerns
the separate outlier project. The scientific aim is a validated CpG-specific
application of known measurement-error correction, not a claim of new
generic mathematics or an improved item 5 similarity ranking.

## What is complete

- Earlier Tasks 53–65 and the one-donor Task 66 plasma provenance pilot
  are preserved. Do not rerun or change their frozen endpoints.
- **Task 67 development analysis is complete.** After coarse depth/density
  matching within 1-Mb blocks, the estimated near–far covariance drop
  shrinks by 54.7% in colon and 71.5% in blood. Support is sparse, and both
  intervals for squared-error improvement cross zero. This does not
  establish improved decay accuracy. See `Results/Task67_DistanceDecay_Report.md`,
  `Results/Task67_DecayContrasts.csv` and `Results/Fig48_AdjustedDistanceDecay.png`.
  Four known-target simulation settings gave conservative coverage for
  linear contrasts; nonlinear error and low-depth eligibility calibration
  remain incomplete. These are post-result development analyses.
- **Task 68 metadata and code preparation are complete.** There are 52
  distinct sperm donors/runs, 26 in each exposure group; full raw data total
  742.8 GB. `Results/Task68_SpermPilotSamples.tsv` selects two donors by
  smallest compressed size within exposure group. A separate processed-report
  depth screen is saved in `Results/Task68_SpermDepthCheck*.csv`; its cohort
  projection is optimistic and does not replace fragment-level counts.
- Trim Galore 0.6.10 (commit `4edff97d22f3837d42a29e4afbfaeb6e07ffb11b`),
  cutadapt 5.2, Bismark 0.24.2, Bowtie2 and samtools are available locally.
  `.venv-shared-fragment` has Python 3.11, numpy and pysam.
- Full initial UCSC hg38 FASTA is downloaded and verified against published
  compressed MD5 `1c9dcaddfa41027f17cd8f7a82c7293b`. Its `.fai` and autosomal
  CpG map are built; the map has **27,852,739 CpGs**. FASTA provenance is
  in `Data/gse165915_sperm_reference/reference_fasta.json`.
- Synthetic checks passed for fragment counting, non-overlapping mate
  linkage, r_A=1 exclusion, coordinates/quality, declared TruSeq clipping,
  checkpoint reuse and rejection of failed/changed outputs. See
  `Scripts/Task68.CoverageChecks.py`, `Scripts/Task68.PilotChecks.py` and
  `Results/Task68_LocalSetup.json` (setup-time provenance, not a lock).

## Exact state at shutdown

Thirteen related processes were terminated and a fresh process listing
confirmed none remained. `Results/Task68_ShutdownRecord.json` records the
shutdown; `Results/Task68_PauseSnapshot.json` records saved chunks and state.
No files were deleted. The live-status JSON now says `paused_by_user`.

| Component | Saved state |
|---|---|
| hg38 FASTA, `.fai`, CpG map | Complete and reusable |
| CT and GA Bowtie2 indices | Interrupted; `.bt2.tmp` files remain; **not usable** |
| Row 1 R1, SRR13602252 | All 39 chunks saved, 5,184,397,030 bytes |
| Row 1 R2 | 20 complete chunks, 2,684,354,560 bytes; four partial chunks also retained |
| Assembled FASTQ checksums | Not yet checked; `DOWNLOAD_COMPLETE.json` absent |
| Row 1 trimming/alignment/counts | Not started; `PILOT_COMPLETE` absent |
| Row 2, SRR13602261 | Not started |
| WGBS statistical outcomes | None computed; no external protocol lock yet |

The initial slow sequential download's `R1.fastq.gz.partial` is retained but
unused. The ranged downloader uses `R1.download_parts/` and
`R2.download_parts/`; it skips complete chunks, replaces unfinished chunks,
assembles both FASTQs and verifies ENA MD5s. It deliberately retains chunks.

Bowtie2 indexing **cannot continue from the interrupted computation**.
Rerunning the reference-preparation script reuses the verified FASTA/map
but rebuilds the converted genomes/indices. Bismark warns that existing
converted files and incomplete indices will be overwritten. Do not treat
`.bt2.tmp` files as finished indices, and do not launch `LocalPilot.sh`
directly before reference preparation has completed. The orchestration
driver waits for `REFERENCE_COMPLETE.json` and `DOWNLOAD_COMPLETE.json`.

## Where the files are

- First donor: `Data/gse165915_sperm_wgbs/SRR13602252/`.
- Reference: `Data/gse165915_sperm_reference/`; Bismark directory
  `hg38_bismark/`, CpG map directory `hg38_cpg_map/`.
- First-donor status: `LOCAL_PIPELINE_STATUS.json` under the donor folder.
- Logs: `prepare_driver.log`, `bismark_genome_preparation.log` under the
  reference folder; `download_driver.log` and `local_driver.log` under donor.
- Each completed processing stage will have `stepN.done.json`, `.log` and
  `.resources.txt`. A failed stage gets no success marker. A file lock stops
  duplicate processing of a donor.
- On completion, aggregate QC will appear as
  `Results/Task68_LocalPilot_SRR13602252.json`. It does **not exist yet**.
- Current operational details: `Results/Task68_LocalExecution.md`.

## Restart, only after the user requests it

1. Read this handoff; inspect current processes, disk space, status files
   and git diff. Another session can edit this repository. Do not duplicate
   an existing run or overwrite unrelated changes.
2. Resume the first donor's download and rebuild its reference indices.
   Use **new PIDs**, never the historical PIDs in logs. The following is a
   single-shell example, from the repository root:

   ```sh
   .venv-shared-fragment/bin/python Scripts/Task68.PrepareReference.py \
     >> Data/gse165915_sperm_reference/prepare_driver.log 2>&1 &
   SRN_REF_PID=$!
   .venv-shared-fragment/bin/python Scripts/Task68.DownloadPilot.py --row 1 \
     >> Data/gse165915_sperm_wgbs/SRR13602252/download_driver.log 2>&1 &
   SRN_DOWNLOAD_PID=$!
   .venv-shared-fragment/bin/python Scripts/Task68.RunLocal.py --row 1 \
     --reference-pid "$SRN_REF_PID" --download-pid "$SRN_DOWNLOAD_PID" \
     >> Data/gse165915_sperm_wgbs/SRR13602252/local_driver.log 2>&1
   ```

   Network downloads and macOS `/usr/bin/time -l` resource counters required
   tool sandbox escalation this session. Request the relevant tool permission
   if blocked; do not use the cluster. Repeated preprocessing commands should
   verify existing data rather than redownload it.
3. The driver automatically waits for prerequisites, then runs the row-1
   pipeline: TruSeq clipping; directional full-hg38 Bismark alignment;
   paired deduplication; query-name sort; state-blind fragment-coverage QC.
   It stops at `completed_pending_qc_review`, without starting donor 2.
4. Review the first donor's mapping/deduplication, coverage/eligible pairs,
   fragment-length distribution, runtime and working storage. Then proceed
   to row 2 if technically sound. No extra user approval is needed for
   already authorized local work once they request resumption.
5. Full-cohort processing needs a sequential storage/retention plan. The
   machine has 24 GiB RAM; about 683 GiB was free before indexing. Raw reads
   alone exceed free space. Do not download all 52 donors together or delete
   existing project data automatically.
6. Before methylation outcomes, finish actual-design missingness/uncertainty
   calibration and a versioned external lock. Covariance error remains
   primary; absolute distance-decay contrast is secondary. Retain original
   eligibility rules; a secondary finding cannot rescue a failed primary.

## Version control and caveats

Existing commits through `4884d94` include Task 67, Task 68 preparation and
the processed-report depth check. This shutdown checkpoint adds the local
workflow, documentation and handoff. Large data/reference/download files
stay on disk under ignored `Data/`; they are **not backed up by git**.
The local environment also remains on disk. No remote push is requested.
Unrelated untracked reports, archives and the `correlation for cg sites/`
tree were preserved and should not be swept into this checkpoint.

No new WGBS performance, covariance or biological conclusion is available.
The previous session's progress messages are operational status, not results.
