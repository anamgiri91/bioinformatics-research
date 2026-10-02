# Task 66: recoverable physical-fragment identity in the plasma cohort

Access and source audit: 2026-10-01. This follows step 1 of
`plan_shared_read_noise.md`. It audits provenance and public access; it does
not change Task 64's lock or recompute Task 65's endpoint. No raw sequencing
files were downloaded by this audit. A separate pilot may use the access
information below.

## Decision

**Reconstruction is feasible from public paired FASTQ, but is not feasible
from the existing aggregated mHap files alone.** All 46 retained healthy
controls have one public paired run. The smallest donor is already the first
provenance-check donor, GSM4502069. Its two compressed FASTQs total 348 MB.
The full cohort totals 31.285 GB compressed, before trimming, alignments,
reference files and temporary storage. No submitted BAM/CRAM was advertised
in the ENA run records checked; this is an availability observation, not
proof that no alignment exists elsewhere.

The published producer also needs correcting in future prose: mHapBrowser's
paper specifies **mHapTools v1.1 for extraction**, then mHapSuite v2.0 for
downstream tracks. Reading only the current mHapSuite converter does not
establish the historical producer. The corresponding mHapTools v1.1 source
was now checked and confirms the mate-linkage limitation.

## Public access and smallest pilot

The [GSM4502069 GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4502069)
identifies normal plasma cfDNA, title KRp1-108 / alternate ID Normal_40,
and links SRX8181977. Its sole GEO supplement is a per-site methratio file,
which cannot supply physical-fragment identity.

The [ENA project run report](https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA628686&result=read_run&fields=run_accession,experiment_accession,study_accession,library_layout,read_count,base_count,fastq_bytes,fastq_ftp,fastq_md5,submitted_ftp&format=tsv)
was joined to the exact retained experiments in `Task62_PlasmaSamples.csv`.
The 46 matched experiments each have one PAIRED run and both FASTQ URLs and
MD5s. The project has 300 runs overall; the other 254 were not selected.
The saved selected metadata are in
[Task66_PlasmaFragmentRecovery_ENAManifest.tsv](Task66_PlasmaFragmentRecovery_ENAManifest.tsv).

| Pilot field | Value |
|---|---|
| GEO / experiment / run | GSM4502069 / SRX8181977 / SRR11615795 |
| BioSample / BioProject | SAMN14738804 / PRJNA628686 |
| ENA `read_count` | 6,615,769 |
| ENA `base_count` | 1,687,574,912 |
| R1 compressed bytes | 158,955,755 |
| R2 compressed bytes | 189,241,555 |
| Total compressed bytes | 348,197,310 |
| R1 MD5 | `f1fb244173f47163c9845c28cdbb6fe3` |
| R2 MD5 | `e70ab282cade76f7b762ea8e0a4b2661` |

Exact files: [R1 FASTQ](https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR116/095/SRR11615795/SRR11615795_1.fastq.gz)
and [R2 FASTQ](https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR116/095/SRR11615795/SRR11615795_2.fastq.gz).
These are archive-listed URLs and checksums; this audit did not download
and checksum the raw files. Raw FASTQ inspection must confirm synchronized
identifiers and actual read counts and lengths.

GEO states paired-end 100 bp. The archive totals imply approximately 127.5
bases per end if `read_count` counts paired spots. This inconsistency must
be resolved from the FASTQs; do not adopt the earlier “2 x 128” assertion
as a verified uniform read length.

## What was lost in the mHap conversion

The [mHapBrowser methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC10767976/)
describe Trim Galore preprocessing, BSMAP 2.90 alignment using
`-q 20 -f 5 -r 0 -v 0.05 -s 16 -S 1`, and sambamba duplicate marking for
paired WGBS and targeted bisulfite data. The clipping search used the first
250,000 reads, trying 0–20 bases and reverting to no clipping when the
mapping improvement was below 5%. The exact clipping values and execution
logs for this run were not found. The paper says mHapTools v1.1 extracted
the mHaps and mHapSuite v2.0 made downstream tracks.

Verified [mHapTools v1.1 source](https://github.com/butyuhao/mHapTools/blob/073dba69a1aac367196a59857f128d2205544105/convert.cpp),
tag commit `073dba69a1aac367196a59857f128d2205544105`:

- Lines 297–307 require compatible chromosome/name and intersecting CpG
  ranges to merge mates. For ordered nonempty ranges, the alternative
  both-false clause cannot make disjoint ranges pass.
- Lines 309–341 choose the higher-quality call when mates disagree; a tie
  retains the first stored call.
- Lines 633–634 skip duplicate, unmapped, secondary, supplementary and
  failed-QC alignments.
- Lines 729–741 can emit an unmatched mate after more than 5,000 intervening
  processed alignments. Its effect in these samples remains unmeasured.
- Lines 760–770 emit nonmerging mates separately; lines 994–1001 write
  six fields without read name, mate name or fragment identifier.

The matching [current mHapSuite converter](https://github.com/yoyoong/mHapSuite/blob/10d3601d830862e46c3d47bb148c535e206ec38d/src/main/java/com/Convert.java)
also groups by read name internally, merges intersecting CpG ranges, and
aggregates identical outputs. Its
[output class](https://github.com/yoyoong/mHapSuite/blob/10d3601d830862e46c3d47bb148c535e206ec38d/src/main/java/com/bean/MHapInfo.java)
prints chromosome, start, end, state string, multiplicity and strand.
Neither output preserves the mapping from separately emitted mates to a
common original fragment.

This is an identifiability issue, not a missing parser option. For example,
two left records with states 0 and 1 and two right records with states 0
and 1 could have come from fragment states {00,11} or {01,10}. The two
aggregate marginal files are identical, while their shared covariance has
opposite signs. Coordinates and multiplicities alone cannot select the
correct pairing. Randomly re-pairing them would invent information.

## Chemistry and calling rules for the new pilot

The original GEO record describes bisulfite conversion followed by Swift
Methyl-Seq library preparation, 13 PCR cycles and custom Roche NimbleGen
capture. It does not document UMI sequencing, so coordinate deduplication
must not be described as verified molecule-level UMI deduplication.

The manufacturer's [Accel-NGS Methyl-Seq protocol, Appendix E](https://www.genetargetsolutions.com.au/wp-content/uploads/2019/07/Accel-NGS-Methyl-Seq-DNA-Library-Kit-Manual.pdf)
describes an Adaptase tail averaging eight bases at the beginning of R2,
sometimes encountered at the end of R1 for short inserts, and requires
trimming for methylation analysis. This supports testing tail removal;
it does not identify the historical trims used on this donor.

Current [Bismark library guidance](https://felixkrueger.github.io/Bismark/usage/library-types/)
suggests 5-prime clipping of 10 bases on R1 and 15 on R2 for Swift libraries,
with 3-prime clipping considered and M-bias inspection. Treat these as
candidate pilot settings to document before any new covariance endpoint,
not evidence that this exact cohort was previously processed that way.

The pilot should use CIGAR-aware reference/base mapping, preserve QNAME and
read-group identity, and exclude invalid alignments and duplicate pairs
under a recorded rule. Associate mates even when their callable CpG ranges
do not intersect. Union called positions once per physical fragment; retain
uncalled internal CpGs as missing rather than imputing states. For equal
quality conflicting calls, a declared no-call rule avoids arbitrary
mate-order dependence. Keep strand/dyad handling and counts explicit.

## Concrete pilot and gates before a new cohort test

1. Retrieve only the smallest donor's two FASTQs, verify the archive MD5s,
   and audit normalized matching identifiers, lengths, per-cycle quality
   and sequence composition. Record these results before alignment.
2. Align to a complete hg19 reference, retaining mate identity. A genome
   subset is unsuitable for judging unique alignment. Record software
   versions, reference hashes, trim commands, read filtering and duplicate
   policy. Use a small deterministic subset of whole pairs for an initial
   engineering smoke test if needed; that subset is not the scientific
   validation dataset.
3. On the complete pilot library, measure whole-fragment versus converter
   overlap: single-ended callable units, overlapping/disjoint mate CpG
   ranges, conflict resolution, recovered pairs by gap, and recovered
   joint-state counts. Stratify by fragment span, insert size, depth and
   mapping quality. Compare against the existing donor's mHap as a
   processing diagnostic, allowing documented differences from trimming
   and aligner versions.
4. Assign partitions once per retained whole fragment, including disjoint
   mates. Assert no identifier crosses A/B/C. Verify that N, M and joint
   counts are consistent and every fragment contributes at most once per
   site. This resolves the provenance problem for newly reconstructed
   counts only; it does not retroactively validate the old reference.
5. Before processing all 46 donors or inspecting new covariance endpoints,
   freeze the reconstruction and validation rules as a new post-audit
   analysis. Keep Task 65's published files intact. Estimate full-cohort
   compute and storage from the complete pilot. A corrected reanalysis of
   the same 46 donors is not a newly independent replication.

Remaining unknowns: historical donor-specific trims and converter invocation;
any original unpublicized BAM/CRAM; the fraction of physical mates separated
in the existing files; their effect on the Task 65 reference; actual
duplicate-family dependence; and whether recovered counts retain enough
coverage for the locked eligibility rules. These require the raw-data
pilot. Neither a typical cfDNA length nor zero mHap overlap in a distance
bin resolves them.

## Audit reproducibility

- Cohort membership: existing `Results/Task62_PlasmaSamples.csv`; no new
  selection by pair statistics.
- Project ENA metadata: fetched 2026-10-01, SHA-256
  `e64fc88558da52073a45d7c3288743fca78016049a1f2b0277f55c8046bb9c7c`
  for the full 300-run response including `sra_bytes,sra_ftp` fields.
- mHapTools v1.1 `convert.cpp`: SHA-256
  `87ef3c95f2d0423cbd79fc89673aa494017bddf9bff8a347bfe8fdc3eace764e`.
- mHapSuite `Convert.java` at the pinned commit: SHA-256
  `91bed0ae15b2af6fc86b7b6580e2fd203fe1cfe40e90be50d1e2dd5852bf94ff`.
- Small retrieved source/metadata snapshots were staged at
  `/private/tmp/plasma_fragment_access`; temporary snapshots are not
  durable repository records. The selected ENA manifest is retained beside
  this note.

All source URLs above were accessed on 2026-10-01. The frozen Task 61/64
documents were read but left unchanged; the historical producer correction
is recorded here explicitly.
