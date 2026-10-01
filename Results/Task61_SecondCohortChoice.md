# Task 61: choosing a second external cohort, from another lab

Recorded on 2026-10-01, before any CpG-pair result from these samples was
computed. This answers item 2 of results.md section 16.11. The blood test
(Task 59) used samples from the same study, lab and pipeline as the
development data. This cohort removes that limit.

## What is needed

- **Donors:** at least 30 donors of one tissue.
- **Data:** read-level methylation, with each counted unit an independent
  molecule. That needs deduplicated reads, and mates merged into one
  fragment.
- **Independence:** a different lab and pipeline from GSE233417, and no
  donor shared with the development data.
- **Practical limit:** no raw-read alignment. This machine has no
  aligners, and aligning 30 or more bisulfite libraries would take days.

## Where read-level files from other labs exist

mHapBrowser (Shanghai Institute of Biochemistry and Cell Biology; Nucleic
Acids Research 2024, doi 10.1093/nar/gkad881) holds read-level files
(`.mhap`) for 5,852 public human samples (hg19). It reprocessed each one
from raw reads with its own pipeline: BSMAP v2.90 alignment, sambamba
duplicate marking for paired-end WGBS and targeted data, and mHapSuite
`convert`.

The mHapSuite source code (`Convert.java`) was read for this record. It
does three things that matter here:
- it skips reads flagged as duplicates;
- it merges the two mates of a fragment into one haplotype when their CpG
  ranges overlap;
- it collapses identical haplotypes into one line with a count.

## Candidates with 30 or more samples of one tissue

Taken from mHapBrowser's own sample table.

| Series | What it is | Assay | Usable donors | Why not, or why |
|---|---|---|---|---|
| **GSE149438** | **Plasma cfDNA, EpiPanGI Dx (Goel lab; Kandimalla et al. 2021, PMID 34465601)** | **Targeted bisulfite capture, paired 2 x 128, NovaSeq** | **46 healthy controls** | **Chosen: see below** |
| GSE107729 | Sorted brain neurons (NeuN+), schizophrenia study (Yi lab) | WGBS | 25 controls (53 with patients) | Fewer than 30 controls. Longer fragments mean many mates are not merged. Files would total about 100 GB or more |
| RRBS series (e.g. SRP200298 blood, GSE87893 adipose, GSE72606 sperm) | Various | RRBS | 40 to 200 | mHapBrowser does not deduplicate RRBS, so PCR copies remain. In simulation these bias the correction (results.md 16.4) |
| CCLE | Cancer cell lines | RRBS | many | Cell lines, not donors; not deduplicated |

## The choice and its limits

**Chosen:** the 46 healthy-control plasma samples of GSE149438.

**Why:**
- **A different lab, assay and pipeline:**
  - the lab is the Goel lab, not the GTEx study's;
  - the assay is targeted capture of cfDNA, not RRBS;
  - the processing is mHapBrowser's, not wgbstools.
- **Deduplicated paired-end data.** cfDNA fragments are short (about 167 bp)
  against 2 x 128 reads, so the mates usually overlap and are merged.
- **Enough donors** (46 healthy controls), with typical depth of about 16
  reads at CpGs outside islands and 12 inside.

**Format check on one file**, done as a provenance check (plan phase A),
before the choice was final. The file was SRX8181977, a healthy control:
- 1,045,174 lines, carrying 1,453,442 molecules.
- Every autosomal line starts and ends on an hg19 CpG (1-based position of
  the C).
- Every haplotype string has exactly one 0 or 1 per CpG in its span.
- Spans reach 288 bp, so merged mates are present.
- No pair statistic was computed.

**Limits, stated now:**
- **Unmerged mates.** If the mates' CpG ranges do not overlap, a fragment
  becomes two lines. Those lines are not seen as shared, and the split can
  send them to different parts, so the reference is no longer free of
  shared reads for those pairs.
  - For a fragment shorter than 256 bp, this can only affect a pair farther
    apart than the overlap of the two mates.
  - For longer fragments, pairs that straddle the gap between the mates can
    be affected.

  A sensitivity analysis on close pairs is therefore declared in the lock.
- **Deduplication by position.** cfDNA fragments often share ends.
  Deduplicating by position may merge different molecules, which lowers
  depth but does not create shared reads.
- **Library chemistry.** The Swift Methyl-Seq library adds a short tail to
  the start of read 2. Whether mHapBrowser trimmed it is not stated.
  Untrimmed tails act like call errors, which simulation showed bias the
  correction slightly.
- **Biology close to blood.** cfDNA from healthy people comes mostly from
  blood cells. This tests a different lab, assay and pipeline more than a
  different tissue.
- **Panel.** Only CpGs inside the capture targets are covered, so there
  will be fewer pairs than in Task 59.
- **No donor IDs.** Duplicates within the 46, and any overlap with the
  development donors, are checked from the data in Task 62 before the lock.

**If the screen fails its own controls, or leaves fewer than 30 samples,
this second external test is declared unmet.** No other set will be
substituted.
