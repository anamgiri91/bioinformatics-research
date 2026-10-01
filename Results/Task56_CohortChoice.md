# Task 56: choosing the external cohort for the shared-read noise test

Recorded on 2026-10-01, from GEO metadata only, before any CpG-pair result
from these samples was computed. This follows plan section 8 of
`plan_shared_read_noise.md`.

## What the plan requires

- **Donors:** at least 30 donors of one tissue.
- **Data:** read-level (fragment) information, with enough coverage for the
  three-way split.
- **Independence:** no donor in common with the development data (the 29
  GTEx colon donors), and no duplicate submissions or reused libraries.

## Candidates

The search covered GEO series of human bisulfite sequencing whose files
include read-level formats (`.pat`, `.mhap`, BAM).

| Series | What it is | Read-level | Build | Raw reads | Why not, or why |
|---|---|---|---|---|---|
| GSE233417, colon | development data | `.pat` | hg19 | EGA (controlled) | development only |
| GSE233417, other GTEx tissues | RRBS, GTEx donors | `.pat` | hg19 | EGA | they share donors with colon (esophagus 5 of 30, blood vessel 4 of 28, adipose 9 of 22), so none has 30 independent donors |
| **GSE233417, white blood cells** | **RRBS, 57 samples, not GTEx** | **`.pat`** | **hg19** | **EGA** | **chosen; see below** |
| GSE186458 | WGBS atlas of sorted cells | `.pat` | hg38 | GEO/EGA | no single cell type from 30 donors |
| GSE262273 | serum cfDNA, capture, liver transplant | `.pat` | hg19 | not public | only 27 pre-transplant samples; post-transplant samples mix the graft's DNA |
| GSE200092 | serum cfDNA, capture, radiation | `.pat` | hg19 | not public | several samples per patient, fewer than 30 patients |
| GSE212391 | methylation haplotypes, 11 cancers | `.mhap` | not checked | | about 10 per cancer |
| GSE316309 | WGBS, early lung cancer | `.mhap` | hg19 | | 26 tumours of two types and about 13 normals |

## The choice and its limits

**Chosen:** the 57 white-blood-cell RRBS samples of GSE233417.

**Why:**
- They are the only public set found that meets the donor count and has
  read-level files.
- They are not GTEx samples, so they should be different people from the
  colon donors.
- They use the same RRBS protocol, UMI deduplication and `.pat` format as
  the development data, so the pipeline needs no changes.

**Limits, stated now:**
- **Same lab and pipeline.** This is independence of donors, not of the lab
  or processing.
- **No donor IDs.** Independence is checked from the data in Task 56: blood
  samples matching a colon donor are dropped, and duplicates among the
  blood samples are dropped.
- **Different tissue.** Blood differs from colon, so the correction's
  benefit may differ for biological reasons.
- **Raw reads are controlled-access.** The same provenance limits as in
  section 16.2 of results.md apply.

**If the screen fails its own positive controls, or leaves fewer than 30
blood samples, external validation is declared unmet.** No other set will
be substituted.
