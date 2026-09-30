# Task 47 addendum: the public cohort for section 5, chosen by the frozen rules

Recorded on 2026-09-29, before any methylation value from this dataset
was downloaded or looked at. Only GEO metadata was used.

## How the candidates were found

An NCBI GEO search (human, "methylation profiling by high throughput
sequencing", RRBS or WGBS in any field) returned 849 series. Of these,
463 have at least 15 samples. The frozen rules prefer a source with
read-level files, so the three series that list read-level file types
were checked first.

| Series | What it is | Read-level files | Build | 15+ people of one tissue? |
|---|---|---|---|---|
| GSE233417 | RRBS of GTEx noncancer tissues, plus white blood cells | `.pat` (fragment level) | hg19 | yes, 14 tissues |
| GSE186458 | WGBS atlas of sorted cell types (Loyfer 2023) | `.pat` | hg38 | no: its largest groups mix cell types or are cfDNA |
| GSE212391 | methylation haplotypes of 11 cancers | `.mhap` | not checked | no: about 10 per cancer |

## The choice

- **Dataset: GSE233417.** It meets every hard rule:
  - human RRBS;
  - 15 or more samples of one tissue from different people;
  - per-CpG files;
  - open access;
  - read-level `.pat` files for section 6.

  It is on hg19, not the preferred hg38. No candidate had both hg38 and
  15 or more people of one tissue.
- **Tissue: colon.** The rule is the tissue with the most distinct donors
  and one sample per donor, judged from GTEx donor IDs in the sample
  titles:

  | Tissue | Samples | Donors |
  |---|---|---|
  | colon | 29 | 29 |
  | blood vessel | 29 | 28 |
  | esophagus | 38 | 30 |
  | white blood cells | 57 | no donor IDs |

  The esophagus has more donors, but 8 of them give two samples, which
  probably come from two sites of the esophagus. The white blood cells
  have no donor IDs, so it cannot be checked that they are different
  people.

## What still has to be checked, and what happens if it fails

- **Read counts.** The GEO description says the `.beta.bed.gz` files
  hold "the methylation levels of all CpGs". If they have no read
  counts, the counts will be rebuilt from the `.pat` files, which record
  every read.
- **Annotation build.** The data are hg19. So the reliability flags use
  the hg19 versions of the same tracks: Bismap k50 and RepeatMasker from
  UCSC hg19. This is an adaptation required by the data, not a change
  to the test.
