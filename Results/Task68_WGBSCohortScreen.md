# Task 68: screening WGBS cohorts for an independent test (metadata only)

Recorded on 2026-10-01. This is step 3 of the next steps in
`plan_shared_read_noise.md`: choose an independent whole-genome bisulfite
cohort from metadata, with recoverable whole-fragment identity and enough
donors and coverage.

- **Who fetched what.** Codex fetched the GEO and ENA metadata in
  `Results/Task68_Metadata/`. Claude added the run report for the sperm
  cohort (`PRJNA698569_ena.tsv`) and wrote this record.
- **What was not done.** No sequencing data were downloaded for any
  candidate, and nothing was computed from methylation values.
- **This is not a lock.** The choice and its protocol still need to be
  frozen before any data are processed.

## What a candidate needs

- At least 30 donors of one tissue or cell type.
- Public paired raw reads, so whole fragments can be rebuilt by read name.
  Processed files are not enough (Task 66).
- Enough depth for the three-way split: at least 6 reads at both sites in
  at least 20 donors, for many CpG pairs.
- No overlap with the cohorts already used: GSE233417 (colon and blood) and
  GSE149438 (plasma). A different lab from both is preferred.

## Candidates

| Series | What it is | Donors of one type | Raw reads | Size | Verdict |
|---|---|---|---|---|---|
| **GSE165915** (PRJNA698569) | Sperm WGBS, men from the Faroe Islands (George Washington University) | **52** | **Public, paired 2 x 150** | about 90 M pairs each (about 9x); **743 GB** of FASTQ | **Best candidate.** It needs a cluster |
| GSE173787 (PRJNA727170) | CD19 B cells from blood, multiple sclerosis study | up to 133, patients and controls mixed | Public, paired | about 165 M pairs each; 817 GB | Second choice. Disease status is mixed, and it needs a cluster |
| GSE186888 | White blood cells, WGBS, "CNVS-NORM" donors | 46 | **None linked in GEO** (no SRA relation) | not known | Fails: whole fragments cannot be rebuilt from GEO |
| GSE107729 (PRJNA421218) | Sorted brain neurons and oligodendrocytes, schizophrenia study | 53 NeuN+ (25 controls) | Public, paired | about 463 M pairs each; **15.4 TB** | Too large; few controls of one type |
| GSE186458 | Human methylome atlas, sorted cells | at most 36 per group, many cell types | Public | about 30x each | Same lab and pipeline as GSE233417; not independent |
| CNP0003513, CNP0005464 (CNGB) | cfDNA WGBS projects | not read | the portal returned only a page shell | not known | Not assessed |

## The choice, if a cluster run goes ahead

**GSE165915 (sperm, 52 men)** best meets the requirements:
- one tissue;
- more than 30 donors;
- public paired reads;
- moderate depth;
- a lab and population unrelated to the cohorts already used.

**Limits, stated now:**
- **Sperm is unusual.** Most of the genome is highly methylated, so many
  CpG pairs will be nearly constant. Constant pairs carry little
  shared-read covariance, so the test's power will rest on the variable
  pairs.
- **Depth is moderate.** At about 9x, the split leaves about 3 reads per
  part. Many CpG pairs will fail the eligibility rule, so the pair count
  could be far below that of the earlier tests.
- **Exposure groups.** The cohort was chosen for high organochlorine
  exposure. That is a population trait, not a mixture of tissues, but it
  should be recorded as a covariate.

**Not feasible on this laptop.** It means 743 GB of downloads and about
4.7 billion read pairs to align. At typical Bismark speeds that is several
hundred CPU hours, so it is a job for the leap2 cluster. On the laptop, a
pilot of one or two donors could check depth and fragment structure
first. That would cost about 15 GB and a few hours per donor.

## Before any cluster run

1. Freeze the reconstruction rules from the Task 66 pilot: trimming,
   alignment, duplicate policy, mate linkage and conflict handling.
2. Write the locked protocol: the Task 64 endpoint, with partitions
   assigned per whole fragment. Calibrate the interval for this design
   (depth about 9x, sperm-like spread) in simulation first.
3. Hash the sample list, code and reference, and commit the lock before
   processing any methylation data.
