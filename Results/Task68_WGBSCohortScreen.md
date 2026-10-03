# Task 68: screening WGBS cohorts for an independent test (metadata only)

**Update, 2026-10-02:** the user selected GSE165915, deferred the plasma
redo, and subsequently instructed us not to use the cluster. Follow the
[local-only execution plan](Task68_LocalExecution.md). The [preparation record](Task68_ClusterExecution.md)
supersedes the provisional coverage and processing assumptions below.
Published usable coverage is roughly **4×**, whereas the earlier ~9× was
computed from raw bases. The sperm libraries require TruSeq-specific
trimming, not the Swift plasma settings. Two counts-only pilots and a
52-donor array were prepared; no cluster job was submitted. The array is
not the current execution route. Counts-only feasibility precedes the final outcome
protocol lock; no sperm methylation outcomes have been inspected.

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
| **GSE165915** (PRJNA698569) | Sperm WGBS, men from the Faroe Islands (George Washington University) | **52** | **Public, paired 2 x 150** | about 90 M pairs each; published usable coverage ~4×; **743 GB** of FASTQ | **Selected.** Split-coverage feasibility is untested |
| GSE173787 (PRJNA727170) | CD19 B cells from blood, multiple sclerosis study | up to 133, patients and controls mixed | Public, paired | about 165 M pairs each; 817 GB | Second choice. Disease status is mixed, and it needs a cluster |
| GSE186888 | White blood cells, WGBS, "CNVS-NORM" donors | 46 | **None linked in GEO** (no SRA relation) | not known | Fails: whole fragments cannot be rebuilt from GEO |
| GSE107729 (PRJNA421218) | Sorted brain neurons and oligodendrocytes, schizophrenia study | 53 NeuN+ (25 controls) | Public, paired | about 463 M pairs each; **15.4 TB** | Too large; few controls of one type |
| GSE186458 | Human methylome atlas, sorted cells | at most 36 per group, many cell types | Public | about 30x each | Same lab and pipeline as GSE233417; not independent |
| CNP0003513, CNP0005464 (CNGB) | cfDNA WGBS projects | not read | the portal returned only a page shell | not known | Not assessed |

## The selected cohort

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
- **Depth may be limiting.** At roughly 4× usable coverage, equal splitting
  averages only about 1.3 reads per part. Eligibility requires at least two
  reads at both sites in every part, so feasibility must be measured from
  counts rather than inferred from raw sequencing yield.
- **Exposure groups.** The cohort was chosen for high organochlorine
  exposure. That is a population trait, not a mixture of tissues, but it
  should be recorded as a covariate.

**The complete raw cohort cannot be retained on this laptop's current free
disk.** It means 743 GB of downloads and about 4.7 billion read pairs to
align. The selected two-donor pilot totals 23.35 GB and can be processed
sequentially. Measure actual runtime and working storage locally before
designing sequential full-cohort processing; no cluster will be used.

## Current execution gates

1. Verify local tools, reference and storage; run the two
   metadata-selected counts-only pilots with TruSeq-specific trimming,
   paired deduplication and whole-fragment partitions.
2. After pilot QC, measure coverage eligibility across all 52 donors.
   Keep methylation states out of the feasibility summaries.
3. Calibrate the interval for the observed coverage/missingness design and
   sperm-like spread. Finalize conflict handling and use the Task 64
   covariance-error endpoint, with the Task 67 absolute contrast secondary.
4. Hash the sample list, code, settings and reference, and commit the final
   lock before computing or inspecting methylation outcomes. Do not claim
   this metadata screen or the counts-only pilot is a completed external test.
