# Task 60: prior art for the shared-read noise correction (plan phase H)

Checked on 2026-10-01. This covers every method in section 11 of
`plan_shared_read_noise.md`, plus the nearest work found in other fields.
Not finding a matching paper does not prove the method is new.

**How each was checked.**
- **Full text, or the code itself:** Saito & Suyama, MethylPCA, MPCI,
  Jajoo et al., Saccenti et al., MeasurementError.cor (vignette), CS-CORE
  (source code) and Simon & Coop.
- **Abstract or a published summary only:** Guo et al., DNAmBERT, epiG,
  dSOMNiBUS, Buonaccorsi et al., Buffalo & Coop (as Simon & Coop describe
  it), molecular cross-validation and data thinning.
- **Section 4 (the three further areas):** CompDTUme was checked in its
  source code. sleuth, SOMNiBUS and Affinito et al. were checked from their
  abstracts or published summaries.

## The exact question

Has anyone used each donor's joint states on shared fragments (the
n00, n01, n10 and n11 of two CpGs read together) to estimate the sampling
covariance those shared reads add, and then subtracted it from the
covariance of neighbouring CpGs across people?

## 1. Methylation methods

| Work | What it does | How it differs |
|---|---|---|
| Saito & Suyama 2015, *Epigenetics* 10(12):1093-1098 | Linkage coefficient D = p11 - p_i p_j for two adjacent CpGs "residing in a single sequencing read" | Describes heterogeneity within one sample; no correction across people |
| Guo et al. 2017, *Nature Genetics* 49:635-642 (MHL) | Co-methylation of adjacent CpGs from fully methylated read substrings | Within-sample read structure; no correction across people |
| Nazer, Mohammadzade & Mehrmohamadi 2026, *PLoS Computational Biology* (MPCI) | Weighted Manhattan similarity across reads and consecutive CpGs | Within-sample; no correction across people |
| DNAmBERT 2026, *Briefings in Bioinformatics* 27(5):bbag455 | A transformer that learns from DNA sequence and read-level methylation haplotypes to detect cancer in cfDNA | A prediction model; no covariance estimate |
| Vincent et al. 2017, *Genome Biology* (epiG) | Infers haplotypes of methylation and SNPs from bisulfite reads, checked against independent experiments | Validates states within a sample, not covariance across people |
| Chen et al. 2013, *BMC Bioinformatics* 14:74 (MethylPCA) | States that "in MBD-seq neighboring CpGs will be highly correlated because they are largely covered by the same DNA fragments" | Merges sites with r > 0.9 into blocks; does not estimate or remove the shared-fragment term |
| Zhao et al. 2024, *Statistics in Medicine* (dSOMNiBUS) | Regression for differentially methylated regions, allowing for read depth, experimental error and extra-binomial dispersion | Its target is covariate effects on regional methylation, not covariance between CpGs |
| Jajoo et al. 2023, Research Square preprint | Genome-wide map of correlated methylation | Does not discuss shared-fragment inflation |

## 2. Measurement-error corrections

| Work | What it does | How it differs |
|---|---|---|
| Buonaccorsi et al. 2016, *Statistics in Medicine* 35(22):3987-4007 | Corrects regression for binomial read-count error in one methylation proportion used as a predictor | The variance part of our correction; one site, so no covariance between sites |
| Saccenti, Hendriks & Smilde 2020, *Scientific Reports* 10:438 (author correction 2023) | Corrects the Pearson correlation for measurement error under several error models, including an error term shared by both variables | The error variances and covariance must be supplied. The authors say there is "no simple and immediate approach" to estimate them in the complex cases |
| Ding & Gentleman 2003; Bioconductor MeasurementError.cor | Fits the correlation of true values by maximum likelihood under bivariate normality, from point estimates and their standard errors. It also fits one correlation between the two variables' errors | Needs normality and supplied standard errors; one error correlation for all samples; no use of joint read states |
| Su et al. 2023, *Nature Communications* 14:4846 (CS-CORE) | Co-expression from single-cell counts. The variance model removes a Poisson term; the covariance model has no error term | Counts of two genes come from different molecules, so there is no shared-read covariance to remove. This is the variance part only |
| Buffalo & Coop 2020, *PNAS* 117(34):20672-20680; Simon & Coop 2024, *PNAS* 121(9) | Covariance of allele-frequency changes over time. Two adjacent intervals share one sampled frequency, which biases their covariance. They add back that frequency's sampling variance, p(1 - p)/(n - 1) | The same logic: remove shared sampling noise from a covariance. But there the shared noise is one binomial variance, known from p and n. Here it depends on how the two sites co-vary within fragments, so it needs the joint read states |

## 3. Splitting molecules to check a method

| Work | What it does | How it differs |
|---|---|---|
| Batson, Royer & Webber 2019, bioRxiv 10.1101/786269 (molecular cross-validation) | Splits each cell's molecules at random into two groups, which act as two independent draws. One fits a denoiser, the other checks it | The same idea as the A/B/C benchmark. Here the split is by fragment, in three parts, and the reference is a B-by-C cross covariance |
| Neufeld, Dharamshi, Gao & Witten 2024, *JMLR* 25(57):1-35 (data thinning) | General theory for splitting an observation into independent parts that sum to it, for Poisson, binomial and other convolution-closed distributions | Our fragment split is a case of this. It gives independent parts only if fragments are independent molecules, so leftover PCR copies would break it |

## 4. Further search: the three remaining areas

Done later on 2026-10-01, after the external test.

| Area | Work | What it does | How it differs |
|---|---|---|---|
| Other count data with shared reads | **Young, Van Buren & Rashid 2024, *Biostatistics* 25(2):559 (CompDTUme)** | Differential transcript usage. Each sample's error covariance between transcript proportions is estimated from the quantifier's inferential replicates (Gibbs or bootstrap). The mean of these matrices is subtracted from the across-sample covariance (`UpdatedCovAlt <- SigmaTildeAltNewModeling - mean.withinhat` in the package code). Results with a negative variance are left undefined | **The same moment correction, including the off-diagonal error covariance that ambiguous reads create.** The error covariance comes from resampling a quantification model, not from a closed-form count estimate. The target is a test of transcript usage, not co-methylation |
| Other count data with shared reads | Pimentel et al. 2017, *Nature Methods* 14:687-690 (sleuth) | Uses bootstrap replicates to separate inferential variance from biological variance | The variance part only |
| Models that may be the same | Zhao et al. 2021, *Biometrics* 77(2):424 (SOMNiBUS), and dSOMNiBUS (2024) | Hierarchical binomial regression for regions, with smooth effects along the genome, read depth and error terms | Model latent levels along the genome; no use of the joint read states of two CpGs across donors was found |
| Co-methylation with distance | Affinito et al. 2020, *Genomics* 112(1):144-150 | Co-methylation of nearby CpGs falls with distance, in ultra-deep targeted bisulfite data | Describes co-methylation; does not correct for shared reads |
| Tool documentation | coMethDMR, comb-p, DMRcate | Their inputs are per-site betas, counts or p-values | These inputs carry no joint read states, so they cannot make this correction. Fragment-level tools (wgbstools, the mHap tools) compute within-sample measures |

Not searched further: allele-specific expression, pooled sequencing and
metagenomics. These may hold more examples of the same correction.

## 5. The estimator itself is textbook

Checked in the claims audit below. The per-donor term is

    c_k = (r n11 - (n10 + n11)(n01 + n11)) / (N_i N_j (r - 1))
        = [r / (N_i N_j)] x [(r n11 - n1. n.1) / (r (r - 1))]

The second factor is Weir's unbiased estimator of the linkage
disequilibrium coefficient D from a 2 x 2 table of r haplotypes (Weir
1979, *Biometrics* 35(1):235-254), applied to the r shared reads. The first
factor is the share of the two sites' reads that are shared.

The fact that shared sampling units add covariance in proportion to their
overlap is the standard result for estimates from overlapping samples. It
is used, for example, for GWAS meta-analysis with overlapping subjects (Lin
& Sullivan 2009, *American Journal of Human Genetics* 85(6):862-872).
**So the formula is not new.**

## What was not found

No paper was found that does both of these for methylation:

1. recognises that reads shared by two neighbouring CpGs inflate their
   covariance across people; and
2. removes it, using each donor's joint read states.

The subtraction is not new (CompDTUme does it for transcripts), and the
estimator is not new (section 5).

## A bounded statement of the possible contribution

Rewritten after the claims audit (section 6):

> Reads shared by neighbouring CpGs inflate their covariance across people
> in bisulfite sequencing. In colon, about 40% of the observed covariance
> is this noise. In the test data, with a third of the reads, the observed
> covariance is 2.4 to 3.5 times the noise-free value.
>
> A textbook correction removes it: the overlap fraction times Weir's
> unbiased D, taken from each donor's shared reads. This is a count-based
> special case of correlated measurement-error correction.
>
> It was checked on reads that share no fragments (a split in the spirit of
> molecular cross-validation and data thinning). There were two tests,
> each locked in advance, one of them in another lab's data. In both, the
> corrected covariance matched the noise-free reference and the squared
> error fell by 48%. The correction does not improve correlation estimates
> or the ranking of similar pairs.

What is not new:

- **The estimator.** It is the overlap fraction times Weir's (1979) unbiased
  D. The overlap principle is standard, for example in Lin & Sullivan
  (2009).
- **Subtracting the mean within-sample error covariance, off-diagonal terms
  included, from the across-sample covariance.** CompDTUme does this for
  transcript proportions, with error covariances from resampling.
- **The general idea of removing shared sampling noise from a
  covariance.** Population genetics does this for allele-frequency changes
  (Buffalo & Coop).
- **Correcting a correlation for correlated error, once the error
  covariance is known.** Saccenti et al. do this.
- **Correcting the variance for binomial read-count error.** Buonaccorsi et
  al. do this for methylation, and CS-CORE for single-cell counts.
- **Measuring within-read linkage.** Saito & Suyama, and Guo et al., do
  this.
- **Noting that neighbouring sites share fragments.** MethylPCA does this
  for MBD-seq, and merges such sites.
- **Splitting molecules to evaluate a method.** Molecular
  cross-validation and data thinning do this.

What may be new is narrow, and it is not a method:
- **The bias itself.** The finding is that reads shared by neighbouring
  CpGs inflate their covariance across people in bisulfite sequencing, and
  by how much. No paper reporting this was found.
- **Its effect on a known result.** It changes the apparent fall of
  co-methylation with distance.
- **The demonstration that a textbook correction removes it.** This was
  checked on reads that share no fragments, in two tests locked in
  advance: one in the same study and one in another lab's data.

The error covariance that Saccenti et al. call hard to estimate is easy
here, because the reads record it directly.

## 6. Audit of every novelty claim (2026-10-01)

Done after the second external test, on request. The table covers each
claim in results.md section 16, in this file, in LITERATURE.md, and in the
answer given in the session ("did we find something meaningful?").

| Claim | Where | Verdict | Evidence |
|---|---|---|---|
| Subtracting the mean within-sample error covariance is not new | results.md 16.9; this file | **Holds** | CompDTUme package code (`UpdatedCovAlt <- SigmaTildeAltNewModeling - mean.withinhat`); Buffalo & Coop, as described in Simon & Coop 2024 |
| "What may be new is the closed-form estimate from joint read states" | results.md 16.9 and 16.13; this file; the session answer | **Wrong, corrected** | It is the overlap fraction times Weir's unbiased D (section 5). The overlap principle is standard (Lin & Sullivan 2009) |
| No paper uses joint read states to correct the across-donor covariance of neighbouring CpGs | results.md 16.9; this file | **Holds, as a search result only** | Further checks found nothing (next rows) |
| DMR and smoothing methods already handle this | (a search engine's summary, not a claim of ours) | **Not supported** | The dmrseq vignette notes that neighbouring CpGs are correlated but gives no read-level cause. LuxUS has only replicate and cytosine random effects. BSmooth and DSS model spatial correlation of true levels |
| cfDNA work already shows fragment effects on co-methylation | (a search engine's summary) | **Different quantity** | FinaleMe (Liu et al. 2024, *Nature Communications*) relates mean methylation to fragment length within samples. It does not discuss covariance across samples |
| Part of the apparent fall of co-methylation with distance is technical | results.md 16.3; the session answer | **Holds in our data; not found reported** | Colon: corrected covariance nearly flat within 200 bp. Plasma: observed covariance falls 85% over 200 bp, the noise-free reference about 47% |
| How much of that fall is technical "depends on fragment length" | the session answer | **Too strong, softened** | The cohorts differ in tissue, assay and depth at once, so fragment length is not isolated. What the data show is that the noise term ends where molecules end: it is exactly 0 beyond 150 bp in cfDNA (about 167 bp fragments), and still present at 200 bp in RRBS |
| Shared reads make the difference between two sites less noisy | results.md 16.11 | **Not a novelty claim** | It is the identity Var(e_i - e_j) = v_i + v_j - 2c |
| The three-way fragment split | results.md 16.5 | **Not claimed new** | Molecular cross-validation (2019) and data thinning (2024) are cited as precedent |
| epimutacions fixes its floor at one constant across all CpGs | LITERATURE.md | **Holds** | The package vignette's defaults: `beta$diff_threshold` = 0.1 and `quantile$offset_abs` = 0.15, one value for every CpG |
| The item 5 similarity scores are not new | report.md; Item5_Novelty_Assessment.md and the other Codex notes | **Holds; these make no positive claim** | Their citations were not re-checked one by one |

## What a reviewer could still find

- **More examples of the subtraction** in allele-specific expression,
  pooled sequencing or metagenomics. These would not change the
  conclusion, which already treats the subtraction as known.
- **A methylation model that uses joint read states across donors,** for
  example a bivariate binomial mixed model of the 2 x 2 read table. None
  was found.
