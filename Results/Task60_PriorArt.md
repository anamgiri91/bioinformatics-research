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

## What was not found

No paper was found that does both of these for methylation:

1. estimates, per donor, the covariance that shared fragments add to two
   CpGs' betas, in closed form from their joint read states; and
2. subtracts its mean from the across-donor covariance and squared
   disagreement of neighbouring CpGs.

The subtraction itself is not new: CompDTUme does it for transcripts.

## A bounded statement of the possible contribution

> The correction is a count-based special case of correlated
> measurement-error correction. Per-donor joint read states give an
> unbiased estimate of the shared-read sampling covariance between adjacent
> CpGs, and removing it gives unbiased covariance and squared disagreement
> across donors under stated assumptions. A fragment-split benchmark, in
> the spirit of molecular cross-validation and data thinning, with a
> bootstrap interval calibrated in simulation, tests it against a
> reference that shares no reads.

What is not new:

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

What may be new is narrow. The first part is the estimate of the error
covariance, the piece Saccenti et al. point to. For bisulfite reads, each
donor's joint read states give an unbiased closed-form estimate, with no
resampling model. The second part is applying the correction to
neighbouring CpGs and checking it on independent fragments. This is closer
to a useful application than to a new statistical method.

## What a reviewer could still find

- **More examples of the subtraction** in allele-specific expression,
  pooled sequencing or metagenomics. These would not change the
  conclusion, which already treats the subtraction as known.
- **A methylation model that uses joint read states across donors,** for
  example a bivariate binomial mixed model of the 2 x 2 read table. None
  was found.
