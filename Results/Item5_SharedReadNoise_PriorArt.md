# Item 5: does CpG-specific error correction with independent validation already exist?

Literature assessment, 2026-09-30. This supplements
[the shared-fragment simulation and derivation](Item5_SharedReadNoise_Design.md).
It does not change the frozen score, hypotheses, analyses, or validation results.

**Related implementations and independent assay validation exist. An exact match
to our proposed between-person covariance correction was not verified.** The
earlier phrase “the remaining challenge is a reliable CpG-specific implementation
and independent validation” was too broad if interpreted as a claim that the
field lacks these generally. The evidence supports a narrower research question,
not a demonstrated gap or a claim of priority.

## What would count as an exact match?

The proposed target is the covariance, across people, of two latent methylation
proportions. It differs from within-molecule CpG linkage, methylation-state
calling, differential methylation, and prediction of age or disease.

For donor k, the candidate explicitly uses site coverages N_ik and N_jk, the
number r_k of fragments calling both sites, and joint methylation states on that
shared subset. Under the observation model in the design note it estimates:

\[
\widehat C_{ij}=s_{ij}-\frac{1}{n}\sum_k\widehat c_k,
\qquad
\widehat c_k=\frac{r_k^2}{N_{ik}N_{jk}(r_k-1)}
\left(\widehat p_{11,k}-\widehat p^R_{i,k}\widehat p^R_{j,k}\right).
\]

The shared-fragment expression requires r_k >= 2; zero overlap and one-fragment
overlap need the separate handling described in the design note. A practically
equivalent method need not use this exact algebra, name, or estimator. A joint
likelihood estimating the same latent covariance with the same sampling
dependence could also qualify.

The substantive questions are whether such a method is implemented, handles
weak biological variance and uncertainty reliably, and validates recovery of
between-person covariance using evidence not supported by the same fragments.

## Closest primary sources

### MethylPCA — Chen et al., 2013

The paper distinguishes neighboring CpG correlation from shared DNA fragments
versus biological correlation. Its two-stage grouping was applied to 1,500
MBD-seq samples. The authors report independent-sample, targeted-pyrosequencing
replication of top association findings, with details referred elsewhere.
This is a strong conceptual precedent, but grouping enrichment-assay signals
does not establish the RRBS sampling-covariance correction above. Nor does
replication of associations establish unbiased covariance recovery.
[Paper](https://link.springer.com/article/10.1186/1471-2105-14-74).

### Binomial measurement-error correction — Buonaccorsi et al., 2016

This paper directly addresses bisulfite methylation proportions with binomial
error and varying sequencing effort. It develops modified SIMEX and explores
beta-binomial regression calibration, with simulations and a methylation–age
example. The accessible abstract does not establish shared-fragment correction
of a CpG-pair covariance. Full methods and software availability were not
verified in this search, so that narrower comparison remains provisional.
[Primary abstract and citation](https://pubmed.ncbi.nlm.nih.gov/27264206/).

### epiG — Vincent et al., 2017

epiG implements read-pattern inference accounting for experimental errors and
biological noise. Its CpG-state benchmark uses independent Infinium measurements;
SNP calls are compared with SNP6.0 arrays. The methylation benchmark selects
extreme array betas and optimizes classification thresholds on that benchmark.
This is independent-assay evidence for state inference, not external-cohort
validation of corrected between-person covariance.
[Paper](https://link.springer.com/article/10.1186/s13059-017-1168-4).
The authors provide an [R implementation](https://github.com/vincent-dk/epiG)
and [analysis scripts](https://github.com/vincent-dk/using-epiG).

### dSOMNiBUS — Zhao et al., 2024

The method models mismeasured multivariate binomial outcomes, read depth,
dispersion, and within-person dependence through random effects, with spatially
varying covariate effects. Its target is differential methylation. Validation
includes simulations and a rheumatoid-arthritis sequencing application; these
do not establish our covariance-recovery target.
[Paper](https://onlinelibrary.wiley.com/doi/10.1002/sim.10149).

The [official SOMNiBUS vignette](https://bioconductor.org/packages/release/bioc/vignettes/SOMNiBUS/inst/doc/SOMNiBUS.html)
documents per-site methylated counts, total counts, position, sample ID, and
covariates. It does not list pairwise fragment-overlap or joint-state counts.
**Inference from that interface:** it is not a direct implementation of our
overlap-specific correction. Dependence modeling alone does not establish that
the model separates this sampling covariance from biological covariance.

### MeasurementError.cor — an existing general implementation

The Bioconductor manual documents estimates of both measurement-error correlation
(`cor.me`) and underlying correlation (`cor.true`), using observations and their
standard errors. Its motivating setting is microarray expression, not shared
binary CpG fragments. It demonstrates that implementing a distinction between
error and latent correlation is already established, and supplies a potential
statistical baseline rather than a ready-made CpG solution.
[Official manual](https://bioconductor.org/packages/release/bioc/manuals/MeasurementError.cor/man/MeasurementError.cor.pdf).

## Implications for our claim and evaluation

We cannot claim to introduce CpG measurement-error correction, recognition of
fragment-induced correlation, latent-correlation estimation, or independent
assay validation. These have substantial precedents above. A CpG specialization
of established measurement-error methods may be useful without constituting a
new correlation coefficient.

The unresolved question from this search is narrower: **does explicitly
estimating shared-fragment sampling covariance improve recovery of latent
between-person CpG relationships, beyond existing corrections, with stable
uncertainty and independent validation?** No inspected source established all
of those properties together. That is a search finding, not proof of absence.

For this project, distinguish three validation claims:

1. Independent assay measurements can check CpG methylation estimates.
2. Independent people can check generalization of a locked method.
3. Fragment-disjoint estimates or independent library replicates can help
   test whether apparent CpG covariance is supported by shared sampling noise.

None alone proves all three. Fragment splitting must keep each physical
fragment together across its CpGs. It reduces depth and removes shared sampling
only under suitable independence assumptions; batch effects, PCR dependence,
mapping errors, and cell mixture require separate consideration. The existing
simulation verifies moment identities under its model, not those assumptions
in real data. These are proposed evaluation requirements, not completed tests.

An appropriate current description is: “We investigate a shared-fragment-aware
adaptation of measurement-error correction for between-person CpG covariance.
Its practical advantage and exact relationship to prior estimators remain to
be established.” Independent validation is evidence of utility, not novelty by
itself.

## Search scope and limits

The targeted search covered bisulfite/binomial measurement error, corrected
CpG co-methylation, shared/overlapping fragments, latent correlation, and
read-level methylation models, followed by primary papers and official software
documentation. Representative query concepts were “bisulfite binomial
measurement error correlation,” “CpG shared fragments biological correlation,”
and “methylation covariance sequencing noise correction.”

This was not a systematic review with a preregistered search protocol. The
2016 paper was assessed from its primary abstract; software was not installed
or benchmarked. Documented availability does not establish present-day
installability or suitability. Overlapping paired-end deduplication addresses
double counting a site; it should not be confused with dependence between two
distinct sites measured on one correctly counted fragment.
