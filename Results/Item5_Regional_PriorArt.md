# Item 5: regional dependence, missingness and novelty

Targeted primary-literature check, 2026-09-29. This report supports development;
it is not an exhaustive priority search or evidence that a new method works.

The strongest candidate contribution is a carefully defined **small-sample
regional dependence test with patient-influence diagnostics**, rather than a
new name for a weighted pairwise correlation. Pooling nearby CpGs, weighting
by genomic separation, and pooling nonlinear dependence statistics each have
strong precedents. Any novelty claim must identify and validate a property
that those precedents do not already supply.

## Six closest methodological precedents

| Primary paper | Established mechanism | Consequence for our claim and benchmark |
|---|---|---|
| [Sofer et al. (2013), A-clustering, Bioinformatics](https://doi.org/10.1093/bioinformatics/btt498) | Agglomerative clustering constrained to adjacent CpGs/clusters; correlation-derived distances, alternative linkage rules, and an optional base-pair separation restriction. Its second-stage exposure analysis uses GEE. | Adjacent correlation plus physical proximity is established. Use the clustering stage as a regional comparator; a phenotype-free analysis should not claim to reproduce its exposure test. |
| [Gatev et al. (2020), CoMeBack, Bioinformatics](https://doi.org/10.1093/bioinformatics/btaa049) | Defines co-methylated regions using across-person correlation, genomic proximity, and background CpGs, including sites not measured on the array. The first genomic-chain stage is independent of the methylation values. | Correlation plus distance plus density is established. Our RRBS input positions are measured/retained positions, not a complete reference CpG map; an approximation using those positions must not be called the original CoMeBack method. |
| [Gomez et al. (2019), coMethDMR, Nucleic Acids Research](https://doi.org/10.1093/nar/gkz590) | Within a predefined region, identifies contiguous co-methylated subregions using correlation of a CpG with the sum of the remaining CpGs (the rdrop statistic). A random-coefficient mixed model then tests association with a continuous phenotype. | Pooling weak site information through a regional summary is established. We can compare rdrop-based regional coherence without needing a phenotype. Its official package also contains explicit missingness filtering, so merely accepting NA values is not a novelty claim. [Package manual](https://master.bioconductor.org/packages/release/bioc/manuals/coMethDMR/man/coMethDMR.pdf). |
| [Nustad et al. (2022), Modeling dependency structures in 450k DNA methylation data, Bioinformatics](https://doi.org/10.1093/bioinformatics/btab774) | Models spatial methylation dependence using a Gaussian random field with Matérn covariance; the selected smoothness gives exponential distance decay. Uses Bayesian spatial inference. | An exponential base-pair kernel is established. An empirical dependence score differs from a generative covariance model, but that distinction alone is not a methodological advance. |
| [Yao, Zhang and Shao (2018), Testing Mutual Independence in High Dimension via Distance Covariance, JRSS B](https://doi.org/10.1111/rssb.12259) | Constructs a pooled test from pairwise distance covariance, including a banded-dependence extension. Targets weak, distributed nonlinear dependence. The discussion explicitly suggests permutation critical values for small-sample size distortion. | A sum of distance-covariance terms, a local graph/band restriction, and permutation calibration are close existing ideas. This is the strongest objection to claiming the regional sum itself as new. |
| [Pfister et al. (2018), Kernel-Based Tests for Joint Independence, JRSS B](https://doi.org/10.1111/rssb.12235) | dHSIC tests joint independence of several possibly multivariate variables using product kernels; the paper develops permutation, bootstrap and approximate calibration. | A multivariate nonlinear regional test with patient permutations is established. A local pairwise sum is less general: it can miss higher-order dependence with independent pairs. dHSIC should be a comparator on common complete observations when computationally practical. |

Additional close statistical prior art is [Jin and Matteson (2018)](https://doi.org/10.1016/j.jmva.2018.08.006),
which develops multivariable dependence measures including sums of squared
distance covariances. These sources make a broad "first regional distance
correlation" claim untenable without much stronger qualification.

## A specific research target worth testing

**Target:** Detect a spatially coherent region whose dependence across patients
survives removal of any one patient, at a calibrated region-level false-positive
rate, using beta values and positions alone. Quantify how often that decision
changes when the available observations are reduced by realistic missingness.

This target addresses two concrete weaknesses of the present dataset: only
18 patients are available, and a single unusual patient can produce an extreme
pairwise correlation. It is not yet established that this target is absent
from all previous work. Its value would come from a demonstrated improvement
in reproducibility at matched error rates, with a precise inferential result.
It deliberately prioritizes distributed support over detecting a true effect
present in only one patient; that tradeoff must be explicit.

A possible statistic, offered as a research design rather than a validated
formula, is

\[
T_R^{\rm stable}=\min_{k\in\{0,1,\ldots,n\}}
 \frac{\sum_{(i,j)\in E_R^{(-k)}}w_{ij}^{(-k)}
       \widehat{\operatorname{dCov}}^{\,2}_{U,ij,(-k)}}
      {\sum_{(i,j)\in E_R^{(-k)}}w_{ij}^{(-k)}}.
\]

Here k=0 means the full sample; the other terms omit one patient. E_R is a
prespecified local CpG graph; U denotes the unbiased, U-centered estimator;
negative estimated values are retained. Weights can depend on genomic gaps
and observation counts, with all rules fixed before testing. The statistic
has units of squared beta difference and is **not a correlation coefficient**.
Normalization or studentization would introduce further choices needing
separate study. Undefined edges and regions require a fixed eligibility rule.

The minimum is a generic worst-case deletion construction, not a claimed new
mathematical operation. Its use could reduce dependence on one patient, but
can also reduce power and alter rankings. Compare it directly to the same
regional statistic without the minimum; otherwise an apparent gain cannot
be attributed to the influence rule. Do not interpret a stable region as
causal regulation: shared cell mixture, batch or other patient factors remain
possible explanations.

## Calibration details that determine whether the method is defensible

1. **Preserve shared-site structure.** A single permutation replicate must
   assign one shuffled patient vector to each CpG and reuse that vector for
   every edge containing it. Independently shuffling each edge creates
   incompatible realizations for shared CpGs. Use independent permutations
   across CpG vectors under the regional mutual-independence null. Applying
   the same permutation to all CpGs leaves their dependence unchanged.
2. **Repeat the statistic completely.** Recompute the deletion minimum, all
   data-dependent normalizations, and any methylation-based region selection
   for each null replicate. Regions fixed using coordinates alone avoid the
   last selection issue. A selected region's unadjusted fixed-region p-value
   does not account for having searched for that region.
3. **State the null accurately.** Independent CpG-wise permutations justify
   a mutual-independence null under exchangeable patients. A sum over selected
   edges may have no power against an alternative involving only untested
   edges or purely higher-order dependence. It is not a universal dependence
   coefficient merely because its null is global.
4. **Missingness needs an assumption, not an imputation trick.** Permuting
   observed values within each site's fixed observed-patient set preserves
   the mask. It is valid only under the required conditional independence
   and exchangeability of the observed values given that mask; MCAR with
   independent identically distributed patients is a sufficient setting.
   Patient-specific distribution shifts or informative dropout can violate
   this. Neither pairwise completion nor permutation automatically solves
   informative missingness. Start validation on complete observations and
   then simulate fixed-mask MCAR and explicit violation scenarios.
5. **Multiplicity must match the advertised result.** A calibrated p-value
   for one predefined region is not genome-wide error control. If max-T is
   used for a scanned family, use coherent full-data permutations across
   overlapping regions. Its strongest simple interpretation is control
   under the complete independence null; strong control under mixed nulls
   needs additional justification.
6. **Do not overinterpret rare patterns.** A leave-one-patient-stable signal
   is not proof of independent replication. The same samples still contribute
   to all deletion fits. External data or a truly independent patient split
   is needed for a separate validation claim.

These are mathematical design requirements and inferences for our proposed
workflow, not claims that each cited CpG method implements them. General
permutation consistency for distance covariance, HSIC and dHSIC is studied by
[Rindt, Sejdinovic and Steinsaltz (2021)](https://doi.org/10.1002/sta4.364).

## Minimal evidence needed for an advance

- Freeze genomic windows/graphs and parameters before evaluating the test.
  Distinguish true consecutive input sites from bridges over missing sites.
- Compare weighted and unweighted regional dCov, Pearson/Spearman regional
  summaries, rdrop, and the existing pairwise Manhattan and raw-xi baselines.
  Use existing software only when input requirements are met; otherwise label
  an implementation as a simplified comparator.
- Use the same null replicates and report uncertainty on rejection rates.
  Include independent continuous values, tied/sparse values, independent
  values with correlated masks, weak distributed regional dependence,
  nonlinear dependence, one-patient spikes and structured dropout.
- Measure power at matched actual false-positive rates and region-boundary
  recovery when truth is known. Held-patient absolute agreement answers a
  different question and naturally rewards Manhattan agreement.
- On the real matrix, report regional support, eligible patients/edges and
  deletion sensitivity. Without labels or independent data, real-data
  examples demonstrate feasibility rather than biological superiority.

**Defensible current wording:** "We investigate a spatially constrained,
patient-influence-aware regional dependence test for sparse small-cohort
methylation data." The influence mechanism and missingness assumptions need
both prior-art checking and successful evaluation before calling the method
novel. A negative result against a simpler comparator should be retained.
