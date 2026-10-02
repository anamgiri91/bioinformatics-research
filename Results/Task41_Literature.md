# Item 5: literature and defensible interpretation

Reviewed 2026-09-28 for `Data/NN.hg38.18P.forw.chr22.w.header.txt` and the measures documented in `Scripts/Task39.RRBS.StatesAndNeighbours.Sep28.2026.R`. This is a targeted literature check, not an exhaustive novelty review. No external data were downloaded.

**Item 5 is feasible as a comparison of paired methylation agreement and between-patient dependence.** Each genomic pair contributes two vectors aligned across the same 18 patients. The ten levels retain ordering but discard within-bin differences. Agreement answers whether each patient has similar methylation at the two sites; correlation answers whether variation across patients is associated. These are distinct scientific questions.

| Primary source | What it establishes and why it matters here |
|---|---|
| [Chatterjee, *A New Coefficient of Correlation*, JASA 116 (2021), 2009–2022](https://doi.org/10.1080/01621459.2020.1758115); [author manuscript](https://arxiv.org/abs/1909.10140) | The population coefficient measures directional functional dependence. The paper specifies random predictor-tie breaking and the response-tie denominator. It already proposes max symmetrization and averaging over all increasing tie rearrangements. These operations are established methods. |
| [Shi, Drton and Han, *On the power of Chatterjee rank correlation*, Biometrika 109 (2022), 317–333; author manuscript](https://arxiv.org/abs/2008.11619) | Demonstrates limitations in local power against certain alternatives. Its asymptotic results do not by themselves rank all methods at 18 patients, but they rule out assuming that a general dependence measure automatically beats simpler correlations. |
| [Q. Zhang, *On the asymptotic null distribution of the symmetrized Chatterjee's correlation coefficient*, Statistics & Probability Letters 194 (2023), 109759](https://doi.org/10.1016/j.spl.2022.109759); [author manuscript](https://arxiv.org/abs/2205.01769) | Studies the joint/null behavior of the two directions and their maximum. A univariate directional null distribution cannot simply be assigned to the maximum. |
| [W. Zhang et al., *Predicting genome-wide DNA methylation using methylation marks, genomic position, and DNA regulatory elements*, Genome Biology 16 (2015), 14](https://doi.org/10.1186/s13059-015-0581-9) | Direct precedent for comparing neighboring assayed CpGs using between-individual Pearson correlation and mean squared Euclidean differences, stratified by genomic distance/context and compared with background pairs. Their absolute correlation loses direction; retaining signed correlation is useful here. |
| [Clifford et al., *Comparison of clustering methods for investigation of genome-wide methylation array data*, Frontiers in Genetics 2 (2011), 88](https://doi.org/10.3389/fgene.2011.00088) | Compares distance/linkage choices for methylation clustering. No method wins across all simulated settings. This supports scenario-based comparison rather than announcing a universally best distance. Their sample-clustering task is different from the present site-pair task. |
| [Guo et al., *Identification of methylation haplotype blocks aids in deconvolution of heterogeneous tissue samples and tumor tissue-of-origin mapping from plasma DNA*, Nature Genetics 49 (2017), 635–642](https://doi.org/10.1038/ng.3805); [PubMed record](https://pubmed.ncbi.nlm.nih.gov/28263317/) | Establishes tightly coupled adjacent CpGs and methylation haplotype blocks. Read-level co-methylation and across-patient correlation are different observation units; marginal beta ratios do not recover phased methylation haplotypes. |
| [Jajoo et al., *A first-generation genome-wide map of correlated DNA methylation demonstrates highly coordinated and tissue-independent clustering across regulatory regions*, Research Square (2023), version 1](https://doi.org/10.21203/rs.3.rs-2852818/v1) | A primary **preprint**, explicitly labeled unreviewed in its [PubMed record](https://pubmed.ncbi.nlm.nih.gov/37333260/), describing correlated methylation units. It is relevant prior art for correlation-based region discovery; it should not be cited as a published Nature Communications article. |
| [Nazer, Mohammadzade and Mehrmohamadi, *MPCI: A novel metric for quantifying DNA methylation patterns in NGS data*, PLOS Computational Biology 22 (2026), e1014076](https://doi.org/10.1371/journal.pcbi.1014076) | Published March 24, 2026. Uses weighted Manhattan similarities across reads and consecutive CpGs in binary methylation haplotypes. Thus Manhattan-based adjacent-CpG similarity has explicit recent prior art, although MPCI needs read-level patterns unavailable in this beta matrix. |

**Recommended analysis choices (statistical interpretation for this dataset):**

1. Use `1 - mean(abs(level_x - level_y))/9` as an interpretable primary ordinal similarity, alongside exact agreement and within-one-level agreement. Report raw-beta mean absolute difference and RMSE to expose bin-boundary artifacts. These are familiar distance transformations, not new metrics. Hamming/identity matching is a useful analogy to genetic state matching; a methylation decile is not a diploid genotype or a phased allele.
2. Keep signed Pearson, Spearman, and both xi directions as complementary dependence summaries. For example, `y = 1 - x` can be strongly dependent while disagreeing substantially. A perfect constant match deserves maximal agreement even when correlation is undefined. Do not replace undefined correlation with zero.
3. Retain the true adjacent-assayed-site indicator and genomic gap. Consecutive complete sites can skip an incomplete site; that selection is a sensitivity analysis, not a new definition of genomic adjacency. Distinguish assayed adjacency from consecutive CpGs in the entire reference genome.
4. Define each measure's constant-vector behavior separately. Pearson/Spearman need both vectors to vary. For quadratic kappa, one constant vector versus a varying vector, or two different constants, gives zero whenever expected disagreement is positive; identical constants give an undefined denominator. Excluding all constant pairs is an optional analysis policy and should be documented.
5. Xi requires a nonconstant response; constant predictor alone does not make its denominator undefined. Use the paper's tied formula, preserve negative finite-sample estimates, and report tie sensitivity. The no-response-tie maximum at `n = 18` is `16/19`, not 1; this ceiling does not transfer to tied vectors. [Chatterjee manuscript](https://arxiv.org/abs/1909.10140)
6. For exploratory inference, permute the patient pairing of one vector, retaining both marginal vectors and all ties, and recompute the complete statistic, including max symmetrization. This tests patient correspondence, not excess genomic proximity. A separate matched nonneighbor comparison addresses locality. Permutations require patient exchangeability; family, batch, or other structure may require restricted permutations. Use multiple-testing correction for any declared discoveries.
7. Evaluate identity, constant matches/mismatches, shifted levels, inverse relationships, nonlinear dependence, independent pairs with matched marginals, and bin-boundary perturbations. Add leave-one-patient-out stability. These are diagnostic checks, not independent biological validation. Neighboring pairs share sites, so do not treat every genomic pair as an independent replicate for uncertainty estimates.

**Feasibility and a possible contribution.** The provided ratios support all of the agreement/distance summaries and descriptive dependence comparisons above. They cannot distinguish read-depth noise from biological variation without methylated/total read counts, or establish molecular haplotype coupling without read-level data. Shared tissue composition or batch can also generate across-patient association.

A defensible candidate contribution is a reproducible comparison showing when agreement and dependence disagree in sparse, heavily tied, 18-patient RRBS data, with explicit adjacency, constant-vector handling, marginal-preserving nulls, and sensitivity to discretization. Stable boundaries or discordant pairs could become biological hypotheses for later validation. Correlated adjacent CpGs, ten-bin coding, Hamming/Manhattan/Euclidean distances, max xi, and tie averaging should not be individually claimed as inventions.

Targeted searches for Chatterjee/xi with methylation/CpG did not identify a directly matching application in the results inspected. That absence does not establish novelty. An external cohort, read counts, and functional annotation would strengthen any eventual method or biological claim; they are not prerequisites for completing item 5 descriptively.

## Addendum: defining and stress-testing a bounded ordinal score

The clarified item 5 allows defining similarity scores **or** a new correlation. An explicit, useful score can therefore complete the task without claiming a general new correlation. The following derivations concern fixed empirical marginals and a uniform random permutation of patient labels in one vector. They are mathematical calculations for the proposed scores, not empirical results from the RRBS file.

For paired decile vectors `x,y` of length `n`, define:

```text
D    = mean_i |x_i - y_i|                         observed paired distance
D0   = mean_(i,j) |x_i - y_j|                     expected distance after permutation
Dmin = mean_i |sort(x)_i - sort(y)_i|             smallest attainable distance
Dmax = mean_i |sort(x)_i - reverse(sort(y))_i|    largest attainable distance
a    = D0 - Dmin
b    = Dmax - D0
S    = 1 - D/9                                  absolute methylation similarity
```

The extrema follow by exchanging crossed pairings on a line. `Dmin` is the empirical one-dimensional Wasserstein-1 distance; the distribution-function representation is established in [Vallander, *Calculation of the Wasserstein distance between probability distributions on the line*, Theory of Probability and Its Applications 18 (1974), 784–786](https://www.mathnet.ru/eng/tvp4387).

**Exact relationship to prior methods.** With disagreement weights `|i-j|/9`, linearly weighted kappa is

```text
kappa_L = (D0-D)/D0,
kappa_max = a/D0,
kappa_min = -b/D0,                               when D0 > 0.
```

Therefore the initially proposed piecewise score

```text
C = (D0-D)/a  if D <= D0,
    (D0-D)/b  if D > D0
```

is exactly `kappa_L/kappa_max` on its positive branch and `kappa_L/(-kappa_min)` on its negative branch. This is an attainable-bound normalization of weighted kappa. [Rapallo, *Analysis of the Weighted Kappa and Its Maximum with Markov Moves*, Psychometrika (2022)](https://doi.org/10.1007/s11336-022-09844-y) studies weighted kappa with fixed margins, its maximum, and nonuniqueness under linear weights.

Moreover, [Safak, *Min-Mid-Max Scaling, Limits of Agreement, and Agreement Score* (2020 preprint), Definition 5](https://arxiv.org/abs/2006.12904), already defines the same piecewise transformation of agreement relative to its attainable minimum, chance expectation, and attainable maximum. Substitution of `A=1-D/9` reproduces the proposed branches. Safak's worked measure uses categorical agreement, so this establishes direct prior art for the normalization construction, rather than proving that this exact ten-level RRBS application was published. Neither an invention claim nor a new-correlation label is supported.

**Why the piecewise candidate should not be the primary score.** It maps the expected distance to zero, but its permutation expectation need not be zero. Consider 18 patients with each vector containing seventeen 1s and one 10. A random permutation aligns the two 10s with probability `1/18`, giving `C=+1`; all other permutations give `C=-1`. Consequently `E_perm(C)=-16/18`, even under random pairing. This is especially problematic in sparse RRBS profiles.

A preferable single-denominator variant is

```text
C_star = (D0-D)/max(a,b),                        if max(a,b) > 0;
C_star = NA, with an unidentifiable-L1 flag,      otherwise.
```

For every pairing, `Dmin <= D <= Dmax`, so `C_star` lies in `[-1,1]`. Its conditional permutation expectation is exactly zero because its denominator is unchanged by permutation and `E_perm(D)=D0`. Its actual attainable range is `[-b/max(a,b), a/max(a,b)]`; both endpoints need not be attainable at unit magnitude. In the rare-high example it equals `+1` when highs align and `-1/17` otherwise, restoring mean zero. The smallest one-sided exact permutation p-value in that example is still `1/18`, so unit score does not establish statistical evidence.

`C_star` remains exactly `kappa_L / max(kappa_max,-kappa_min)`. It is a transparent rescaling of an established coefficient, not evidence of a newly invented dependence measure. Fixed-margin permutation rankings and p-values coincide with those from `D0-D`, since the positive denominator is fixed.

**Limits common to both candidates.** If supports are separated, for example `x=(1,2)` and `y=(8,9)`, all pairings have distance 7. Thus `Dmin=D0=Dmax=7` and both denominators vanish despite both vectors varying. Linear distances cannot identify patient alignment here. This is consistent with Rapallo's Proposition 4 concerning unchanged linear kappa under certain fixed-margin moves. Small denominators can also amplify small absolute changes, so retain `a`, `b`, sample size, occupied levels, and leave-one-patient-out sensitivity beside the score.

The mechanism is visible through nine threshold indicators: let `A_k=1(x>k)` and `B_k=1(y>k)` for `k=1,...,9`. Direct expansion gives `D0-D = 2 sum_k Cov_emp(A_k,B_k)`, using divisor `n`. Thus these scores aggregate agreement across matched thresholds; they do not detect every form of dependence and positive/negative contributions can cancel. Decomposition of linear kappa into embedded binary tables is established by [Vanbelle and Albert, *A note on the linearly weighted kappa coefficient for ordinal scales*, Statistical Methodology 6 (2009), 157–163](https://doi.org/10.1016/j.stamet.2008.06.001).

**Concrete deliverable.** Define `S` as the main similarity and `C_star` as a supplementary permutation-centered alignment contrast, with degeneracy flags and signed Pearson/Spearman/xi comparators. Assess sparse-profile null calibration, bin-boundary sensitivity, absolute agreement versus alignment, and stability after removing one patient. The useful contribution is an explicit, tested scoring procedure for this dataset with documented limitations. Whether that application offers a publishable methodological advance remains an open validation and literature question.
