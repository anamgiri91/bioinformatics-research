# Manhattan similarity, genomic distance and additional factors

Primary-source review, 2026-09-29, for the user's proposed extension of item 5.

## Direct answer

Yes. There are close published precedents for both parts: Manhattan-based
methylation similarity and combining methylation similarity/correlation with
base-pair distance and biological context. I did not establish that a single
paper implements the exact proposed combination of ten-level Manhattan
similarity, statistical distance correlation, and genomic distance in the
18-sample RRBS setting. Absence of such a match in this targeted search is not
a novelty finding.

| Primary paper | What is actually combined | Relationship to this project |
|---|---|---|
| [Nazer et al., MPCI, PLOS Computational Biology (2026)](https://doi.org/10.1371/journal.pcbi.1014076) | Weighted Manhattan similarities across sequencing reads and consecutive CpGs; methylated and unmethylated patterns contribute to a regional score. | Closest direct Manhattan precedent. Uses binary read-level haplotypes, not an 18-patient beta matrix. Our file cannot reconstruct its input. |
| [Gatev et al., CoMeBack, Bioinformatics (2020)](https://doi.org/10.1093/bioinformatics/btaa049) | Across-person methylation correlation, physical probe separation and background genomic CpG density. | Closest precedent for correlation plus genomic context. It uses a reference-genome CpG chain (at least one CpG every 400 bp), endpoints under 2 kb apart, and a correlation requirement. It is not simply Manhattan times a distance penalty. |
| [SMART, Nucleic Acids Research (2016; online 2015)](https://doi.org/10.1093/nar/gkv1332) | Methylation specificity state, entropy-based similarity, Euclidean similarity and separation of at most 500 bp; segment merging also limits intervening CpGs. | Particularly close to the idea of combining methylation state, similarity and base-pair distance. These are joint criteria rather than the exact scalar proposed here. |
| [Nustad et al., Bioinformatics (2022; online 2021)](https://doi.org/10.1093/bioinformatics/btab774) | Models residual methylation using a Gaussian random field with a Matérn covariance; the chosen smoothness gives exponential decay with base-pair separation. Examines technical and biological associations with the spatial parameters. | Direct precedent for a fitted correlation-versus-genomic-distance model. Its within-person spatial process differs from across-person pairwise distance correlation. |
| [Sofer et al., A-clustering, Bioinformatics (2013)](https://doi.org/10.1093/bioinformatics/btt498) | Correlation-based clustering of neighboring methylation sites, with optional maximum physical separation. | Earlier baseline for spatially constrained co-methylation regions. It demonstrates that adjacency plus correlation is established. |
| [Hvitfeldt et al., Methcon5, JCO Clinical Cancer Informatics (2020)](https://doi.org/10.1200/CCI.19.00109) | Manhattan distances between samples across CpGs in a gene/region, with resampling to assess conservation and account for region-size differences. | Useful for how to normalize and compare regional conservation; its comparison unit differs from two sites across patients. |
| [Zhang et al., Genome Biology (2015)](https://doi.org/10.1186/s13059-015-0581-9) | Neighboring methylation, genomic distance, CpG-island/coding annotations and regulatory features including DNase accessibility, TF binding and histone marks. | Precedent for combining several factors to predict methylation. This is a supervised prediction model, not a new pairwise correlation. |
| [Hay et al., Nature Communications (2023)](https://doi.org/10.1038/s41467-023-40845-2) | A neighboring-site similarity score: fraction of 14 clonal lines in which the closest CpG differs by at most 10 percentage points; analyzed alongside local CpG density and methylation fidelity. | Close to the agreement-within-one-level idea, but uses a raw-beta tolerance and mouse clonal lines. Binned within-one agreement is not exactly the same tolerance rule. |

## Three distances that must remain distinguishable

1. **Methylation Manhattan distance:** for the same patients at CpGs i and j,
   `D_beta = mean_k |beta_ik-beta_jk|`, or the corresponding ten-level distance.
   It describes absolute agreement. The current `S=1-D_level/9` is a standard
   normalized distance similarity.
2. **Genomic distance:** `g_ij=|position_i-position_j|`, one base-pair separation
   for the site pair. It describes proximity and is constant across patients.
3. **Statistical distance correlation:** a dependence statistic computed from
   patient-to-patient distance matrices for the two methylation vectors. It
   detects more than monotone relationships. It does not automatically include
   genomic coordinates.

For one CpG pair there are 18 paired observations but only one genomic gap.
Adding that same gap to every patient's measurement does not add an observable
within-pair covariate. Its proper role is a prior, a candidate-pair restriction,
a context feature across pairs, or a component of a prioritization score.

## A concrete combination to test, not an established new method

A transparent spatial affinity could be written

\[
A_{ij}=\exp(-g_{ij}/\ell)\,
\{\alpha S^{raw}_{ij}+(1-\alpha)C_{ij}\}\,Q_{ij},
\]

where `S_raw=1-D_beta`, `C` is a bounded dependence summary, and `Q` describes
the reliability/support available for this pair. This equation is **our
illustrative candidate**, not a formula attributed to any paper above. Its
weights, length scale and quality definition need a stated endpoint and
development/validation separation. No values were fitted in this literature
review. It is a composite affinity, not automatically a correlation, calibrated
p-value, or probability of co-regulation.

An additional computational pilot is implemented separately in
[`Item5.SpatialAffinity.Pilot.R`](../Scripts/Item5.SpatialAffinity.Pilot.R).
It fixes alpha=0.5, C=ordinary raw-beta distance correlation, Q=1, and ell=100 bp,
with 50/200-bp sensitivity analyses. Constants retain an undefined dependence
component. Its protocol is recorded before the new outcomes are computed;
this is an exploratory follow-up on an already examined cohort, not a
preregistered or independent validation. The pilot evaluates ranking for
absolute agreement in a held-out patient, which is one specific endpoint and
cannot validate all forms of dependence or biological co-regulation.

There is a useful mathematical constraint: multiplying a pair statistic by
any fixed positive factor depending only on gap/context does not change that
pair's patient-permutation ordering. For constant w>0,
`w*T_perm >= w*T_observed` iff `T_perm >= T_observed`. Thus genomic weighting
can change **which pairs rank highly relative to other pairs**, but cannot by
itself improve a fixed pair's permutation evidence. The same issue defeated
the C* rescaling argument in Task 42. Spatial information must contribute to
a regional endpoint or a validated cross-pair model to offer an actual gain.

A second, potentially more informative endpoint is **excess similarity relative
to genomic context**: ask whether a pair is more alike than comparable pairs
at similar gap, baseline methylation and dispersion. That requires an estimated
background on separate genomic blocks. Its result is a context-relative score,
not automatically an independence p-value. It also has statistical prior art;
novelty would need a narrower demonstrated improvement.

## Factors worth evaluating and what is currently available

| Factor | Why it could matter | Availability/qualification |
|---|---|---|
| Base-pair gap | Local co-methylation decays with distance | Available from input positions. Main complete consecutive pairs span 2–205 bp, so they cannot estimate a general kilobase-scale decay. |
| CpG density and intervening CpGs | Distinguishes sparse gaps from dense CpG neighborhoods | Input-row density can be measured. Reference-genome density requires a complete, build-matched CpG map; the existing sequence cache contains only short windows. |
| CpG island/shore and promoter context | Similarity can depend on genomic function | A local chr22 island file exists, but its genome-build provenance needs verification before joining to hg38. |
| Baseline level and dispersion | Two low/constant sites can agree without showing patient-specific dependence | Available, but state labels and similar means are not biological ground truth. |
| Paired sample support, ties and missingness | Prevents confusing high scores with strong evidence | Available. One shared rare patient can give maximal correlation with exact alignment probability 1/18. |
| Sequencing coverage and conversion quality | Determines measurement uncertainty | Not present in the supplied ratio-only matrix; must not be reconstructed by guessing denominators. |
| Cell composition, batch, genotype | Can create shared patient patterns without direct coupling | Needs sample metadata or compatible measurements; genomic distance does not remove these effects. |
| Regulatory marks | Could improve biologically informed regional prediction | Needs tissue- and build-matched annotations and a validation endpoint. |

## Implications for the ongoing candidate work

The earlier `correlation for cg sites` work contributes useful failures to avoid:
shared imputation can generate dependence; rare-event calibration requires
care; asymmetric quadratic cross-terms matter for U-shaped links; a combined
test must include its own selection/multiple-component calibration. These are
lessons, not inventions. Its TCGA tumor-purity model and frozen outlier protocol
are not silently applied to this RRBS cohort.

For the next candidate, keep three comparisons explicit: Manhattan agreement
alone, dependence alone (including ordinary raw-beta xi, which Task 43 omitted),
and the proposed contextual combination. Evaluate a specified regional or
held-out predictive endpoint before claiming the additional factors help.
Compare against CoMeBack/SMART's design principles and against simpler models
using exactly the same data. This is a credible research direction; adding
several established quantities together does not on its own establish novelty.

## Completed pilot

The fixed combination has now been evaluated on our complete adjacent pairs.
For ranking pairs by agreement in a held-out patient, mean absolute beta
difference was 0.01635 with Manhattan alone and 0.03990 with the 100-bp weighted
blend (lower is better). These are descriptive internal results, not a
biological validation. The simple combination is not supported as an
improvement for this endpoint. See the [full pilot report](Item5_Spatial_Pilot.md)
for the candidate formula, all component comparisons, eligibility, limitations
and reproducible outputs.

## Explicit formula developed next

The follow-up [formula derivation and experiment](Item5_Formula.md) implements
`F=S*[exp(-g/ell)+(1-exp(-g/ell))*E]`, where E is positive excess squared dCor
above its exact conditional permutation baseline. It preserves agreement as
an upper bound and changes the dependence contribution with genomic distance.
At the initial 100-bp scale, held-patient disagreement was 0.01986: better than
the earlier blend (0.03990), but worse than Manhattan alone (0.01635).
The component mathematics is established; this is an evaluated candidate
adaptation, with no novelty or independent-validation claim.
