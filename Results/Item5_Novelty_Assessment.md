# Novelty assessment of the proposed item-5 formula

Reviewed 2026-09-29 after the user asked whether the formula is novel.

**Assessment: methodological novelty is not established. Do not describe this
as a novel correlation coefficient.** The strongest currently defensible
description is a proposed spatially informed composite similarity score,
assembled from established components and evaluated internally on this RRBS
matrix. We have not located an exact published CpG implementation in the
targeted search; that does not establish originality.

## A closer mathematical equivalence

Our formula can be rewritten exactly as

\[
F=S\{w+(1-w)E\}=S(w+E-wE)=T\{S,P(w,E)\},
\]

where T(a,b)=ab is the product operator and P(a,b)=a+b-ab is the standard
probabilistic-sum operator. These are established fuzzy conjunction/disjunction
constructions, explicitly used in the primary paper
[Generalized fuzzy Petri nets as pattern classifiers](https://www.sciencedirect.com/science/article/abs/pii/S0167865599000732).
They also appear explicitly in the more recent primary research article
[Counterfactuals in fuzzy relational models](https://doi.org/10.1007/s10462-024-10996-9).
The identification of our gate with those operators is an algebraic deduction,
not a claim that either paper studies CpG methylation or our exact inputs.

This gives a precise interpretation: agreement AND (proximity OR dependence),
under these particular fuzzy operators. It does not turn S, w, E or F into
calibrated probabilities. The label "probabilistic sum" names the operator.

## What is already established

| Component | Assessment |
|---|---|
| Normalized Manhattan agreement S | Standard distance-to-similarity transformation; methylation Manhattan precedents include MPCI and Methcon5. |
| Squared sample distance correlation q | Established dependence statistic. |
| Conditional permutation expectation q0 | Established RV/centered-matrix expectation; its application to the distance matrices is explained in the formula derivation. |
| Subtracting a chance baseline and dividing by 1 minus that baseline | An established chance-adjustment construction; positive-part clipping does not establish a new dependence estimator. |
| Exponential genomic-proximity weight | Established spatial-decay construction, including methylation covariance modeling. |
| Combining S, w and E through product/probabilistic-sum operators | Standard aggregation structure, as the algebraic equivalence above shows. |
| This exact choice of inputs and application to the 18-patient RRBS matrix | A proposed adaptation; exact publication priority has not been established. |

Primary references for the methylation and dependence components are in the
[literature review](Item5_Manhattan_Spatial_PriorArt.md) and
[formula derivation](Item5_Formula.md). The component precedents and the
aggregation identity weaken a claim of new mathematical machinery, while
leaving open the possibility of an application-specific contribution.

## Novelty and usefulness are separate questions

A method can be original without outperforming a baseline, and a useful
implementation can use established mathematics. Our present experiment gives
no advantage over Manhattan alone for held-patient absolute agreement:
0.01986 versus 0.01635 mean absolute beta difference (lower is better).
The candidate improved on our earlier blend, which is insufficient to
establish an advantage over published or simple established methods.

The next defensible contribution would need a specific unresolved target,
such as regional recovery under sparse observations, a justified treatment
of measurement uncertainty when read counts exist, or a demonstrated
finite-sample property. It would also need comparison with appropriate
established methods and separate validation. Neither adding more factors nor
finding no exact equation in a search is enough by itself.

**Current claim wording:** "We propose and evaluate a spatially informed
composite similarity score combining methylation agreement and dependence."
Avoid "first," "novel correlation coefficient," and "superior method" on the
current evidence.
