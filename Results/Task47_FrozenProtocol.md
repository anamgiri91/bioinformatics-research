# Task 47: the frozen score and test plan

Fixed on 2026-09-29, before any new data were scored with a similarity
score. Every choice below was made from the Task 46 results only. The
code for the score is `Scripts/Task47.FrozenScore.Sep29.2026.R`. The
file hashes are in `Results/Task47_FrozenProtocol.json`.

## Why this file exists

In Task 46 many variants were tried on the same 18 RRBS patients. Some
were chosen after seeing results. A result found that way can be luck.
The fix is to write the score and the tests down first, then run them
on data the score has never seen. If a test fails, we report it.

## 1. The frozen score

    F+ = S x [w + (1 - w) x E+]

- **S, agreement.** 1 minus the mean absolute difference of the two
  sites' beta values over the reference samples.
- **E, co-methylation.** Codex's term: squared distance correlation
  above its exact chance level, scaled to 0 to 1.
- **E+, the sign guard.** E counts only when Spearman is above 0.
  Distance correlation cannot see the sign, so an inverse pair would
  otherwise get full credit.
- **w, distance.** exp(-gap / 200 bp). 200 bp is the decay length of
  neighbour Spearman fitted in Task 46 (197 bp).

## 2. The reliability flags

These are reported next to the score. They are never multiplied into
it, because in Task 46 that did not help.

| Flag | Sequencing data (RRBS, WGBS) | 450K array |
|---|---|---|
| Not unique | Bismap k50 C-to-T uniqueness track | Zhou mask `M_nonuniq` |
| Poor mapping | Bismap k50 multi-read mappability below 1 | Zhou mask `M_mapping` |
| SINE | RepeatMasker class SINE (hg38) | the same |
| General array mask | not used | Zhou `M_general` (reported only) |

A pair is flagged when either site is flagged.

## 2b. Test design shared by sections 3 and 5

These are the Task 46 rules, kept the same:

- Scores use only the training patients. Two outcomes are measured on
  each held-out patient:
  - agreement error, |x - y|, lower is better;
  - prediction gain, how much a straight line fitted on the training
    patients improves the guess of each site from the other, higher is
    better.
- Each score is judged on its top 10% of pairs. Ties are split evenly.
- **Formula families.** Each is traced over the length scales 0, 3, 5,
  10, 20, 30, 50, 75, 100, 150, 200, 300, 500, 1,000, 3,000, 10,000,
  100,000 bp and infinity:
  - Codex F;
  - F with one weight for every pair (no distance);
  - F+;
  - F guarded (sign guard, plus the Task 46 spread shrinkage m / (m + 3));
  - the planned New C (asin(sqrt) S, shrunk Spearman, distance, the
    constant rule);
  - New C with dCor in place of shrunk Spearman.
- **Same agreement.** For each comparison, take six evenly spaced error
  levels across the middle 80% of the range where the two curves
  overlap (fold means). Count the fold-by-level cells where the first
  formula has more gain.
- **Two halves.** A result has to hold in the whole genome and in each
  half. For the array the halves are the odd and even autosomes. For
  sequencing data they are the two halves of each chromosome by
  position.

## 3. Step 2: does Task 46 hold on new data? (TCGA-BRCA 450K)

- **Data.** Primary: the 53 "Alive" normal breast samples. Second: the
  32 "Dead" normal samples. Neither has been used to build or tune a
  similarity score.
- **Pairs.** The 380,355 filtered probes on autosomes, placed at their
  hg38 CpG position from the Zhou HM450 manifest (CpG_beg + 1). Each
  probe is paired with the next probe on the same chromosome.
- **Folds.** 10 patient folds, seed 20260929.

| Hypothesis, from Task 46 | Task 46 result | It holds if |
|---|---|---|
| H1 the distance weight helps F | F won 66 of 85 | F beats F-no-distance in more than half the cells |
| H2 F beats the planned New C | New C won 23 of 81 | New C wins less than half |
| H3 dCor beats shrunk Spearman inside New C | shrunk won 8 of 77 | shrunk wins less than half |
| H4 the sign guard is free | at most 0.00001 | F+ and F differ by less than 0.0005 in both outcomes at every length scale |
| H5 the spread guard costs gain | guarded won 12 of 85 | guarded wins less than half |
| H6 flagged pairs behave worse | worse at matched spread | for "not unique", "poor mapping" and SINE, flagged pairs have higher error and lower gain at matched spread (20 bins), and are worse in at least 12 of 20 bins on both |

Each hypothesis is reported for both cohorts. The Alive cohort is the
primary one.

## 4. Step 3: do neighbour-supported outliers replicate better?

This is the professor's item 1(2): can neighbours help separate
biological outliers from technical ones?

- **Flags.** The Task 13 external-reference OutlierMeth flags for the 53
  Alive normal samples, on autosomal probes whose flags are not missing.
- **Primary set.** Flags with |beta - cohort median| >= 0.10, the
  project's recommended floor, taken over the 53 normals. Flags below
  the floor are already shown to be threshold artefacts.
- **Neighbours.** Probes with flags are taken in hg38 order. The
  neighbours of a probe are the previous and next such probe on the same
  chromosome. A neighbour counts as flagged if OutlierMeth flagged it in
  the same sample and direction. The floor applies only to the flag
  being judged, not to the neighbour.
- **Similarity.** F+ for each neighbour pair, computed on the 32 Dead
  normal samples only. These are other people, so the flagged person
  cannot inflate it.
- **Outcome.** A normal flag replicates if the same patient's tumour
  sample has an OutlierMeth flag at the same probe in the same direction
  (Task 13 tumour flags). The chance level for that flag is the share of
  the other 52 patients' tumours with a same-direction flag at that
  probe. Excess replication is replicated (1 or 0) minus the chance
  level.

| Rule | Selects a flag when |
|---|---|
| Run rule (the `epimutacions` definition, the established comparison) | it is in a run of 3 or more consecutive flagged probes in the same sample and direction, each step 1,000 bp or less |
| Adjacency rule | at least one neighbour is flagged, 1,000 bp or less away |
| Similarity rule (ours) | at least one flagged neighbour has F+ >= 0.5, at any distance |

**Similarity support** is the highest F+ among the flagged neighbours,
or 0 if none are flagged. It is used for ranking.

- **Primary test.** Let K be the number of flags the run rule selects.
  Compare the mean excess replication of the top K flags by similarity
  support (ties split evenly) with that of the run-rule flags.
  - Uncertainty: 1,000 bootstrap resamples of patients, seed 20260929,
    with a 95% percentile interval.
  - **Success:** the whole interval is above 0.
- **Secondary.** These are reported, but the claim does not rest on them:
  - the three rules at their own sizes;
  - the top K by size alone (|beta - median|), to see whether support
    adds anything beyond size;
  - all flags without the floor;
  - excluding Zhou `M_general` probes;
  - tumour flags replicated in the matched normal;
  - similarity cut-offs of 0.3 and 0.7.
- **Caveat.** Replication in the matched tumour detects outliers that
  belong to the person, such as constitutive changes. A real outlier
  found only in the normal tissue, for example from cell mix, would not
  replicate. So a flag that fails to replicate is not proven technical.

## 5. Step 2, second part: a public sequencing cohort with read counts

- **Choosing the data.** The dataset is chosen by these rules, before
  any of its methylation values are looked at:
  - human RRBS or WGBS;
  - at least 15 samples of one tissue or cell type from different
    people;
  - methylated and total read counts per CpG;
  - open access.

  hg38 is preferred, and so is a source that also has read-level files
  for section 6.
- **Sites.** Keep CpGs covered by at least 10 reads in every sample. If
  fewer than 10,000 truly consecutive pairs remain, use 5.
- **Tests.** Repeat H1 to H6 exactly as in section 3, with the sequencing
  flags. Then test a new one:
  - **H7.** Weighting each sample by the smaller of the two sites' read
    counts, inside S, gives more gain at the same agreement than F+
    does. It holds if the weighted version wins more than half the cells
    in both halves.

## 6. Step 4: do close CpGs share reads?

- **Data.** Read-level files (BAM, or `wgbstools` pat files) for at
  least one sample from the section 5 source. If that source has none, use
  another open source of the same tissue.
- **Linkage.** Take each CpG pair within 200 bp with at least 10 reads
  covering both sites. Linkage is the observed share of those reads with
  the same state at both sites, minus the share expected if the sites
  were independent: p_i p_j + (1 - p_i)(1 - p_j), with p from the same
  reads.
- **H8.** Mean linkage is above 0, and it falls with distance.
- **H9, the link to Task 46.** Look only at the low-spread pairs (one
  site has no sample more than 0.1 from its median; these are the pairs
  the spread guard penalised). Among them, pairs with above-median
  linkage have higher held-out prediction gain in the section 5 cohort than
  pairs below the median, in both halves. This needs the section 5 cohort
  and the read-level data to cover the same CpGs.

## 7. Rules for changes

- Nothing above changes after new results are seen. That covers the
  score, the flags, the data choices, the outcomes and the success
  rules.
- If something must change, for example because of a bug, it is
  reported as a deviation with the reason. Results are then shown both
  ways where possible.
- Extra analyses are allowed, but they are labelled exploratory.
