# Task 51: a bug in the "not unique" flag, found after the first run

Recorded 2026-09-29. This fix came after a result was seen, so it is
reported here as a change, not hidden.

## What went wrong

The first complete run of Task 51 marked 1,439,615 of 1,574,581 pairs (91%)
as "not unique". On the development RRBS data the same Bismap track flagged
only 4.7% of pairs. By chromosome, the share was 3.8% on chr1 and exactly
1.000 on every other chromosome.

The cause is in how the track was read. rtracklayer's bigBed import was given
query ranges across all chromosomes at once, and it returned features for the
first chromosome only. We checked this directly: a query over chr1, chr2,
chr21 and chr22 came back with chr1 features only. So every site outside chr1
looked like it was in no unique region.

- The mappability bigWig, read the same way, was not affected: its flag
  shares are 17% to 31% on every chromosome.
- The SINE flag does not use rtracklayer, so it was not affected.

## The fix

The uniqueness track is now read one chromosome at a time. The script stops
if any chromosome has fewer than half of its sites unique. On chr22 the cohort
sites are 93.8% unique, close to the development data's 95.6%.

## What it affects

- Only the "not unique" part of H6.
- H1 to H5 and H7 do not use the flags. Their results come from the same
  score, folds and seed, so they are unchanged.

The first run's reliability table is not kept.

## Historical values recovered from the handoff

The user supplied the before/after rows on September 30. The following records
are transcribed from that handoff, not recovered from an original output file.
They preserve the three affected rows for review; the other flags were unchanged.

| Subset | Initial nonunique pairs | Error difference | Gain difference | Bins worse error | Bins lower gain | Initial holds |
|---|---:|---:|---:|---:|---:|---|
| All | 1,439,615 | −0.00087 | −0.00019 | 9/20 | 13/20 | No |
| First | 717,600 | −0.00317 | +0.00039 | 3/20 | 5/20 | No |
| Second | 722,015 | +0.00247 | −0.00100 | 18/20 | 20/20 | Yes |

The corrected rows are in `Task51_Reliability_SpreadMatched.csv`. The first two
H6 aggregate verdicts changed from false to true; the second-half verdict was
already true. The score comparisons H1–H5 and H7 are unchanged.
