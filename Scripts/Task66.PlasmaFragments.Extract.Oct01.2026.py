#!/usr/bin/env python3
# ================================================================
# Task 66 pilot, step 3a - rebuild whole plasma fragments for one donor
#           from Bismark paired-end alignments, two ways
# Date: Oct 1, 2026
#
# Codex's access note (Results/Task66_PlasmaFragmentRecovery_Access.md)
# found that the mHap files keep a fragment's two mates as separate
# records when their CpG ranges do not overlap (mHapTools v1.1). This
# script rebuilds the fragments of one donor from its raw reads and
# writes them two ways, so the only difference is mate linkage:
#
#   linked    one record per fragment: the calls of both mates, always
#             joined (uncalled CpGs in the gap between mates are '.')
#   mhaprule  the mHapTools rule: mates are joined only when their CpG
#             ranges intersect; otherwise each mate is its own record
#
# RULES (fixed before any output was seen)
#   Pairs: consecutive mates with the same name in Bismark's deduplicated
#     paired-end BAM, both on the same autosome; unmapped, secondary,
#     supplementary, QC-failed and duplicate records are skipped.
#   Calls: Bismark's XM tag, CpG context only (z = unmethylated,
#     Z = methylated). Reference positions come from the CIGAR. A
#     top-strand (XG:CT) call at 0-based position p is the CpG whose C is
#     at 1-based p + 1; a bottom-strand (XG:GA) call at p is the CpG whose
#     C is at 1-based p. The call must land on the hg19 CpG list.
#   Overlap: a CpG called by both mates keeps its call if they agree and
#     becomes missing ('.') if they disagree, in both outputs, so linkage
#     is the only difference between them. As in mHapTools v1.1, whether
#     mates merge is decided on their raw CpG ranges, before conflicts
#     are blanked (mHapTools then keeps the higher-quality call).
#
# OUTPUTS (prefix given on the command line; under Data/, gitignored)
#   <prefix>_linked.pat.gz, <prefix>_mhaprule.pat.gz
#     chr, global CpG index of the first CpG, pattern (C/T/.), count
#   <prefix>_fragment_stats.json
#
# Run: .venv-shared-fragment/bin/python Scripts/Task66.PlasmaFragments.Extract.Oct01.2026.py <bam> <prefix>
# ================================================================
import collections, gzip, json, re, sys, time
import numpy as np
import pysam

bam_path, prefix = sys.argv[1], sys.argv[2]
AUTO = ["chr%d" % k for k in range(1, 23)]
POS, OFF, total = {}, {}, 0
for c in AUTO:
    POS[c] = np.fromfile("Data/hg19_seq/cpg_int32/%s.int32" % c, dtype="<i4")
    OFF[c] = total; total += len(POS[c])
CALL = re.compile(r"[zZ]")


def mate_calls(r):
    """{local CpG index (0-based): '1' or '0'} for one mate."""
    xm = r.get_tag("XM")
    hits = [(m.start(), m.group()) for m in CALL.finditer(xm)]
    if not hits:
        return {}, 0
    shift = 1 if r.get_tag("XG") == "CT" else 0
    cig = r.cigartuples
    if len(cig) == 1 and cig[0][0] == 0:                 # plain match: offset arithmetic
        ref = [r.reference_start + q for q, _ in hits]
    else:                                                # indels or clipping: walk the CIGAR
        q2r = dict(r.get_aligned_pairs(matches_only=True))
        ref = [q2r.get(q) for q, _ in hits]
    p = POS[r.reference_name]; out = {}; bad = 0
    xs = np.array([x + shift if x is not None else -1 for x in ref])
    ix = np.searchsorted(p, xs)
    for (q, ch), x, i in zip(hits, xs, ix):
        if x >= 0 and i < len(p) and p[i] == x:
            out[int(i)] = "1" if ch == "Z" else "0"
        else:
            bad += 1
    return out, bad


def record(calls, lo, hi):
    return "".join({"1": "C", "0": "T"}.get(calls.get(i), ".") for i in range(lo, hi + 1))


st = collections.Counter(); ins = collections.defaultdict(collections.Counter)
linked, mhap = collections.Counter(), collections.Counter()
t0 = time.time()
with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam:
    prev = None
    for r in bam.fetch(until_eof=True):
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_qcfail or r.is_duplicate:
            st["records skipped (flags)"] += 1; continue
        if prev is None or prev.query_name != r.query_name:
            if prev is not None:
                st["unpaired records"] += 1
            prev = r; continue
        a, b = prev, r; prev = None
        st["pairs"] += 1
        if a.reference_name != b.reference_name or a.reference_name not in POS:
            st["pairs off autosomes or split"] += 1; continue
        ca, bad_a = mate_calls(a) if a.has_tag("XM") else ({}, 0)
        cb, bad_b = mate_calls(b) if b.has_tag("XM") else ({}, 0)
        st["calls off the CpG list"] += bad_a + bad_b
        raw_a, raw_b = (min(ca), max(ca)) if ca else None, (min(cb), max(cb)) if cb else None   # mHapTools tests raw ranges
        conflict = [i for i in set(ca) & set(cb) if ca[i] != cb[i]]
        st["conflicting overlap calls"] += len(conflict)
        for i in conflict:
            ca[i] = cb[i] = "."
        both = dict(ca); both.update(cb)
        both = {i: v for i, v in both.items() if v != "."}
        if not both:
            st["pairs with no CpG call"] += 1; continue
        chrom = a.reference_name; isize = abs(a.template_length); bin_ = min(isize // 50 * 50, 500)
        lo, hi = min(both), max(both)
        linked[(chrom, OFF[chrom] + lo + 1, record(both, lo, hi))] += 1
        ka = {i: v for i, v in ca.items() if v != "."}; kb = {i: v for i, v in cb.items() if v != "."}
        if ka and kb:
            la, ha, lb, hb = min(ka), max(ka), min(kb), max(kb)
            if raw_a[0] <= raw_b[1] and raw_b[0] <= raw_a[1]:
                kind = "mates' CpG ranges intersect (merged)"
                mhap[(chrom, OFF[chrom] + lo + 1, record(both, lo, hi))] += 1
            else:
                kind = "mates' CpG ranges disjoint (split)"
                mhap[(chrom, OFF[chrom] + la + 1, record(ka, la, ha))] += 1
                mhap[(chrom, OFF[chrom] + lb + 1, record(kb, lb, hb))] += 1
        else:
            kind = "only one mate has calls"
            mhap[(chrom, OFF[chrom] + lo + 1, record(both, lo, hi))] += 1
        st[kind] += 1; ins[kind][bin_] += 1
        if st["pairs"] % 500000 == 0:
            print("%d pairs (%.1f min)" % (st["pairs"], (time.time() - t0) / 60), flush=True)


def write(cnt, path):
    with gzip.open(path, "wt") as h:
        for (c, i, pat), n in sorted(cnt.items(), key=lambda kv: (AUTO.index(kv[0][0]), kv[0][1])):
            h.write("%s\t%d\t%s\t%d\n" % (c, i, pat, n))


write(linked, prefix + "_linked.pat.gz"); write(mhap, prefix + "_mhaprule.pat.gz")
out = {"bam": bam_path, "counts": dict(st), "insert_size_bins_50bp": {k: dict(sorted(v.items())) for k, v in ins.items()},
       "linked_records": sum(linked.values()), "mhaprule_records": sum(mhap.values()),
       "minutes": round((time.time() - t0) / 60, 1)}
json.dump(out, open(prefix + "_fragment_stats.json", "w"), indent=1)
print(json.dumps(out["counts"], indent=1)); print("done")
