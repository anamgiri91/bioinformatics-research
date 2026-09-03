"""Builds Meeting_2026-09-03_Outlier_Methylation.docx from the Results/ CSVs.

Answers the 2026-09-03 pre-meeting email: the two literature questions, the
two analysis questions, and the request for new ideas on improving the
flagging. Regenerate after any analysis rerun so the document cannot drift
from the numbers it cites.

    python3 -m venv .venv && .venv/bin/pip install python-docx
    .venv/bin/python Scripts/build_meeting_doc_sep03.py

Figure numbering in the document is by order of APPEARANCE and does not match
the Fig<N>_*.png filenames, which are numbered by order of creation:
    Figure 1 -> Fig9_ContraryConcentration.png
    Figure 2 -> Fig11_StateFloors.png
    Figure 3 -> Fig8_PEnvelope.png
    Figure 4 -> Fig12_MethodScores.png
    Figure 5 -> Fig10_TransformAmplification.png

Numbers are read from the CSVs where a small helper makes that clean and are
otherwise quoted from results.md, which is itself generated from them. Any
number appearing here also appears in results.md sections 10-13.
"""

import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from docx import Document
from docx_kit import Kit, base_styles, NAVY, GREY, RED, CENTER

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(REPO, "Results")
OUT = os.path.join(REPO, "Meeting_2026-09-03_Outlier_Methylation.docx")


def load(name):
    with open(os.path.join(RES, name)) as fh:
        return list(csv.DictReader(fh))


def fmt_p(v):
    """0.0001 rather than 1e-04, which is what the package documentation says."""
    return ("%.5f" % float(v)).rstrip("0")


def sig(v):
    """Five decimals throughout, so the columns line up as numbers."""
    return "%.5f" % float(v)


def pick(rows, **kw):
    """First row matching every key=value pair, as strings."""
    for r in rows:
        if all(str(r.get(k, "")) == str(v) for k, v in kw.items()):
            return r
    raise KeyError(kw)


doc = Document()
base_styles(doc)
k = Kit(doc, RES)

# ============================ TITLE ============================
p = k.para("Outlier DNA Methylation in TCGA-BRCA", size=21, bold=True,
           align=CENTER, space_after=2)
p.runs[0].font.color.rgb = NAVY
k.para("Progress report for the meeting of 3 September 2026", size=12,
       align=CENTER, space_after=2).runs[0].font.color.rgb = GREY
k.para("Anam Giri", size=10.5, align=CENTER, space_after=16)\
    .runs[0].font.color.rgb = GREY

k.callout("What this report is",
 "It answers the four questions in the pre-meeting email — where the contrary "
 "flags come from, whether a state-pooled threshold works, whether the "
 "OutlierMeth parameters can be tuned, and what the Borealis and epimutacions "
 "packages do — and then reports what came out of answering them: a correction "
 "to how the external reference has been described, a measured replacement for "
 "the one assumption the reproducibility argument rested on, and five candidate "
 "fixes scored against each other.\n\n"
 "Every number is reproducible from a named script and CSV in the project "
 "repository, listed in the appendix. Where a claim is not supported, it is "
 "labelled as such rather than softened.", 'E8EEF7')

# ============================ SUMMARY ============================
doc.add_heading("Summary", 1)

k.para("**The four questions, answered.**")
k.table(["#", "Question asked", "Answer"],
 [["1", "For the 100 normal chr22 sites: are the contrary flags from all CG "
        "sites or 1–2? From 1–3 samples, or can every sample have a chance?",
   "**Concentrated on sites, spread across samples.** Half the sites of a state "
   "carry all of its contrary flags and the top two hold 57–86%. But 31 of 53 "
   "normal samples carry at least one (§1)"],
  ["2", "Try a “bio-stat” flag: the >96th percentile of the whole state pool, "
        "8×53 for the L sites",
   "**It works, with one change.** Pool the deviation from each site's median, "
   "not raw beta, and use the whole chromosome rather than 424 cells. It removes "
   "the n = 53 degeneracy and makes the flag rate exactly state-independent (§2)"],
  ["3", "Can the “parameter” setting be adjusted to make the external "
        "reference less sensitive?",
   "**No.** The p argument has two usable settings on this reference, not four, "
   "and even the strictest leaves 100% of the H-site flags below a beta shift of "
   "0.10 (§3)"],
  ["4", "Does OutlierMeth report or identify outlying samples? And check "
        "Borealis and epimutacions",
   "**It does not** — it returns a CpG × sample matrix and nothing else. Both "
   "other packages carry the missing piece: an explicit magnitude requirement (§4)"]],
 widths=[0.3, 3.0, 3.3])
k.source("Full detail: results.md §§10–13 and LITERATURE.md in the repository.")

k.callout("One correction before anything else",
 "Previous reports described the external reference as **2,015 samples across 25 "
 "tissue types**. That is the combined TCGA-GEO panel. This project loads "
 "tcga.rda and calls flagMeth(…, reference = tcga), which is **747 TCGA normal "
 "samples across 21 tissue types**.\n\n"
 "The code was never ambiguous; the write-up was wrong, and it has been corrected "
 "throughout. It matters, and it makes the external arm harder to defend rather "
 "than easier: at n = 747, three of the four significance levels the package "
 "offers collapse onto the reference maximum (§3).", 'FDF0F0')

k.para("**What changed in the conclusions.** Nothing already reported is "
       "overturned. Three things are now settled that were open:")
k.bullet("The reproducibility argument no longer rests on a literature value. "
         "The noise level is measured from this cohort's own adjacent probes, "
         "and the assumed 0.01 sits inside the measured range.")
k.bullet("Tuning the package's parameters is ruled out as an alternative to "
         "adding a magnitude floor — not argued, computed.")
k.bullet("The flat 0.10 floor survived a head-to-head against four "
         "elaborations of it, including the 0.05 proposed in the email.")

k.pagebreak()

# ============================ 1 ============================
doc.add_heading("1. Where the contrary flags come from", 1)

k.para("A **contrary flag** is one that runs against the state the CpG already "
       "sits in: a “too low” call at a site that is unmethylated in everyone, or "
       "a “too high” call at a site that is fully methylated in everyone. The "
       "email listed the four such cells in the 100-CpG normal window and asked "
       "where they come from.")

k.para("The counts reproduce exactly, with one difference: the H cell is **10 "
       "flags across 11 H sites**, not 11 flags.")

d29 = load("Task29_Window100_ContraryDecomposition.csv")
rows = []
for st, dr in (("L", "-1"), ("LM", "-1"), ("H", "+1"), ("HM", "+1")):
    r = pick(d29, tissue="Normal", **{"methy.state": st, "direction": dr})
    rows.append([st, dr, r["flags"], r["sites.in.state"],
                 "**%s**" % r["sites.with.flag"], r["max.per.site"],
                 "**%s%%**" % r["top2.site.share.pct"],
                 "%s / 53" % r["samples.with.flag"], r["max.per.sample"],
                 "%s%%" % r["top3.sample.share.pct"]])
k.table(["state", "dir", "flags", "sites in state", "sites with ≥1",
         "most at one site", "top-2 sites hold", "samples with ≥1",
         "most in one sample", "top-3 samples hold"], rows,
        widths=[0.42, 0.32, 0.42, 0.62, 0.62, 0.62, 0.66, 0.68, 0.62, 0.66],
        font=8.5)
k.source("Results/Task29_Window100_ContraryDecomposition.csv · "
         "Scripts/Task29.ContraryFlagDecomposition.Sep03.2026.R")

k.para("**Sites: concentrated.** Half the sites of a state carry every contrary "
       "flag it has. A Monte Carlo test against uniform allocation rejects at "
       "p = 0.0008 for the normal LM cell and below 0.0001 for the tumour one. "
       "In tumour, a single probe — cg15668074 — carries 19 of the 29 LM “too "
       "low” flags on its own.")

k.para("**Samples: spread.** No sample holds more than 1–7 flags in any one "
       "cell, and pooling all four cells:")
k.table(["tissue", "contrary flags", "samples carrying ≥1", "top 3 hold",
         "uniform expectation", "leaders"],
 [["Normal", "85", "**31 / 53**", "27.1%", "5.7%", "N15 = 9, N48 = 9, N17 = 5"],
  ["Tumour", "76", "**37 / 53**", "19.7%", "5.7%", "T33 = 6, T36 = 5, T14 = 4"]],
 widths=[0.7, 0.9, 1.1, 0.8, 1.1, 1.9])
k.source("Results/Task29_Window100_ContraryBySample.csv")

k.para("So a **majority** of samples can and do carry a contrary flag. The lead "
       "N15 and N17 hold is real but modest. N48 is new — it was not in the "
       "N15 / N17 / N14 group identified earlier — and a check against the "
       "chromosome-wide burden table confirms it belongs there: N48 ranks **4th "
       "of 53** on chr22, so this is a fourth high-burden sample rather than a "
       "quirk of the window.")

k.figure("Fig9_ContraryConcentration.png",
 "**Figure 1.** The two halves of the question have opposite answers. **Left:** "
 "flags per site, sites ranked, one line per cell — only sites that carry at "
 "least one flag are drawn, and the legend gives how many of the state's sites "
 "that is. **Right:** the same flags counted per sample, all four cells pooled.")

k.callout("Why those sites, then? Not biology.",
 "The external threshold knows nothing about this cohort, so the number of "
 "samples a site flags is decided by how far that threshold happens to sit from "
 "this cohort's own distribution at that site. Recovering the thresholds from "
 "the flag matrix and measuring that distance in cohort-MAD units, the "
 "correlation with contrary-flag count is **negative in all six testable state × "
 "tissue cells** (Spearman −0.95 to −0.12), at a median distance of only "
 "**1.0–1.9 MAD**.\n\n"
 "The concentrating sites are the geometrically unlucky ones. Probe masking will "
 "remove some of them, but the pattern will reappear at the next-closest sites, "
 "so both numbers should be reported.")

k.pagebreak()

# ============================ 2 ============================
doc.add_heading("2. One threshold per state instead of one per CpG", 1)

k.para("The “bio-stat” proposal is a direct answer to a real problem. An "
       "empirical percentile over 53 values is an order statistic: at n = 53 the "
       "99th percentile falls between the 52nd and 53rd sorted values, so exactly "
       "one sample can exceed it, at every significance level. A percentile over "
       "8 × 53 = 424 cells, or over 117,183 for all the chromosome's L sites, is "
       "an actual estimate, and a site can then have zero outliers — which is "
       "what most sites should have.")

k.para("**Two ways to pool, and the difference decides whether it works.**")
k.table(["", "what is pooled", "the threshold is"],
 [["bio-stat, absolute", "raw beta values of the state",
   "one absolute beta value per state — the email's version literally"],
  ["bio-stat, deviation", "beta minus that site's own cohort median",
   "one **effect size** per state, re-centred at each site"]],
 widths=[1.3, 2.2, 3.1])

k.para("States are wide — L spans beta 0.005 to 0.10 — so under the absolute "
       "rule the top of the pooled tail is dominated by whichever **sites** sit "
       "highest within the state. It ends up flagging whole sites rather than "
       "unusual samples. The deviation version removes each site's own level "
       "first and asks the intended question. At the same 4% flag rate on tumour "
       "R sites the two give a median beta shift of **0.016** and **0.465** "
       "respectively, on the same cells.")

k.table(["tissue", "variant", "sites carrying a flag", "median beta shift",
         "under 0.10", "stability", "state-rate ratio"],
 [["Normal", "absolute", "40.9%", "0.079", "63.9%", "0.823", "1.00"],
  ["Normal", "**deviation**", "**63.5%**", "**0.137**", "**38.6%**", "0.804", "1.00"],
  ["Tumour", "absolute", "41.7%", "0.078", "60.4%", "0.779", "1.00"],
  ["Tumour", "**deviation**", "**63.5%**", "**0.242**", "**30.2%**", "**0.865**", "1.00"]],
 widths=[0.7, 0.9, 1.1, 1.0, 0.8, 0.8, 1.0], highlight_rows=(1, 3))
k.source("chr22, q = 0.96 · Results/Task30_StatePooled_Scores.csv · "
         "Scripts/Task30.StatePooledThreshold.Sep03.2026.R")

k.para("**One practical warning.** The threshold estimated from the window's 424 "
       "L cells is 0.084; from all 117,183 chr22 L cells it is 0.102 — 18% "
       "higher. The small pool biases it low. Use the chromosome, not the window.")

k.para("**What it fixes and what it does not.** It fixes the state-rate bias "
       "completely — every state gets exactly 1 − q of its cells flagged by "
       "construction — and it fixes the degeneracy. It does **not** fix the "
       "magnitude problem, because the 96th percentile of the L and H states' own "
       "deviations is only 0.037 and 0.046. A state-pooled percentile is still a "
       "percentile.")

k.callout("The number the email asked about: is 0.05 the right floor?",
 "The per-state thresholds above are exactly a data-calibrated version of the "
 "proposed |beta − median| > 0.05 rule. Read off each state's own deviation "
 "distribution at the 96th percentile (chr22, normal): **L 0.037, LM 0.142, "
 "M 0.144 / 0.179, HM 0.125, H 0.046, R 0.306 / 0.341**.\n\n"
 "So **0.05 is about right at L and H, and three to seven times too lenient "
 "everywhere else.** No single constant can be right at every state — but §5 "
 "shows that a flat 0.10 still beats the per-state version in practice.")

k.figure("Fig11_StateFloors.png",
 "**Figure 2.** The floor each state's own data picks, against the three "
 "constants under discussion. Each state gets one bar per tail it has a rule "
 "for. The R state needs six to eight times the floor that L and H need.")

k.pagebreak()

# ============================ 3 ============================
doc.add_heading("3. Can the parameters be adjusted?", 1)

k.para("OutlierMeth has exactly two settings: **reference**, one of four packaged "
       "panels, and **p**, one of four significance levels. The answer is no, on "
       "two independent grounds, and neither needed the reference file itself — "
       "which is not available on this machine.")

doc.add_heading("3.1  The p setting has less resolution than it appears", 2)

k.para("The package builds every threshold with R's ordinary quantile function. "
       "A quantile at level p is a real estimate of a tail only when the "
       "reference has roughly 1/p samples or more. The reference here has 747:")

d31 = load("Task31_PLevelResolution.csv")
rows = []
for r in d31:
    if r["panel"] != "tcga":
        continue
    rows.append([fmt_p(r["p"]), r["position.h"],
                 ("**%s**" % ("1st–2nd" if r["distinct.from.max"] == "FALSE"
                              else r["rank.from.top"] + "th")),
                 r["expected.n.above"],
                 "yes" if r["distinct.from.max"] == "TRUE" else "**no**"])
k.table(["p", "position in the sorted 747", "rank from the top",
         "reference samples expected above", "a real tail estimate?"], rows,
        widths=[0.7, 1.6, 1.2, 1.7, 1.3], highlight_rows=(1, 2, 3))
k.source("Results/Task31_PLevelResolution.csv · "
         "Scripts/Task31.ParameterHeadroom.Sep03.2026.R")

k.para("The three strictest levels all interpolate between the **two most "
       "extreme reference samples**. The two strictest differ from each other by "
       "6% of one gap between samples. Reaching p = 0.0001 honestly would need a "
       "reference of 10,002 samples; the largest panel the package ships has "
       "2,015. **So the p argument offers two usable settings on this reference, "
       "not four.**")

k.callout("This is the same problem, now on the other side",
 "The self-referential arm of this project was already known to be degenerate "
 "for exactly this reason at n = 53. What is new is that the same order-statistic "
 "limit reaches the **external** arm at the strict end. It is a property of "
 "estimating an extreme quantile from a finite sample, not a property of either "
 "dataset.")

doc.add_heading("3.2  And a stricter setting would not help anyway", 2)

k.para("Thresholds only move upward as p is tightened, so the set of flagged "
       "measurements shrinks in a nested way: any stricter setting keeps a subset "
       "of the current flags, always dropping the ones closest to the threshold "
       "first. That means every outcome reachable by tuning p can be worked out "
       "from the results already in hand, without the reference file. The "
       "strictest setting that still flags anything keeps one measurement per "
       "CpG per direction.")

d31b = load("Task31_PEnvelope.csv")


def env(tissue, state, direction, r):
    for x in d31b:
        if (x["tissue"] == tissue and x["methy.state"] == state
                and x["direction"] == direction
                and x["keep.top.per.cpg"] == r):
            return x
    raise KeyError((tissue, state, direction, r))


rows = []
for lbl, st, dr in (("H sites, “too high”", "H", "+1"),
                    ("L sites, “too low”", "L", "-1"),
                    ("R sites, “too high” (for contrast)", "R", "+1")):
    a = env("Normal", st, dr, "NA")
    b = env("Normal", st, dr, "1")
    rows.append([lbl, "p = 0.01, as run", a["flags"], a["median.abs.dbeta"],
                 "**%s%%**" % a["pct.under.0.10"],
                 "**%s**" % a["surviving.floor.0.10"]])
    rows.append(["", "strictest reachable p", b["flags"], b["median.abs.dbeta"],
                 "**%s%%**" % b["pct.under.0.10"],
                 "**%s**" % b["surviving.floor.0.10"]])
k.table(["cell", "setting", "flags", "median beta shift", "under 0.10",
         "surviving a 0.10 floor"], rows,
        widths=[1.7, 1.3, 0.7, 1.0, 0.8, 1.1], highlight_rows=(0, 1))
k.source("chr22, normal tissue · Results/Task31_PEnvelope.csv")

k.para("**The ceiling is structural.** At H sites the median methylation across "
       "the cohort is 0.934, so there is only 0.066 of scale left above it. Just "
       "**0.24%** of all H measurements sit 0.10 or more from their site's median "
       "— against 38.4% at R sites. There is nothing there for a stricter "
       "threshold to find.")

k.figure("Fig8_PEnvelope.png",
 "**Figure 3.** For each contrary cell, the median beta shift of its flags at "
 "the setting used against the strictest any setting of p could reach. If "
 "tightening p helped, the dark bar would be taller. The dashed line is 0.10, "
 "the smallest shift worth calling real on this platform.")

k.para("**In one sentence: p changes how many flags you get, and never how big "
       "they are.**")

k.pagebreak()

# ============================ 4 ============================
doc.add_heading("4. Does OutlierMeth identify outlying samples? And the two packages", 1)

k.para("**No.** The package exports four functions and the unit of every one of "
       "them is the individual measurement, not the sample.")
k.table(["function", "what it returns"],
 [["referenceMeth(beta)", "the per-probe threshold table"],
  ["flagMeth(beta, reference, p)", "a CpG × sample matrix of −1 / 0 / +1"],
  ["deltMeth(…)", "the same cells as beta − threshold"],
  ["relMeth(…)", "the same, divided by the remaining room to the boundary"]],
 widths=[2.0, 4.3])
k.para("There is no sample-level call, no burden statistic and no "
       "multiple-testing control. Every per-sample table in this project is a "
       "column sum we compute ourselves. If “outlying sample” is a wanted output "
       "it has to be defined here, with its own null — and the batch, purity and "
       "cell-composition confounds excluded first.")

k.callout("A trap worth flagging: relMeth points the wrong way",
 "deltMeth returns beta − threshold, which is one comparison away from the "
 "magnitude floor this project recommends — the machinery is already there. "
 "relMeth divides that by the room left to the boundary, which at H sites is "
 "0.037. Measured on chr22 normal tissue, the median deltMeth at H sites is "
 "**0.0022** and the median relMeth is **0.0746**: a **34-fold inflation**, and "
 "11-fold at HM.\n\n"
 "A user reaching for relMeth to get an interpretable effect size would get the "
 "flags whose direction the state made inevitable coming back as the largest "
 "effects in the dataset.", 'FDF0F0')

doc.add_heading("4.1  Borealis", 2)
k.para("Outlier detection for **bisulfite sequencing read counts**, aimed at rare "
       "disease diagnosis. It fits a beta-binomial per CpG — methylated reads out "
       "of total reads, with a cohort mean and a dispersion — and tests each "
       "sample against that fitted distribution. Output carries a p-value, an "
       "effect size, and a direction field whose third value is “change too small "
       "to call”.")
k.para("**Why it matters here.** It is the model-based version of what this "
       "project keeps concluding is needed: dispersion is in the model, so a "
       "fixed shift is not equally significant everywhere; read depth carries its "
       "own uncertainty; and “too small to call” is a built-in class rather than "
       "something bolted on afterwards.")
k.para("**Why it cannot be used directly.** It needs read counts. A 450k beta "
       "value is an intensity ratio with no underlying count, so a beta-binomial "
       "cannot be fitted without inventing one. The transferable idea is the "
       "shape, not the code.")

doc.add_heading("4.2  epimutacions", 2)
k.para("Six outlier methods for array beta values, run either against a "
       "reference panel or leave-one-out within a cohort. An **epimutation** is "
       "defined as a run of at least 3 outlier CpGs within 1 kb — the definition "
       "this project already adopted after N37's eleven flags turned out to be "
       "three events.")
k.table(["method", "parameters and defaults"],
 [["quantile", "window_sz = 1000, qsup = 0.995, qinf = 0.005, **offset_abs = 0.15**"],
  ["beta", "pvalue_cutoff = 1e-06, **diff_threshold = 0.1**"],
  ["manova / mlm", "pvalue_cutoff = 0.05"],
  ["iForest", "outlier_score_cutoff = 0.7, ntrees = 100"],
  ["mahdist", "nsamp = \"deterministic\""]],
 widths=[1.4, 4.9])
k.callout("The single most useful thing in the literature review",
 "**Both** of the threshold-based methods ship a mandatory magnitude floor — "
 "0.15 on the quantile method, 0.10 on the beta method. Neither will report a "
 "statistically extreme but biologically trivial measurement.\n\n"
 "And note the structure: epimutacions' quantile method **is the same rule as "
 "flagMeth** — an empirical percentile of a reference panel — with one extra "
 "term. The entire difference between the two packages at that method is the "
 "absolute-difference floor. That is independent precedent for the fix this "
 "project arrived at from the data.", 'E8EEF7')

k.pagebreak()

# ============================ 5 ============================
doc.add_heading("5. New ideas, tested", 1)

k.para("With tuning ruled out, any fix has to add a term the package does not "
       "have. Five candidates were built and scored against each other at a "
       "single flag rate, so none of them can look better simply by flagging "
       "less.")

doc.add_heading("5.1  First: the noise level is now measured, not assumed", 2)

k.para("Every stability figure this project has reported perturbs the data by a "
       "random amount with a standard deviation of 0.01, because that is what the "
       "literature says about this array's technical error. That was the one "
       "assumption the whole reproducibility argument rested on. Neighbouring "
       "probes replace it with a measurement: for two probes within 100 bp, "
       "whatever systematic difference exists between the two positions is the "
       "same in every sample and cancels, so what is left bounds the technical "
       "noise from above.")

d32 = load("Task32_MeasuredNoise.csv")


def noise(tissue, pair_set):
    return next(x for x in d32 if x["tissue"] == tissue
                and x["pair.set"] == pair_set and x["methy.state"] == "ALL")


rows = []
for ps in ("<= 50 bp", "<= 100 bp", "<= 200 bp", "<= 500 bp",
           "distant (>1 Mb, null)"):
    n, t = noise("Normal", ps), noise("Tumor", ps)
    lbl = ps.replace("<=", "≤")
    rows.append([("**%s**" % lbl) if ps == "<= 100 bp" else lbl,
                 sig(n["sigma.p10"]), sig(n["sigma.median"]),
                 sig(t["sigma.p10"]), sig(t["sigma.median"])])
k.table(["probe pairs", "normal: tightest 10%", "normal: median",
         "tumour: tightest 10%", "tumour: median"], rows,
        widths=[1.7, 1.2, 1.2, 1.2, 1.2], highlight_rows=(1,))
k.source("chr22 · Results/Task32_MeasuredNoise.csv · "
         "Scripts/Task32.ImprovedFlagging.Sep03.2026.R")

k.para("Close pairs are about three times tighter than distant ones, so the "
       "bound is measuring probe-level noise rather than biology. **The assumed "
       "0.01 sits inside the measured range and is conservative relative to the "
       "median, so every previously reported stability result stands.** The "
       "caveat is that this bounds technical noise plus any genuine divergence "
       "between neighbouring probes, so it is an upper bound.")

doc.add_heading("5.2  The scores", 2)

d32s = load("Task32_ImprovedMethodScores.csv")
LABEL = {
    "ext.percentile": "external reference, as published",
    "ext.floor.med.0.05": "+ floor at 0.05 (proposed in the email)",
    "ext.floor.med.0.10": "+ floor at 0.10",
    "ext.floor.state": "+ per-state calibrated floor",
    "ext.floor.noise": "+ floor scaled by measured noise",
    "ext.delt": "+ floor on distance to threshold (deltMeth)",
    "bio.stat.dev": "bio-stat, deviation — **no reference panel at all**",
    "mad.beta": "median ± k·MAD on beta",
    "asin.mad": "median ± k·MAD on the angular transform",
    "mad.mvalue": "median ± k·MAD on M-values",
}
for tissue, disp in (("Normal", "Normal tissue"), ("Tumor", "Tumour")):
    sub = [x for x in d32s if x["tissue"] == tissue]
    sub.sort(key=lambda x: -float(x["stability.jaccard"]))
    target = next(x for x in sub if x["method"] == "ext.floor.med.0.10")["flag.rate.pct"]
    rows = []
    for x in sub:
        matched = abs(float(x["flag.rate.pct"]) - float(target)) < 0.01
        rows.append([LABEL[x["method"]],
                     x["flag.rate.pct"] + ("" if matched else " †"),
                     x["stability.jaccard"], x["median.abs.dbeta"],
                     x["pct.trivial"] + "%", "**" + x["pct.contrary"] + "%**"])
    k.para("**%s** — all rows at %s%% flag rate except those marked †."
           % (disp, target))
    k.table(["method", "rate", "stability", "median beta shift",
             "trivial flags", "contrary flags"], rows,
            widths=[2.6, 0.6, 0.8, 1.0, 0.8, 0.9], highlight_rows=(0,), font=8.5)
k.source("chr22 · Results/Task32_ImprovedMethodScores.csv. **Stability** is the "
         "overlap between the flags on clean data and on data perturbed by noise "
         "of sd 0.01, with every self-derived threshold re-estimated. **Trivial** "
         "is the share of flags below a beta shift of 0.10. **Contrary** is the "
         "share running against the site's state. † flags at its own, higher, "
         "rate — so its loss is not an artifact of sensitivity.")

k.figure("Fig12_MethodScores.png",
 "**Figure 4.** The benchmark on two axes; up and to the right is better — "
 "reproducible flags that correspond to a real change. Point size is the share "
 "of flags running against the state.")

doc.add_heading("5.3  Five conclusions", 2)
k.numbered("**Over half the external reference's normal-tissue flags run against "
           "the state.** 52.9% on chr22, 27.0% in tumour. A 0.10 floor cuts that "
           "to 21.0% and 6.5%. This is probably the clearest one-line statement "
           "of the problem the project has.")
k.numbered("**The plain constant floor still wins.** Requiring the measurement to "
           "sit at least 0.10 from the cohort median, on top of the external "
           "flag, is the most stable rule tested in both tissues. The "
           "elaborations do not beat it.")
k.numbered("**0.05 is measurably worse than 0.10, even with a rate advantage.** "
           "It scores 0.748 / 0.871 against 0.793 / 0.907 while flagging *more*, "
           "and leaves 39% / 23% of its flags below a shift of 0.10. If one "
           "constant has to be chosen, 0.10 is better; epimutacions uses 0.15.")
k.numbered("**The state-pooled rule is the best option that needs no reference "
           "panel, by a wide margin** — 0.779 / 0.833 against 0.429 / 0.678 for "
           "the standard robust alternative, a perfectly state-independent flag "
           "rate, and no contrary flags at all by construction. It needs no "
           "external panel and no assumption that a pan-tissue reference applies "
           "to breast. This is the suggestion from the email, and it is the "
           "strongest thing in the table after the constant floor.")
k.numbered("**Scaling the floor by the measured noise makes things worse — a "
           "negative result.** The noise is smallest exactly where the state is "
           "most compressed, so the floor at L sites lands at 0.023 and keeps the "
           "trivial flags. Noise scaling and interpretability are different "
           "requirements: a shift of 0.02 at an L site may be many standard "
           "deviations, and it is still not a finding.")

doc.add_heading("5.4  A question that is now settled: what about M-values?", 2)
k.para("An earlier report found M-values — the textbook correction for this "
       "array's non-constant variance — to be the **least** stable method tested, "
       "and blamed the logarithmic transform for magnifying noise at the "
       "extremes. That was a hypothesis. Adding the **angular transform**, which "
       "stabilises variance in the same way but is bounded and so magnifies far "
       "less, tests it directly.")
k.table(["tissue", "transform", "magnification where the flags are", "stability",
         "contrary flags"],
 [["Normal", "raw beta", "1.00", "**0.429**", "**0.1%**"],
  ["Normal", "angular", "1.48", "0.326", "2.5%"],
  ["Normal", "M-value", "2.26", "**0.192**", "**12.1%**"],
  ["Tumour", "raw beta", "1.00", "**0.678**", "**1.3%**"],
  ["Tumour", "angular", "1.37", "0.546", "5.8%"],
  ["Tumour", "M-value", "2.17", "**0.429**", "**14.6%**"]],
 widths=[0.8, 1.1, 2.1, 1.1, 1.1])
k.para("Monotone in both tissues, on both measures: **the more a transform "
       "stretches the extremes, the less reproducible its flags and the more of "
       "them run against the state.** The hypothesis holds. The practical "
       "conclusion is not “use a gentler transform” — the angular one is gentler "
       "and still loses to untransformed beta — but **do not transform; add a "
       "magnitude requirement instead**.")

k.figure("Fig10_TransformAmplification.png",
 "**Figure 5.** Magnification of noise at the measurements each method flags, "
 "against reproducibility (left) and against the share of flags running counter "
 "to the state (right). Both panels, both tissues, all at the same flag rate.")

k.pagebreak()

# ============================ 6 ============================
doc.add_heading("6. What to decide", 1)

k.numbered("**Which fix is the recommendation?** The external reference plus a "
           "0.10 floor is the most stable rule and keeps the published method. "
           "The state-pooled deviation rule is only slightly behind, needs no "
           "reference panel at all, and cannot produce a contrary flag by "
           "construction. They answer different questions — “how do I make this "
           "published method usable” against “what should replace it” — and the "
           "paper can carry both, but only one can lead.")
k.numbered("**Should the bio-stat arm be kept?** If it is to be a detector it is "
           "ready and it works. If it is a contrast, it should be labelled that "
           "way. Leaving it as a third unlabelled column repeats the confusion "
           "the earlier state-aware and self-referential columns caused.")
k.numbered("**Probe masking, and what to expect from it.** It is still blocking "
           "for any biological claim. But §1 predicts it will *move* the contrary "
           "flags rather than remove them, because the concentrating sites are "
           "concentrated for a geometric reason. Both the before and after "
           "numbers should be reported.")
k.numbered("**Per-sample burden.** N15, N17, N14 and now N48 lead at every "
           "scale. Batch, detection-p, tumour purity and cell composition all "
           "remain unexcluded, and a published outlier-burden study found burden "
           "strongly confounded by exactly those. Until they are excluded, "
           "per-sample burden cannot be interpreted.")

doc.add_heading("Appendix — where every number comes from", 1)
k.para("Each script has a matching console transcript in Scripts/. Tables "
       "naming a sample together with a genomic position are individual-level "
       "and stay out of the shared repository; see CLASSIFICATION.md.")
k.table(["§", "script", "output"],
 [["1", "Task29.ContraryFlagDecomposition.Sep03.2026.R",
   "Task29_Window100_ContraryDecomposition / _BySite / _BySample.csv"],
  ["2", "Task30.StatePooledThreshold.Sep03.2026.R",
   "Task30_StatePooled_Thresholds / _QSweep / _Scores.csv"],
  ["3", "Task31.ParameterHeadroom.Sep03.2026.R",
   "Task31_EffectSizeCeiling / _PEnvelope / _PLevelResolution / _DeltVsRelMeth.csv"],
  ["4", "— (LITERATURE.md, sources cited there)", "—"],
  ["5", "Task32.ImprovedFlagging.Sep03.2026.R",
   "Task32_MeasuredNoise / _FloorsApplied / _ImprovedMethodScores.csv"],
  ["figs", "Task34.Figures.Sep03.2026.R", "Results/Fig8–Fig12*.png"]],
 widths=[0.4, 2.6, 3.3], font=8.5)

doc.save(OUT)
print("written:", OUT)
