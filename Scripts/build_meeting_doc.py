# Builds Meeting_2026-08-28_Outlier_Methylation.docx from the Results/ CSVs
# and PNGs. Regenerate after any analysis rerun so the document cannot drift
# from the numbers it cites.
#
# REQUIRES python-docx, which will not install into the system python on macOS
# (PEP 668). Use a virtualenv:
#
#     python3 -m venv .venv
#     .venv/bin/pip install python-docx
#     .venv/bin/python Scripts/build_meeting_doc.py
#
# Figure numbering in the document is by order of APPEARANCE and does not match
# the Fig<N>_*.png filenames, which are numbered by order of creation:
#     Figure 1 -> Fig7_100CpG_MethodBySite.png
#     Figure 2 -> Fig6_100CpG_WindowMap.png
#     Figure 3 -> Fig2_FlagMagnitudeByState.png
#     Figure 4 -> Fig5_N37_LocalCluster.png
#     Figure 5 -> Fig3_ThresholdGeometry.png
#     Figure 6 -> Fig4_FlagStability.png
# Fig1_FlagRateByState.png is generated but not used in the document.

import os
from docx import Document
from docx.shared import Pt, Inches, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_BREAK
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

REPO = "/Users/anamgiri/Desktop/bioinformatics-research"
OUT  = os.path.join(REPO, "Meeting_2026-08-28_Outlier_Methylation.docx")
FIG  = lambda n: os.path.join(REPO, "Results", n)

NAVY = RGBColor(0x1F, 0x38, 0x64)
GREY = RGBColor(0x60, 0x60, 0x60)
RED  = RGBColor(0xA8, 0x28, 0x28)

doc = Document()

# ---- base styles ----
st = doc.styles['Normal']
st.font.name = 'Calibri'; st.font.size = Pt(10.5)
st.paragraph_format.space_after = Pt(7)
st.paragraph_format.line_spacing = 1.12
for lvl, sz, col in ((1,17,NAVY),(2,13.5,NAVY),(3,11.5,NAVY)):
    h = doc.styles['Heading %d' % lvl]
    h.font.name='Calibri'; h.font.size=Pt(sz); h.font.bold=True; h.font.color.rgb=col
    h.paragraph_format.space_before=Pt(15 if lvl<3 else 11); h.paragraph_format.space_after=Pt(5)

for s in doc.sections:
    s.top_margin=Inches(0.85); s.bottom_margin=Inches(0.85)
    s.left_margin=Inches(0.9); s.right_margin=Inches(0.9)

def para(text=None, bold=False, italic=False, size=None, color=None,
         align=None, space_after=None, style=None):
    p = doc.add_paragraph(style=style)
    if align is not None: p.alignment = align
    if space_after is not None: p.paragraph_format.space_after = Pt(space_after)
    if text:
        for i, chunk in enumerate(text.split('**')):
            if not chunk: continue
            r = p.add_run(chunk)
            r.bold = bold or (i % 2 == 1)
            r.italic = italic
            if size: r.font.size = Pt(size)
            if color: r.font.color.rgb = color
    return p

def bullet(text, bold=False):
    p = doc.add_paragraph(style='List Bullet')
    p.paragraph_format.space_after = Pt(3)
    for i, chunk in enumerate(text.split('**')):
        if not chunk: continue
        r = p.add_run(chunk); r.bold = bold or (i % 2 == 1)
    return p

def numbered(text):
    p = doc.add_paragraph(style='List Number')
    p.paragraph_format.space_after = Pt(3)
    for i, chunk in enumerate(text.split('**')):
        if not chunk: continue
        r = p.add_run(chunk); r.bold = (i % 2 == 1)
    return p

def shade(cell, hexcolor):
    tcPr = cell._tc.get_or_add_tcPr()
    sh = OxmlElement('w:shd'); sh.set(qn('w:fill'), hexcolor); tcPr.append(sh)

def table(headers, rows, widths=None, font=9, highlight_rows=()):
    t = doc.add_table(rows=1, cols=len(headers))
    t.style = 'Table Grid'; t.alignment = WD_TABLE_ALIGNMENT.CENTER
    hdr = t.rows[0].cells
    for i, h in enumerate(headers):
        hdr[i].text = ''
        p = hdr[i].paragraphs[0]; r = p.add_run(h)
        r.bold = True; r.font.size = Pt(font); r.font.color.rgb = RGBColor(0xFF,0xFF,0xFF)
        p.paragraph_format.space_after = Pt(1)
        shade(hdr[i], '1F3864')
    for ri, row in enumerate(rows):
        cells = t.add_row().cells
        for i, v in enumerate(row):
            cells[i].text = ''
            p = cells[i].paragraphs[0]; p.paragraph_format.space_after = Pt(1)
            txt = str(v)
            for j, chunk in enumerate(txt.split('**')):
                if not chunk: continue
                r = p.add_run(chunk); r.font.size = Pt(font); r.bold = (j % 2 == 1)
            if ri in highlight_rows: shade(cells[i], 'FFF2CC')
    if widths:
        for i, w in enumerate(widths):
            for row in t.rows: row.cells[i].width = Inches(w)
    doc.add_paragraph().paragraph_format.space_after = Pt(3)
    return t

def source(text):
    p = doc.add_paragraph(); p.paragraph_format.space_after = Pt(11)
    r = p.add_run(text); r.font.size = Pt(8); r.italic = True; r.font.color.rgb = GREY
    return p

def callout(title, text, fill='F2F2F2'):
    t = doc.add_table(rows=1, cols=1); t.style = 'Table Grid'
    c = t.rows[0].cells[0]; shade(c, fill); c.text = ''
    p = c.paragraphs[0]; p.paragraph_format.space_after = Pt(3)
    r = p.add_run(title); r.bold = True; r.font.size = Pt(10); r.font.color.rgb = NAVY
    p2 = c.add_paragraph(); p2.paragraph_format.space_after = Pt(2)
    for i, chunk in enumerate(text.split('**')):
        if not chunk: continue
        rr = p2.add_run(chunk); rr.font.size = Pt(9.5); rr.bold = (i % 2 == 1)
    doc.add_paragraph().paragraph_format.space_after = Pt(4)
    return t

def figure(fname, caption, width=6.4):
    doc.add_picture(FIG(fname), width=Inches(width))
    doc.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER
    p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.space_after = Pt(12)
    for i, chunk in enumerate(caption.split('**')):
        if not chunk: continue
        r = p.add_run(chunk); r.font.size = Pt(8.5); r.font.color.rgb = GREY
        r.italic = (i % 2 == 0); r.bold = (i % 2 == 1)

def pagebreak():
    doc.add_paragraph().add_run().add_break(WD_BREAK.PAGE)

# ============================ TITLE ============================
p = para("Outlier DNA Methylation in TCGA-BRCA", size=21, bold=True,
         align=WD_ALIGN_PARAGRAPH.CENTER, space_after=2)
p.runs[0].font.color.rgb = NAVY
para("Progress report for the meeting of 28 August 2026", size=12,
     align=WD_ALIGN_PARAGRAPH.CENTER, space_after=2).runs[0].font.color.rgb = GREY
para("Anam Giri", size=10.5, align=WD_ALIGN_PARAGRAPH.CENTER, space_after=16)\
    .runs[0].font.color.rgb = GREY

callout("What this report is",
 "It answers the four questions set in the pre-meeting email, then reports three "
 "things that came out of answering them: a mechanism explaining why the outlier "
 "flags behave as they do, a bug that changes a previously reported result, and a "
 "tested recommendation for what to use instead.\n\n"
 "Every number is reproducible from a named script and CSV in the project "
 "repository. Where a claim is not yet supported, it is labelled as such rather "
 "than softened.", 'E8EEF7')

# ============================ SUMMARY ============================
doc.add_heading("Summary", 1)

para("**The four questions, answered.**")
table(["#", "Question asked", "Answer"],
[["1","Summarise chromosome 22 by methylation state",
  "Done — both tissues, all six states, three methods (§2)"],
 ["2","At high sites, is there ever a “too high” flag? At low sites, a “too low” one?",
  "Yes, routinely — up to 6.3% of measurements at tumour H sites (§3)"],
 ["3","If so: which sample, which site, and why?",
  "Three samples supply a third of them. The pattern is strongly consistent with "
  "the geometry of the threshold, though biology is not excluded (§4)"],
 ["4","Debug the per-sample summary, where N37 came top",
  "N37 ranks 41st of 53 chromosome-wide. Its top rank was an artifact of the "
  "100-site window. It does carry one candidate real event (§5)"]],
 widths=[0.32,3.05,3.15], font=9)

para("**What came out of answering them.**", space_after=4)
numbered("The flags at the methylation extremes are **statistically real but very small**. "
         "97% of flags at high-methylation sites shift methylation by less than 0.10 on a "
         "0–1 scale. This holds genome-wide, not just on chromosome 22.")
numbered("**They are also not reproducible.** Adding measurement noise the size of the "
         "array's own error changes more calls than there were flags to begin with.")
numbered("**A typo was deleting an entire methylation state** from three earlier analyses — "
         "40% of tumour sites and 65% of the tumour flags. Fixed; the corrected table "
         "changes the conclusion.")
numbered("A one-line fix — requiring a **minimum size of difference** before a flag counts — "
         "triples reproducibility. It is the strongest option of seven tested, but the "
         "exact cut-off still needs a biological justification.")

callout("The one decision I need from you",
 "Is the goal of this project to (a) evaluate whether this outlier method works, or "
 "(b) find real biological epimutations in BRCA? The work so far answers (a) well and "
 "is not yet in a position to claim (b), because probe quality-control has not been "
 "done. Committing to one changes what comes next.", 'FFF2CC')

pagebreak()

# ============================ BACKGROUND ============================
doc.add_heading("1. Background — what you need to know to read the tables", 1)

doc.add_heading("The data", 2)
para("TCGA breast cancer (BRCA): 53 normal and 53 tumour samples, measured on the "
     "Illumina 450k array at 380,355 CpG sites. Chromosome 22 holds 6,809 of them. "
     "Each measurement is a **beta value** between 0 (no methylation) and 1 (fully "
     "methylated).")
para("Each CpG site carries a label describing its usual methylation level, assigned "
     "beforehand: **L** (low), **LM** (low-medium), **M** (medium), **HM** (high-medium), "
     "**H** (high), and **R** (the remainder). An “outlier flag” marks one sample at "
     "one site as unusually high (+1), unusually low (−1), or normal (0).")

doc.add_heading("The three methods being compared", 2)
table(["Method","How it decides what counts as an outlier","Why it is here"],
[["bio","At L/LM sites, flag a sample above the 99th percentile of the 53 values. "
  "At H/HM sites, flag below the 1st percentile. No rule at M/R sites.",
  "The state-aware rule proposed in our earlier meetings"],
 ["self","Thresholds calculated from the same 53 samples being tested",
  "How the software was run in Tasks 10–15"],
 ["ext","Thresholds from a separate published panel of 2,015 independent samples",
  "How the software was designed to be used"]],
 widths=[0.6,3.4,2.5], font=9)

callout("Important: two of the three methods cannot show a result",
 "With 53 samples, the 99th percentile falls between the 52nd and 53rd ranked value. "
 "So exactly one sample can ever exceed it — the largest one. The bio and self methods "
 "therefore hand out a fixed number of flags to every site, no matter what the data "
 "says, at every significance level from 0.01 down to 0.00001.\n\n"
 "You can see this directly in the tables that follow: every self row reads "
 "1.89 / 96.23 / 1.89, which is exactly 1/53, 51/53 and 1/53.\n\n"
 "**They are kept in the tables as a negative control** — they show what “no signal” "
 "looks like. Only the ext column can carry a finding.", 'FBE9E9')

figure("Fig7_100CpG_MethodBySite.png",
 "**Figure 1.** The same 100 sites under each method, one bar per site. **self** is a flat "
 "block of exactly two flags at every single site \u2014 that is the arithmetic artifact, not a "
 "result. **bio** gives 0 or 1. Only **ext** varies with the data (0 to 8). This is why only "
 "the ext column can carry a finding.", 6.5)

doc.add_heading("Two things that cause confusion if unstated", 2)
para("**(a) What is being counted changes between sections.**")
table(["Unit","What one count means"],
[["CpG site","One probe on the array. 380,355 in total."],
 ["Measurement","One sample at one site. 53 per site."],
 ["Sample burden","All flags belonging to one sample, added up."],
 ["Event","A run of 3 or more flagged neighbouring sites within 1,000 base pairs, counted once."]],
 widths=[1.5,5.0], font=9)
para("A sample can have eleven flags and one event. That is not a contradiction — it is "
     "one biological event measured by eleven probes.", italic=True, size=9.5)

para("**(b) The external panel does not cover every site, and the gaps are uneven.**")
para("Sites missing from the 2,015-sample panel return no answer, so the ext method is "
     "scored on fewer measurements than the other two. Coverage differs by state, which "
     "matters when comparing states to each other:")
table(["Methylation state","Sites missing from panel — Normal","Tumour"],
[["H","1.05%","0.46%"],["L","2.94%","2.34%"],["HM","3.30%","3.38%"],
 ["R","8.91%","8.27%"],["LM","**14.39%**","**11.47%**"],["M","**16.29%**","**24.82%**"]],
 widths=[2.2,2.4,1.9], font=9)
source("Source: Task26_DenominatorAudit_*.csv  |  Script: Task26.ReviewResponse.Aug28.2026.R")
para("All percentages in this report are calculated over sites that were actually "
     "measurable, so each is internally correct. But comparing raw counts between "
     "methods at LM or M sites compares different denominators.", size=9.5, italic=True)

pagebreak()

# ============================ Q1 ============================
doc.add_heading("2. Question 1 — Chromosome 22 summarised by state", 1)
para("The email's example — “H, 11 CG sites” — identifies the window exactly: the first "
     "100 chromosome-22 sites by position, in which 11 are H-state. Reproduced in full.")
para("**Normal tissue, 100 sites, chr22:11,915,060–17,085,407.** Read the ext rows.")
table(["State","Sites","Method","−1","0","+1","% −1","% 0","% +1"],
[["L","8","bio","0","416","8","0.00","98.11","1.89"],
 ["","","self","8","408","8","1.89","96.23","1.89"],
 ["","","**ext**","**8**","359","4","**2.16**","96.77","1.08"],
 ["LM","14","bio","0","728","14","0.00","98.11","1.89"],
 ["","","self","14","714","14","1.89","96.23","1.89"],
 ["","","**ext**","**21**","718","3","**2.83**","96.77","0.40"],
 ["M","24","bio","0","1272","0","0.00","100.00","0.00"],
 ["","","self","24","1224","24","1.89","96.23","1.89"],
 ["","","**ext**","12","1241","19","0.94","97.56","1.49"],
 ["HM","38","bio","38","1976","0","1.89","98.11","0.00"],
 ["","","self","38","1938","38","1.89","96.23","1.89"],
 ["","","**ext**","16","1952","**46**","0.79","96.92","**2.28**"],
 ["H","11","bio","11","572","0","1.89","98.11","0.00"],
 ["","","self","11","561","11","1.89","96.23","1.89"],
 ["","","**ext**","2","571","**10**","0.34","97.94","**1.72**"],
 ["R","5","bio","0","265","0","0.00","100.00","0.00"],
 ["","","self","5","255","5","1.89","96.23","1.89"],
 ["","","**ext**","5","256","4","1.89","96.60","1.51"]],
 widths=[0.6,0.5,0.65,0.55,0.6,0.55,0.7,0.7,0.7], font=8.5,
 highlight_rows=(2,5,8,11,14,17))
source("Source: Task20_Chr22_index100_StateSummary_Normal.csv  |  Script: Task20.Chr22.TightWindow.Concordance.Aug28.2026.R")
para("Notice every self row: 1.89 / 96.23 / 1.89, identical in all six states. That is the "
     "arithmetic artifact described on page 2, visible in the data.")

figure("Fig6_100CpG_WindowMap.png",
 "**Figure 2.** All 100 sites at once. Top: the methylation landscape \u2014 each vertical bar is "
 "the middle 50% of the 53 samples, the dot is the median, coloured by state. Middle: the state "
 "of each site in normal and tumour tissue, showing how differently they are labelled. Bottom: "
 "how many samples the external reference flags at each site. Tumour flags far more sites "
 "(97 of 100) and far more samples per site than normal (65 of 100).", 6.5)

para("**Tumour tissue, same 100 sites.** The state labels are assigned separately per "
     "tissue, so these sites sit in different states here: R grows from 5 sites to 52.")
table(["State","Sites","−1","0","+1","% −1","% 0","% +1"],
[["L","7","4","307","7","1.26","96.54","2.20"],
 ["LM","8","29","362","33","6.84","85.38","7.78"],
 ["M","13","178","473","38","25.84","68.65","5.52"],
 ["HM","19","204","760","43","20.26","75.47","4.27"],
 ["H","**1**","14","39","0","26.42","73.59","0.00"],
 ["R","52","551","2010","195","19.99","72.93","7.08"]],
 widths=[0.75,0.7,0.7,0.75,0.7,0.85,0.85,0.85], font=9)
callout("Do not read the tumour H row",
 "It contains a single CpG site. “26.42%” means 14 flagged measurements at one probe. "
 "A percentage over one site estimates nothing. The L and LM rows, at 7 and 8 sites, are "
 "also thin. This is a property of the window, not of tumour biology — the tumour "
 "annotation moves most of these sites into R. The genome-wide table in section 7, where "
 "H has 30,080 sites, is the one to quote.", 'FBE9E9')

pagebreak()

# ============================ Q2 ============================
doc.add_heading("3. Question 2 — Are there flags pointing where there is no room?", 1)
para("A “too high” flag at a site that is already highly methylated, or a “too low” flag "
     "at one that is already low, points in the direction with the least room to move. "
     "You asked whether any exist.")
para("**Yes — in every state, in both tissues, at every scale tested.**")
table(["Tissue","State","Flag that contradicts the state","Count","of measurements","Rate"],
[["Normal","L","too low (−1)","39","2,120","1.84%"],
 ["Normal","LM","too low (−1)","5","530","0.94%"],
 ["Normal","HM","too high (+1)","9","1,113","0.81%"],
 ["Normal","H","too high (+1)","16","848","1.89%"],
 ["Tumour","L","too low (−1)","52","2,014","2.58%"],
 ["Tumour","LM","too low (−1)","33","583","**5.66%**"],
 ["Tumour","HM","too high (+1)","47","1,060","**4.43%**"],
 ["Tumour","H","too high (+1)","40","636","**6.29%**"]],
 widths=[0.85,0.6,1.9,0.7,1.35,0.85], font=9)
source("External reference, 100-site tight window (see section 6). Source: Task20_Chr22_tight100_ContraryCounts_*.csv")
callout("A necessary qualification",
 "“Contradicts the state” means the flag points opposite to a six-level summary label — "
 "not that it is biologically impossible. Real variation within a state exists, and a "
 "genuinely under-methylated sample at a low site is perfectly possible.\n\n"
 "What makes these flags questionable is not their direction but their **size**, which is "
 "the subject of the next section.\n\n"
 "Note also that the bio rule returns zero in every row here — by construction it can only "
 "flag in the state-consistent direction, so its zero is not evidence of anything.")

# ============================ Q3 ============================
doc.add_heading("4. Question 3 — Which sample, which site, and why", 1)
doc.add_heading("Which samples", 2)
para("Three samples supply about a third of all state-contradicting flags on chromosome 22, "
     "and they are the same three in both tissues.")
table(["Sample","Contradicting flags","All its flags","% contradicting"],
[["N17","410","1,250","32.8%"],["N15","381","1,297","29.4%"],
 ["N14","369","1,039","35.5%"],["N48","356","496","71.8%"],
 ["N31","301","405","74.3%"],["N27","264","278","**95.0%**"]],
 widths=[1.0,1.75,1.5,1.6], font=9)
source("Source: Task21_Chr22_ContraryDrivers_Normal.csv  |  Script: Task21.Chr22.SampleBurden.Aug28.2026.R")
para("N15, N17 and N14 also lead the flag count at every scale in both tissues, and their "
     "chromosome-22 burden tracks their genome-wide burden closely (rank correlation 0.76 "
     "normal, 0.84 tumour). **This is a property of those whole samples, not of chromosome 22** "
     "— plausibly global methylation loss, tumour purity, or a technical batch effect. "
     "None of these has been tested yet.")

doc.add_heading("Which sites", 2)
para("Of the 69 contradicting flags in the tight window: **none** fall in a CpG island, "
     "**none** have a known genetic variant within 10 base pairs, only 2 are also outliers "
     "by the standard Tukey test, and 48 are isolated with no flagged neighbour nearby. "
     "The usual technical explanations do not apply.")

doc.add_heading("Why — the flagged differences are very small", 2)
para("For every flagged measurement, how far is it from the middle of the cohort? "
     "Measured genome-wide across all 380,355 sites:")
table(["State","Tissue","Sites","Flags","Median difference","% smaller than 0.10"],
[["H","Normal","61,845","78,259","**0.031**","**96.97%**"],
 ["L","Normal","104,540","174,567","**0.042**","81.53%"],
 ["HM","Normal","97,094","114,575","0.119","39.82%"],
 ["M","Normal","36,785","25,087","0.153","16.98%"],
 ["R","Normal","23,455","38,131","**0.269**","3.68%"],
 ["H","Tumour","30,080","68,114","0.028","**96.16%**"],
 ["L","Tumour","76,106","163,788","0.031","85.21%"],
 ["R","Tumour","153,479","1,608,650","**0.213**","19.92%"]],
 widths=[0.6,0.8,0.95,1.15,1.4,1.5], font=9)
source("Source: Task26_GenomeWideMagnitude.csv  |  Script: Task26.ReviewResponse.Aug28.2026.R")
para("**Nearly every flag at a high-methylation site corresponds to a change too small to "
     "matter biologically.** At R-state sites, the same method produces flags eight times "
     "larger. The chromosome-22 figures reproduce genome-wide to within one percentage point.")
figure("Fig2_FlagMagnitudeByState.png",
 "**Figure 3.** How big is a flagged difference? Bars show the median change in "
 "methylation for flagged measurements. The dashed line at 0.10 is a common threshold for "
 "“biologically meaningful”; the dotted line at 0.15 is the cut-off used by the published "
 "epimutacions software. Flags at L and H sites fall far below both.")

pagebreak()

# ============================ Q4 ============================
doc.add_heading("5. Question 4 — Debugging the per-sample summary (N37)", 1)
para("N37 topped the earlier summary with 11 of 71 flags, about eight times the expected "
     "share. Widening the comparison settles it:")
table(["Scale","Sites","Method","Flags","Expected","Rank of 53"],
[["100-site window","100","bio","11","1.3","**1st**"],
 ["100-site window","100","ext","1","2.8","30th"],
 ["tight 100-site window","100","bio","1","1.7","17th"],
 ["whole chromosome 22","6,809","bio","24","108.3","**41st**"],
 ["whole chromosome 22","6,809","ext","114","182.0","25th"]],
 widths=[1.85,0.75,0.75,0.7,0.9,1.05], font=9)
source("Source: Task21_Chr22_SampleBurden_Normal.csv")
para("**Three separate reasons N37 came top, none of them “N37 is an unusual sample”:**")
bullet("**The window was not representative.** Across the whole chromosome N37 carries less "
       "than a quarter of its expected flags and ranks 41st of 53.")
bullet("**The external reference does not see it** — one flag, ranked 30th.")
bullet("**Ten of the flags are one event.** They collapse into three runs of neighbouring "
       "sites; the largest holds seven consecutive probes.")

doc.add_heading("But the underlying event looks real", 2)
para("Seven consecutive sites across 1,740 base pairs, where N37 sits far above everyone else:")
table(["Site","Position (chr22)","N37","Cohort median","Next highest","Margin"],
[["cg01836687","16,601,097","0.553","0.052","0.141","**0.413**"],
 ["cg26822097","16,601,691","0.402","0.039","0.341","0.061"],
 ["cg12431879","16,601,897","0.442","0.089","0.289","0.152"],
 ["cg00332021","16,602,522","0.435","0.048","0.213","0.222"],
 ["cg00816224","16,602,592","0.438","0.046","0.243","0.195"],
 ["cg15554678","16,602,673","0.356","0.015","0.111","0.245"],
 ["cg23572163","16,602,837","0.276","0.065","0.211","0.064"]],
 widths=[1.15,1.3,0.75,1.15,1.05,0.85], font=9)
source("Source: Task21_Chr22_N37_Cluster_Normal.csv (restricted file — contains individual-level values)")
figure("Fig5_N37_LocalCluster.png",
 "**Figure 4.** N37 (red) against the other 52 samples (grey) across this region. The "
 "shaded band marks a stretch with no probes on the array; lines are broken there rather "
 "than drawn across it. Tick marks show actual probe positions.")
callout("This must not be called a real finding yet",
 "The region sits at chr22:16.6 Mb, next to the centromere on a repetitive, poorly-mapped "
 "part of the chromosome. A run of seven elevated probes there is also exactly what a "
 "probe cross-hybridisation artifact looks like. Standard probe-quality filtering has not "
 "been applied to this project.\n\n"
 "**Correct wording: a candidate focal event, pending probe quality control.** What can be "
 "said without qualification is narrower — the per-sample summary counted one continuous "
 "run seven separate times.", 'FBE9E9')

pagebreak()

# ============================ WINDOW + AGREEMENT ============================
doc.add_heading("6. Two supporting checks", 1)
doc.add_heading("The original 100 sites were not close together", 2)
para("The email asked for sites close to each other. The 100 used previously are "
     "consecutive by position in the file, not by distance:")
table(["Window","Total span","Median gap","Largest gap","Neighbours within 1 kb"],
[["Original 100","**5,170,347 bp**","512 bp","**3,410,179 bp**","61 of 99"],
 ["New tight 100","**73,399 bp**","187 bp","13,823 bp","80 of 99"]],
 widths=[1.4,1.35,1.05,1.3,1.5], font=9)
source("Source: Task20_Chr22_WindowComparison.csv")
para("The new window is the 100 consecutive sites with the smallest span on chromosome 22 "
     "(chr22:50,459,312–50,532,711). Both are analysed side by side throughout, so the "
     "earlier result stays reproducible. The choice matters: N37 falls from 1st to 17th "
     "between them.")

doc.add_heading("The state-aware rule and the external reference disagree", 2)
para("Simple percent agreement is misleading here, because more than 95% of all "
     "measurements are “not flagged” by both methods — so any two methods agree about 96% "
     "of the time just by both saying nothing. Two better measures:")
bullet("**Agreement among flagged** — of measurements flagged by at least one method, how "
       "often both gave the same answer.")
bullet("**Cohen's kappa** — a standard agreement score corrected for agreement expected by "
       "chance. 0 means chance-level; 1 means perfect.")
table(["Tissue","Methods compared","Flagged by either","Both agreed","Agreement","Kappa"],
[["Normal","bio vs ext","172","21","**12.2%**","**0.21**"],
 ["Normal","self vs ext","235","63","26.8%","0.41"],
 ["Tumour","bio vs ext","415","28","**6.8%**","**0.11**"],
 ["Tumour","self vs ext","455","99","21.8%","0.33"]],
 widths=[0.8,1.5,1.35,1.05,1.0,0.85], font=9)
source("Source: Task20_Chr22_tight100_Concordance_*.csv")
para("The two approaches **define “outlier” differently** rather than offering two views of "
     "the same underlying set. One caution: with this much imbalance, kappa is sensitive to "
     "the overall flag rates and should be read alongside the agreement column. Both point "
     "the same way here.")

pagebreak()

# ============================ MECHANISM ============================
doc.add_heading("7. Why the flags behave this way", 1)
para("The thresholds the software used can be recovered from its own output, which also "
     "confirms the pipeline is correctly aligned: across 13,618 site-by-tissue combinations "
     "there was **not one case** where a flagged sample's value sat inside the unflagged range.")

doc.add_heading("The outlier zones are extremely lopsided", 2)
para("For each state, how much of the 0–1 methylation scale counts as “too high”, and how "
     "much as “too low”:")
table(["State","“Too high” zone width","“Too low” zone width","Gap between the two samples either side of the line"],
[["L","0.900","**0.021**","0.024"],["LM","0.774","0.049","0.057"],
 ["M","0.318","0.360","0.028"],["HM","**0.090**","0.677","0.013"],
 ["H","**0.037**","0.886","**0.003**"],["R","0.238","0.269","0.052"]],
 widths=[0.65,1.6,1.6,2.5], font=9)
source("Normal tissue, chromosome 22 medians. Source: Task22_Chr22_ThresholdGeometry_Normal.csv")
para("At an H site the “too high” zone is 0.037 wide and the “too low” zone is 0.886 — "
     "24 times larger. Yet **87.8% of the flags that actually occur are “too high”.** "
     "Adjusting for the room available:")
callout("The calculation, in full",
 "There are 1,283 flags at H sites in normal tissue.\n"
 "    “Too high”:  1,283 × 87.8% = 1,126 flags in a zone 0.037 wide  →  30,363 per unit\n"
 "    “Too low”:   1,283 × 12.2% =   157 flags in a zone 0.886 wide  →      177 per unit\n\n"
 "**Flags are 171 times denser in the direction the site's own state already restricts.** "
 "At L sites the same calculation gives 44 times. A method with no directional preference "
 "would give roughly 1.")
figure("Fig3_ThresholdGeometry.png",
 "**Figure 5.** Where the outlier zones sit. Each row is one state. The coloured bands are "
 "**threshold zones, not data** — blue is the range counted as too low, red as too high, "
 "grey as normal. The number beside each red band is its width. Compare the red band at H "
 "(0.037) with the one at L (0.900).")

doc.add_heading("The dividing line falls inside the array's own measurement error", 2)
para("At H sites, the two samples either side of the threshold differ by **0.003**. The "
     "450k array's repeat-measurement error is roughly 0.01–0.03 — several times larger. "
     "Adding artificial noise and re-running the flagging:")
table(["Noise added (sd)","Flags lost","Flags gained","% of original flags lost","Changed calls ÷ original flags"],
[["0.005","1,487","6,833","15.4%","0.86"],
 ["**0.010**","**2,062**","**13,467**","**21.4%**","**1.61**"],
 ["0.020","2,693","23,187","27.9%","2.68"],
 ["0.050","3,485","43,726","36.1%","4.91"]],
 widths=[1.25,1.0,1.15,1.75,1.85], font=9)
source("Normal tissue, chromosome 22, 9,648 original flags, 20 repeats per level. Source: Task22_Chr22_FlagStability_Normal.csv")
para("At noise of 0.01 — within the array's own error — **more calls change than there were "
     "flags to begin with.** Instability follows the compressed states exactly: H 34%, "
     "L 25%, against R 8%.")
figure("Fig4_FlagStability.png",
 "**Figure 6.** How many flags change when realistic measurement noise is added. Each line "
 "is one methylation state. Higher means less reproducible.")

doc.add_heading("This is not a fault in one software package", 2)
para("Tukey's standard outlier rule — a completely different statistic, and the one already "
     "used for the outliers.coef3 column in our source files — shows a state bias too, in a "
     "different direction, and the two rules barely overlap (Jaccard 0.05–0.29).")
para("The reason is that the spread it measures is itself state-dependent: the interquartile "
     "range is 0.008 at L sites but 0.157 at R sites. **Any rule that asks “who is in the "
     "extreme tail?” without asking “by how much?” inherits this problem**, because "
     "methylation values are squeezed against 0 and 1 and vary far less at the ends than in "
     "the middle.")
callout("Two alternative explanations were tested and ruled out",
 "**Was the external panel simply a poor fit for this cohort?** If so, many samples would be "
 "flagged at the same site. They are not — at most 11 of 53, no site reaches half the "
 "cohort, and 0% of normal flags sit in such sites.\n\n"
 "**Was there a data-alignment error?** No — zero inconsistencies across 13,618 checks.")

pagebreak()

# ============================ BUG ============================
doc.add_heading("8. A bug that changes a previously reported result", 1)
para("Three earlier analyses (Tasks 10, 13 and 14) contained this line:")
p = doc.add_paragraph(); p.paragraph_format.left_indent = Inches(0.3)
r = p.add_run('state_order <- c("L", "LM", "M", "HM", "H", "Rc")')
r.font.name='Consolas'; r.font.size=Pt(9.5)
r2 = p.add_run('      ← the data says "R", not "Rc"'); r2.font.size=Pt(9.5); r2.italic=True; r2.font.color.rgb=RED
para("The mismatch meant every R-state site was **silently dropped** from those tables — no "
     "error, no warning, no change in row count to notice.")
table(["Reference","Tissue","R sites dropped","% of all sites","% of all flags lost"],
[["self","Normal","23,455","6.2%","6.2%"],
 ["self","Tumour","153,479","40.4%","40.4%"],
 ["ext","Normal","23,455","6.2%","7.4%"],
 ["**ext**","**Tumour**","**153,479**","**40.4%**","**65.0%**"]],
 widths=[1.1,1.0,1.35,1.3,1.55], font=9, highlight_rows=(3,))
source("Source: Task23_RStateImpact.csv  |  Script: Task23.CorrectedStateTables.Aug28.2026.R")
para("**The last row is the important one.** Task 14's tumour table — the strongest result "
     "we had — was calculated from 35% of the flags. Corrected, genome-wide:")
table(["State","Sites","% too low","% normal","% too high","% flagged"],
[["L","76,106","2.08","95.84","2.08","4.16"],
 ["LM","40,824","3.41","89.66","6.93","10.34"],
 ["M","5,436","8.61","82.63","8.76","17.37"],
 ["HM","74,430","4.55","89.56","5.89","10.44"],
 ["H","30,080","0.50","95.71","3.79","4.29"],
 ["**R**","**153,479**","**10.14**","**78.44**","**11.42**","**21.56%**"]],
 widths=[0.7,1.1,1.0,1.0,1.05,1.05], font=9, highlight_rows=(5,))
source("Tumour, external reference, all 380,355 sites. Source: Task23_StateFlagTable_extref_Tumor.csv")
para("The five original rows reproduce exactly, which confirms the correction is sound. "
     "**The row that had been deleted turns out to have the highest flag rate of any state, "
     "and (from section 4) the largest effect sizes.** The category being thrown away was "
     "the one carrying the signal.")
para("All three scripts are now fixed and carry a check that stops with an error if any "
     "state in the data is missing from the list, so this cannot recur silently.")

pagebreak()

# ============================ RECOMMENDATION ============================
doc.add_heading("9. What to use instead", 1)
para("Seven different outlier definitions were run on the same data. Making that "
     "comparison fair matters more than it sounds: a method that flags fewer measurements "
     "will look more reproducible for free. So every method below was tuned until it "
     "flagged **exactly the same percentage** of measurements (1.03%) before being scored.")
para("Judged on: **reproducibility** (does the same flag survive realistic measurement "
     "noise? higher is better), whether flags at H sites are **as large** as those at R "
     "sites (1.0 would be even-handed), and the share of flags **too small to matter** "
     "(under 0.10).")
table(["Method","Reproducibility","H vs R effect size","% of flags too small"],
[["**External reference + minimum-size rule**","**0.79**","**0.49**","**0%**"],
 ["Tukey interquartile rule","0.47","0.17","54%"],
 ["Median \u00b1 robust SD","0.42","0.16","53%"],
 ["Beta-distribution fit","0.22","0.16","40%"],
 ["Median \u00b1 robust SD on log-odds scale","0.19","0.16","47%"]],
 widths=[2.9,1.2,1.35,1.4], font=9, highlight_rows=(0,))
source("Normal tissue; every method tuned to flag 1.031% of measurements. "
       "Source: Task26_RateMatchedStability_Normal.csv")
para("Tumour behaves the same way: 0.91 for the minimum-size rule, against 0.68, 0.65, "
     "0.54 and 0.43 for the others.")
para("A wider sweep of seven methods, run instead at the external reference's own flag rate "
     "(Task24_MethodBenchmark_Scores.csv), ranks them in the same order and adds two more: "
     "a leave-one-out percentile (0.40) and the current external reference unmodified "
     "(**0.33**, with 61% of its flags too small to matter).")

doc.add_heading("The recommendation, and its limits", 2)
para("**Adding a minimum-size requirement to the existing method is the strongest option "
     "tested.** It means: before a flag counts, the sample must differ from the cohort "
     "middle by at least some fixed amount. This more than doubles reproducibility "
     "(0.33 to 0.79 at matched rates) and removes the too-small flags by definition. It is "
     "one extra line of code, and the published epimutacions software already does it.")
para("**Two objections were tested rather than argued away:**")
para("**Is a cut-off of 0.10 justified?** No — it was chosen after seeing the data. Sweeping "
     "it across a range:", space_after=4)
table(["Minimum size","Flags kept","Reproducibility","H vs R effect size"],
[["0 (current)","100%","0.33","0.12"],["0.05","63%","0.75","0.24"],
 ["0.10","39%","0.79","0.49"],["0.15","20%","0.82","0.53"],
 ["0.20","12%","0.87","not defined"]],
 widths=[1.35,1.15,1.4,1.6], font=9)
source("Source: Task26_FloorSensitivity_Normal.csv")
para("Reproducibility rises steadily with no peak, so **this table cannot choose the value** "
     "— pushed far enough it flags almost nothing and scores near-perfectly. At 0.20 no H-site "
     "flag survives at all. **The cut-off has to be justified biologically**, as the smallest "
     "methylation difference worth calling real on this platform. epimutacions uses 0.15; the "
     "difference between 0.10 and 0.15 here is small.")
para("**Is the advantage just because it flags less?** No — this is why the table above is "
     "rate-matched. The minimum-size rule does naturally flag less (1.03% against 2.82% for "
     "the current method), so every competitor was re-tuned down to its rate before scoring. "
     "It still leads by a wide margin. The advantage belongs to the method, not to its "
     "sensitivity.")
callout("One result that reverses my earlier recommendation",
 "In an earlier draft I recommended converting to the log-odds (M-value) scale — the "
 "textbook correction for this kind of problem. **Tested, it is the worst of the seven "
 "methods** (0.31 against 0.56 for the same rule on the original scale).\n\n"
 "The reason: the log-odds transform stabilises the spread but magnifies noise at the "
 "extremes — a shift of 0.01 near 0.97 becomes a large move — and the extremes are exactly "
 "where the problem lives. Recorded as a negative result rather than dropped.", 'FFF2CC')

pagebreak()

# ============================ LIMITATIONS ============================
doc.add_heading("10. What is not yet established", 1)
para("Listed plainly, because several affect how much weight the findings above can carry.")
table(["Open question","Status"],
[["**Are the 53 normal and tumour samples genuinely from the same patients?**",
  "Cannot confirm. A methylation-based identity check matches the expected pairing in "
  "only 10 of 53 cases (chance would give ~1), so the ordering probably does encode it, "
  "but this is far from proof. Nothing in this report relies on pairing."],
 ["**Is the external 2,015-sample panel technically comparable to our data?**",
  "Not established. Array version, preprocessing, batch structure and the definition of "
  "“tumour-adjacent” are all unverified. This is the largest open risk to the external "
  "results. Gross mismatch is ruled out, but that is not the same as comparability."],
 ["**Are the flagged sites real, or probe artifacts?**",
  "Unknown. Standard filtering for cross-reacting and variant-overlapping probes has not "
  "been applied. This blocks any biological claim, including N37's cluster."],
 ["**Why do N15, N17 and N14 carry so many flags?**",
  "Untested. Candidates in order: technical batch and array quality; cell composition and "
  "tumour purity; genuine global methylation loss. Until the first two are excluded, "
  "per-sample burden cannot be interpreted."],
 ["**Does the threshold-geometry finding hold genome-wide?**",
  "The effect-size result does (shown in section 4). The threshold-geometry table and the "
  "noise experiment are still chromosome-22 only."],
 ["**Is the noise level used in the stability test realistic?**",
  "It comes from published figures, not from this cohort. Technical replicates, or pairs "
  "of adjacent co-methylated probes, would replace the assumption with a measurement."]],
 widths=[2.55,4.3], font=9)

doc.add_heading("11. Proposed next steps", 1)
numbered("**Apply standard probe filtering** (Chen 2013; Zhou 2017). This blocks every "
         "biological claim, so it goes first.")
numbered("**Test what drives the high-burden samples** — batch, detection quality, cell "
         "composition, tumour purity.")
numbered("**Measure this cohort's actual noise level** instead of assuming it, and re-run "
         "the stability analysis.")
numbered("**Extend the threshold-geometry and stability analyses genome-wide.**")
numbered("**Decide the minimum-size cut-off on biological grounds**, then re-run the full "
         "pipeline once with it applied.")

doc.add_heading("Questions for you", 2)
bullet("**Is the goal to evaluate the method, or to find real epimutations?** The work "
       "supports the first; the second needs steps 1 and 2 above first.")
bullet("**Should the state-aware (bio) rule be retired?** As written it cannot produce a "
       "state-contradicting flag, so comparing it against the external reference measures a "
       "difference in definitions rather than a disagreement about the data. It is useful as "
       "a negative control but not as a detector.")
bullet("**What minimum methylation difference should count as real** on the 450k platform? "
       "This is the one number the analysis cannot supply for itself.")

# ============================ APPENDIX ============================
doc.add_heading("Appendix A — where every number comes from", 1)
para("Every table in this report is produced by one named script and saved as one named "
     "CSV. Each script has a matching .Rout transcript recording its inputs, parameters and "
     "console output.")
table(["Script","What it produces","Used in"],
[["Task20.Chr22.TightWindow.Concordance","State summaries, agreement measures, tight window","§2, §3, §6"],
 ["Task21.Chr22.SampleBurden","Per-sample burden at three scales, effect sizes, N37 cluster","§4, §5"],
 ["Task22.Chr22.FlagDefects","Threshold zones, noise stability, Tukey comparison","§7"],
 ["Task23.CorrectedStateTables","Corrected state tables with R restored","§8"],
 ["Task24.OutlierMethodBenchmark","Seven-method comparison","§9"],
 ["Task25.Figures","Figures 3–6","§4, §5, §7"],
 ["Task26.ReviewResponse","Coverage audit, cut-off sweep, rate-matched comparison, genome-wide effect sizes","§1, §4, §9"],
 ["Task27.Full100CpGSheet","The complete 100-site sheet and Figures 1\u20132","§1, §2, App. B"]],
 widths=[2.3,3.35,1.2], font=8.5)
para("**Data protection.** Tables that name both an individual sample and a genomic position "
     "are individual-level data and are excluded from the shared repository. Summary tables "
     "reporting counts per sample, with no genomic coordinates, are included. The full "
     "classification is in CLASSIFICATION.md.", size=9)
para("Companion documents: **results.md** (full technical record), **REVIEW.md** (response to "
     "your written review, plus code findings), **plan.md** (ordered next steps).", size=9)


# ============ APPENDIX B: the full 100-site sheet ============
pagebreak()
doc.add_heading("Appendix B \u2014 the complete 100-CpG sheet", 1)
para("Every site in the window, in genomic order. Flag columns are written as "
     "**too-low / too-high**, counting how many of the 53 samples were flagged in each "
     "direction. N = normal tissue, T = tumour.")
para("Reading it: **self** reads 1/1 on every single row \u2014 that is the arithmetic "
     "artifact. **bio** reads 1/0 or 0/1 or 0/0. Only **ext** varies with the data. The gap "
     "column shows how far the site is from the previous one; note row 7, where the gap is "
     "3.41 million base pairs.", size=9.5)

import csv
rows = list(csv.DictReader(open(os.path.join(REPO,"Results","Task27_Chr22_Full100CpG_SiteTable.csv"))))
def g(r,k):
    v = r.get(k,"")
    return v if v not in ("NA","") else "\u2014"
def fmtgap(v):
    if v in ("NA",""): return "\u2014"
    return format(int(float(v)), ",")
body = []
for r in rows:
    body.append([
        g(r,"site.index"), g(r,"cgID"), format(int(r["pos"]), ","), fmtgap(r["gap.to.prev.bp"]),
        g(r,"state.normal"), g(r,"state.tumor"), g(r,"median.beta.N"),
        "%s/%s" % (g(r,"bio.neg1.N"),  g(r,"bio.pos1.N")),
        "%s/%s" % (g(r,"self.neg1.N"), g(r,"self.pos1.N")),
        "%s/%s" % (g(r,"ext.neg1.N"),  g(r,"ext.pos1.N")),
        "%s/%s" % (g(r,"ext.neg1.T"),  g(r,"ext.pos1.T")),
    ])
table(["#","CpG probe","chr22 position","Gap (bp)","State N","State T",
       "Median \u03b2 (N)","bio N","self N","ext N","ext T"],
      body,
      widths=[0.28,0.86,0.85,0.72,0.5,0.5,0.6,0.48,0.48,0.48,0.48], font=7)
source("Source: Task27_Chr22_Full100CpG_SiteTable.csv  |  Script: Task27.Full100CpGSheet.Aug28.2026.R  |  "
       "The CSV additionally carries the 25th and 75th percentiles of beta for both tissues and the "
       "number of evaluable samples per site. Cohort minimum and maximum are deliberately excluded: "
       "with 53 samples those are single individuals\u2019 values.")

doc.save(OUT)
print("saved:", OUT)
print("paragraphs:", len(doc.paragraphs), " tables:", len(doc.tables))
