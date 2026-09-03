"""Shared python-docx helpers for the meeting reports.

Scripts/build_meeting_doc.py (the 2026-08-28 report) carries its own copy of
these, written inline. It is left alone deliberately - it is committed, it
runs, and its output has already been shown to the supervisor, so there is no
reason to risk it for a refactor. New builders import from here instead.

Every helper takes the Document it writes into, so one process can build more
than one report. Text arguments support **bold** spans, matching the
convention the original builder used.

    python3 -m venv .venv && .venv/bin/pip install python-docx
    .venv/bin/python Scripts/build_meeting_doc_sep03.py
"""

import os
from docx.shared import Pt, Inches, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_BREAK
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

NAVY = RGBColor(0x1F, 0x38, 0x64)
GREY = RGBColor(0x60, 0x60, 0x60)
RED  = RGBColor(0xA8, 0x28, 0x28)
CENTER = WD_ALIGN_PARAGRAPH.CENTER


def base_styles(doc):
    st = doc.styles['Normal']
    st.font.name = 'Calibri'; st.font.size = Pt(10.5)
    st.paragraph_format.space_after = Pt(7)
    st.paragraph_format.line_spacing = 1.12
    for lvl, sz, col in ((1, 17, NAVY), (2, 13.5, NAVY), (3, 11.5, NAVY)):
        h = doc.styles['Heading %d' % lvl]
        h.font.name = 'Calibri'; h.font.size = Pt(sz)
        h.font.bold = True; h.font.color.rgb = col
        h.paragraph_format.space_before = Pt(15 if lvl < 3 else 11)
        h.paragraph_format.space_after = Pt(5)
    for s in doc.sections:
        s.top_margin = Inches(0.85); s.bottom_margin = Inches(0.85)
        s.left_margin = Inches(0.9); s.right_margin = Inches(0.9)


def _runs(p, text, bold=False, italic=False, size=None, color=None):
    for i, chunk in enumerate(text.split('**')):
        if not chunk:
            continue
        r = p.add_run(chunk)
        r.bold = bold or (i % 2 == 1)
        r.italic = italic
        if size:
            r.font.size = Pt(size)
        if color:
            r.font.color.rgb = color
    return p


class Kit:
    def __init__(self, doc, fig_dir):
        self.doc = doc
        self.fig_dir = fig_dir

    def para(self, text=None, bold=False, italic=False, size=None, color=None,
             align=None, space_after=None, style=None):
        p = self.doc.add_paragraph(style=style)
        if align is not None:
            p.alignment = align
        if space_after is not None:
            p.paragraph_format.space_after = Pt(space_after)
        if text:
            _runs(p, text, bold, italic, size, color)
        return p

    def heading(self, text, level=1):
        return self.doc.add_heading(text, level)

    def bullet(self, text, bold=False):
        p = self.doc.add_paragraph(style='List Bullet')
        p.paragraph_format.space_after = Pt(3)
        return _runs(p, text, bold)

    def numbered(self, text):
        p = self.doc.add_paragraph(style='List Number')
        p.paragraph_format.space_after = Pt(3)
        return _runs(p, text)

    @staticmethod
    def shade(cell, hexcolor):
        tcPr = cell._tc.get_or_add_tcPr()
        sh = OxmlElement('w:shd'); sh.set(qn('w:fill'), hexcolor); tcPr.append(sh)

    def table(self, headers, rows, widths=None, font=9, highlight_rows=()):
        t = self.doc.add_table(rows=1, cols=len(headers))
        t.style = 'Table Grid'; t.alignment = WD_TABLE_ALIGNMENT.CENTER
        hdr = t.rows[0].cells
        for i, h in enumerate(headers):
            hdr[i].text = ''
            p = hdr[i].paragraphs[0]
            r = p.add_run(h)
            r.bold = True; r.font.size = Pt(font)
            r.font.color.rgb = RGBColor(0xFF, 0xFF, 0xFF)
            p.paragraph_format.space_after = Pt(1)
            self.shade(hdr[i], '1F3864')
        for ri, row in enumerate(rows):
            cells = t.add_row().cells
            for i, v in enumerate(row):
                cells[i].text = ''
                p = cells[i].paragraphs[0]
                p.paragraph_format.space_after = Pt(1)
                for j, chunk in enumerate(str(v).split('**')):
                    if not chunk:
                        continue
                    r = p.add_run(chunk)
                    r.font.size = Pt(font); r.bold = (j % 2 == 1)
                if ri in highlight_rows:
                    self.shade(cells[i], 'FFF2CC')
        if widths:
            for i, w in enumerate(widths):
                for row in t.rows:
                    row.cells[i].width = Inches(w)
        self.doc.add_paragraph().paragraph_format.space_after = Pt(3)
        return t

    def source(self, text):
        p = self.doc.add_paragraph()
        p.paragraph_format.space_after = Pt(11)
        r = p.add_run(text)
        r.font.size = Pt(8); r.italic = True; r.font.color.rgb = GREY
        return p

    def callout(self, title, text, fill='F2F2F2'):
        t = self.doc.add_table(rows=1, cols=1); t.style = 'Table Grid'
        c = t.rows[0].cells[0]; self.shade(c, fill); c.text = ''
        p = c.paragraphs[0]; p.paragraph_format.space_after = Pt(3)
        r = p.add_run(title)
        r.bold = True; r.font.size = Pt(10); r.font.color.rgb = NAVY
        p2 = c.add_paragraph(); p2.paragraph_format.space_after = Pt(2)
        _runs(p2, text, size=9.5)
        self.doc.add_paragraph().paragraph_format.space_after = Pt(4)
        return t

    def figure(self, fname, caption, width=6.4):
        path = os.path.join(self.fig_dir, fname)
        if not os.path.exists(path):
            raise FileNotFoundError(path)
        self.doc.add_picture(path, width=Inches(width))
        self.doc.paragraphs[-1].alignment = CENTER
        p = self.doc.add_paragraph(); p.alignment = CENTER
        p.paragraph_format.space_after = Pt(12)
        for i, chunk in enumerate(caption.split('**')):
            if not chunk:
                continue
            r = p.add_run(chunk)
            r.font.size = Pt(8.5); r.font.color.rgb = GREY
            r.italic = (i % 2 == 0); r.bold = (i % 2 == 1)

    def pagebreak(self):
        self.doc.add_paragraph().add_run().add_break(WD_BREAK.PAGE)
