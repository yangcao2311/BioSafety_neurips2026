# -*- coding: utf-8 -*-
"""Deterministically regenerate figure2.pptx (9-row main text) and
figure_4.pptx (14-row appendix) from one explicit data table, so the two
figures cannot disagree.  Colour is DERIVED from the value by the legend
rule, never hand-set."""
from pptx import Presentation
from pptx.util import Emu, Pt
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE
from pptx.oxml.ns import qn
import copy, sys

FONT = 'Arial Rounded MT Bold'
INK, WHITE, MUTED, ACCENT = '0F172A', 'FFFFFF', '475569', '7C3AED'
BORDER = 'CBD5E1'
PASS, MILD, MOD, HIGH = 'C6EFCE', 'FFEB9C', 'F4B942', 'FF6666'

LEGEND = [(PASS, 'Pass (<0.50)'), (MILD, 'Mild (0.50–0.65)'),
          (MOD, 'Moderate (0.65–0.80)'), (HIGH, 'High (>0.80)'),
          (None, 'Not evaluated')]

COLS = ['L1', 'L2', 'L3', 'L4', 'L5', 'L6a', 'L6b', 'CRISPR', 'siRNA']
NA = 7          # index where Pipeline B starts

# ---------------------------------------------------------------- data
# value: number used ONLY to pick the tier; label: what is printed.
# src  : where the number comes from (kept for the audit trail).
def C(label, value, src):
    return dict(label=label, value=value, src=src)

FLAG = lambda src: C('flag', 1.0, src)

ROWS = [
 ('Genome', 'Evo-2 7b_base', {
    'L1': C('27.8%', 0.278, 'v10_new.tex:468'),
    'L3': C('18.5%', 0.185, 'v10_new.tex:468')}),
 ('Genome', 'NT-v2', {
    'L1': C('24%',   0.24,  'original figure (unsourced in tex)')}),
 ('RNA', 'RNA-FM', {
    'L2': C('54.5%', 0.545, 'v10_new.tex:339,752'),
    'siRNA': C('flag ★', 1.0, 'v10_new.tex:767 T_siRNA(RNA-FM) CD=2')}),
 ('RNA', 'RiNALMo', {
    'L2': C('39%',   0.39,  'original figure (unsourced in tex)'),
    'siRNA': C('0.61', 0.61, 'original figure (unsourced in tex)')}),
 ('Protein', 'ESM-2', {
    'L3': C('0.81',  0.81,  'original figure (unsourced in tex)'),
    'L4': C('r=0.70', 0.70, 'original figure2 (fig4 said -1.6; unified)'),
    'L5': C('<0.50', 0.49,  'original figure: Pass cell, no numeric value')}),
 ('Protein', 'ProtBERT', {
    'L3': C('54.0%', 0.540, 'v10_new.tex:1053')}),
 ('Protein', 'ProtT5', {
    'L3': C('1.7%',  0.017, 'v10_new.tex:1053')}),
 ('Protein', 'ProteinMPNN', {
    'L3': C('r=0.95', 0.954, 'v10_new.tex:788 fig3 caption (NA H275Y)'),
    'L4': C('r=0.70', 0.70,  'v10_new.tex:788'),
    'L5': C('r=0.66', 0.66,  'v10_new.tex:788'),
    'L6b': C('r=0.875', 0.875, 'v10_new.tex:788')}),
 ('Protein', 'ESM-IF1', {
    'L3': C('r=0.51', 0.51, 'original figure (unsourced in tex)'),
    'L4': C('r=0.70', 0.70, 'original figure (unsourced in tex)'),
    'L5': C('r=0.68', 0.68, 'original figure (unsourced in tex)'),
    'L6b': C('r=0.75', 0.75, 'original figure (unsourced in tex)')}),
 ('Chemical', 'MolFormer-XL', {
    'L4': C('0.81',  0.81,  'original figure (unsourced in tex)')}),
 ('Chemical', 'REINVENT 4', {
    'L4': FLAG('v10_new.tex:760 T4.2 L4 flag')}),
 ('Chemical', 'Pocket2Mol', {
    'L4': FLAG('v10_new.tex:762 T4.4 L4 binding flag')}),
 ('Pathway / LLM', 'Geneformer', {
    'L5':  C('r=0.85', 0.85, 'v10_new.tex:765 T5.1'),
    'L6b': C('r=0.80', 0.80, 'v10_new.tex:765 T5.1')}),
 ('Pathway / LLM', 'GPT-4o', {
    'L3': FLAG('v10_new.tex:756 T3.3 GPT-4o CD=4'),
    'L4': FLAG('v10_new.tex:759 T4.3 L4 binding flag'),
    'CRISPR': C('48/604 ★', 1.0, 'v10_new.tex:773 T_CRISPR')}),
]

FIG2_ROWS = ['Evo-2 7b_base', 'RNA-FM', 'ESM-2', 'ProtBERT', 'ProteinMPNN',
             'ESM-IF1', 'REINVENT 4', 'Geneformer', 'GPT-4o']

def tier(v):
    if v < 0.50:  return PASS
    if v < 0.65:  return MILD
    if v <= 0.80: return MOD
    return HIGH

# ---------------------------------------------------------------- helpers
def hatch(shape):
    """Replace solid fill with a light diagonal pattern = 'not evaluated'."""
    spPr = shape.fill._xPr
    for tagname in ('a:solidFill', 'a:noFill', 'a:pattFill', 'a:gradFill'):
        for e in spPr.findall(qn(tagname)):
            spPr.remove(e)
    patt = spPr.makeelement(qn('a:pattFill'), {'prst': 'ltUpDiag'})
    fg = patt.makeelement(qn('a:fgClr'), {})
    fg.append(fg.makeelement(qn('a:srgbClr'), {'val': 'DCE3EC'}))
    bg = patt.makeelement(qn('a:bgClr'), {})
    bg.append(bg.makeelement(qn('a:srgbClr'), {'val': 'FFFFFF'}))
    patt.append(fg); patt.append(bg)
    ln = spPr.find(qn('a:ln'))
    spPr.insert(list(spPr).index(ln) if ln is not None else len(spPr), patt)

def box(slide, x, y, w, h, text='', size=7.7, color=INK, align=PP_ALIGN.CENTER,
        bold=True, anchor=MSO_ANCHOR.MIDDLE):
    tb = slide.shapes.add_textbox(Emu(x), Emu(y), Emu(w), Emu(h))
    tf = tb.text_frame
    tf.word_wrap = True
    tf.vertical_anchor = anchor
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    for i, line in enumerate(text.split('\n')):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = align
        r = p.add_run(); r.text = line
        r.font.name, r.font.size, r.font.bold = FONT, Pt(size), bold
        r.font.color.rgb = RGBColor.from_string(color)
    return tb

def cell(slide, x, y, w, h, fill, text):
    sh = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Emu(x), Emu(y), Emu(w), Emu(h))
    sh.shadow.inherit = False
    sh.line.color.rgb = RGBColor.from_string(BORDER)
    sh.line.width = Emu(4445)
    if fill is None:
        sh.fill.solid(); sh.fill.fore_color.rgb = RGBColor.from_string('FFFFFF')
        hatch(sh)
    else:
        sh.fill.solid(); sh.fill.fore_color.rgb = RGBColor.from_string(fill)
    tf = sh.text_frame
    tf.word_wrap = True
    tf.vertical_anchor = MSO_ANCHOR.MIDDLE
    tf.margin_left = tf.margin_right = Emu(18288)
    tf.margin_top = tf.margin_bottom = 0
    p = tf.paragraphs[0]; p.alignment = PP_ALIGN.CENTER
    r = p.add_run(); r.text = text
    r.font.name, r.font.size, r.font.bold = FONT, Pt(7.7), True
    r.font.color.rgb = RGBColor.from_string(WHITE if fill == HIGH else INK)
    return sh

def dashed_vline(slide, x, y0, y1):
    ln = slide.shapes.add_connector(1, Emu(x), Emu(y0), Emu(x), Emu(y1))
    ln.line.color.rgb = RGBColor.from_string('334155')
    ln.line.width = Emu(9525)
    ln.line._get_or_add_ln().append(
        ln.line._get_or_add_ln().makeelement(qn('a:prstDash'), {'val': 'dash'}))
    return ln

def hline(slide, x0, x1, y, color='94A3B8', w=6350):
    ln = slide.shapes.add_connector(1, Emu(x0), Emu(y), Emu(x1), Emu(y))
    ln.line.color.rgb = RGBColor.from_string(color)
    ln.line.width = Emu(w)
    return ln

# ---------------------------------------------------------------- builder
def build(path, rows, modality_col, slide_w, slide_h, grid_x, label_x, label_w,
          mod_x=None, mod_w=None):
    CW, SX = 820000, 830000
    CH, SY = 400000, 410000
    GAPB   = 190000                      # gap before the Pipeline B block
    Y0     = 980000

    colx = []
    for i in range(len(COLS)):
        colx.append(grid_x + i * SX + (GAPB if i >= NA else 0))
    grid_r = colx[-1] + CW
    rule_x = colx[NA] - GAPB // 2

    prs = Presentation()
    prs.slide_width, prs.slide_height = Emu(slide_w), Emu(slide_h)
    slide = prs.slides.add_slide(prs.slide_layouts[6])

    # --- pipeline banners -------------------------------------------------
    box(slide, colx[0], 430000, colx[NA - 1] + CW - colx[0], 230000,
        'Pipeline A', size=10.5, color=ACCENT)
    box(slide, colx[NA], 300000, grid_r - colx[NA], 230000,
        'Pipeline B', size=10.5, color=ACCENT)
    box(slide, colx[NA] - 130000, 530000, grid_r - colx[NA] + 260000, 300000,
        'task-specific weights\nnot comparable to A', size=6.6, color=ACCENT)

    # --- column headers ---------------------------------------------------
    for i, name in enumerate(COLS):
        box(slide, colx[i], 690000, CW, 230000,
            'T_' + name if i >= NA else name, size=9.5)

    # --- cells ------------------------------------------------------------
    order = {n: k for k, n in enumerate(COLS)}
    for r, (mod, model, data) in enumerate(rows):
        y = Y0 + r * SY
        box(slide, label_x, y, label_w, CH, model, size=8.5,
            align=PP_ALIGN.RIGHT)
        for i, cname in enumerate(COLS):
            d = data.get(cname)
            if d is None:
                cell(slide, colx[i], y, CW, CH, None, '')
            else:
                cell(slide, colx[i], y, CW, CH, tier(d['value']), d['label'])

    bottom = Y0 + len(rows) * SY

    # --- modality grouping (appendix figure only) -------------------------
    if modality_col:
        groups, cur = [], None
        for r, (mod, model, _) in enumerate(rows):
            if mod != cur:
                groups.append([mod, r, r]); cur = mod
            else:
                groups[-1][2] = r
        for mod, a, b in groups:
            box(slide, mod_x, Y0 + a * SY, mod_w, (b - a + 1) * SY - (SY - CH),
                mod, size=8.5, color=MUTED, align=PP_ALIGN.LEFT)
            if a:
                hline(slide, mod_x, grid_r, Y0 + a * SY - (SY - CH) // 2)

    dashed_vline(slide, rule_x, 690000, bottom - (SY - CH))

    # --- legend -----------------------------------------------------------
    lx = grid_r + 330000
    ly = Y0 + 60000
    for k, (col, text) in enumerate(LEGEND):
        sw = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Emu(lx),
                                    Emu(ly + k * 300000), Emu(220000), Emu(190000))
        sw.shadow.inherit = False
        sw.line.color.rgb = RGBColor.from_string(BORDER); sw.line.width = Emu(4445)
        sw.fill.solid(); sw.fill.fore_color.rgb = RGBColor.from_string(col or 'FFFFFF')
        if col is None:
            hatch(sw)
        box(slide, lx + 300000, ly + k * 300000 - 20000, 1600000, 230000,
            text, size=8.2, align=PP_ALIGN.LEFT)
    box(slide, lx - 20000, ly + 5 * 300000 + 150000, 1900000, 640000,
        '%: identity/divergence\n\n0–1: normalized r_L\n\n★: max in Pipeline B',
        size=8.2, color=MUTED, align=PP_ALIGN.LEFT, anchor=MSO_ANCHOR.TOP)

    prs.save(path)
    return grid_r, bottom

if __name__ == '__main__':
    by_name = {m: (mod, m, d) for mod, m, d in ROWS}
    rows2 = [by_name[n] for n in FIG2_ROWS]

    g, b = build('figure2.pptx', rows2, False,
                 11450000, 4950000, grid_x=1500000,
                 label_x=150000, label_w=1270000)
    print('figure2  grid_right=%d bottom=%d' % (g, b))

    g, b = build('figure_4.pptx', ROWS, True,
                 12300000, 6800000, grid_x=2350000,
                 label_x=1070000, label_w=1200000,
                 mod_x=170000, mod_w=860000)
    print('figure_4 grid_right=%d bottom=%d' % (g, b))
