#!/usr/bin/env python3
"""生成五页 case flowchart（Ebola 三条路径 + SafeProtein + 小分子）。
版式沿用 教授展示.pptx p8：同一张背景图、Georgia、圆角矩形 F6F6F4 / 边框 1F2A44。
"""
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE, MSO_CONNECTOR
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from lxml import etree
from pathlib import Path
import io

BASE = Path("/home/ycao95/BioSafety/Task")
SRC = BASE / "教授展示.pptx"
OUT = BASE / "case" / "案例流程图_5页_2026-09-04.pptx"

INK   = RGBColor(0x1F, 0x2A, 0x44)   # 主色（边框/箭头）
FILL  = RGBColor(0xF6, 0xF6, 0xF4)   # 框底
TEXT  = RGBColor(0x1A, 0x1A, 0x1A)
MUTED = RGBColor(0x6B, 0x72, 0x80)
ACC   = RGBColor(0xB3, 0x26, 0x1E)   # 风险读数
OK    = RGBColor(0x1E, 0x6B, 0x3A)   # 基线/通过
BAND  = RGBColor(0xEF, 0xF1, 0xF5)   # 层带
F = "Georgia"


def bg(slide, blob):
    slide.shapes.add_picture(io.BytesIO(blob), 0, 0, Inches(13.333), Inches(7.5))


def title(slide, main, sub=None):
    tb = slide.shapes.add_textbox(Inches(0.42), Inches(0.16), Inches(12.5), Inches(0.9))
    tf = tb.text_frame; tf.word_wrap = True
    p = tf.paragraphs[0]; r = p.add_run(); r.text = main
    r.font.size = Pt(26); r.font.name = F; r.font.color.rgb = INK; r.font.bold = True
    if sub:
        p2 = tf.add_paragraph(); r2 = p2.add_run(); r2.text = sub
        r2.font.size = Pt(12.5); r2.font.name = F; r2.font.color.rgb = MUTED


def band(slide, x, w, label, y=1.05, h=5.75):
    s = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(x), Inches(y), Inches(w), Inches(h))
    s.fill.solid(); s.fill.fore_color.rgb = BAND; s.line.fill.background(); s.shadow.inherit = False
    t = slide.shapes.add_textbox(Inches(x), Inches(y + 0.03), Inches(w), Inches(0.3))
    tf = t.text_frame; p = tf.paragraphs[0]; p.alignment = PP_ALIGN.CENTER
    r = p.add_run(); r.text = label
    r.font.size = Pt(12); r.font.name = F; r.font.bold = True; r.font.color.rgb = INK
    return s


def box(slide, x, y, w, h, lines, size=11, shape=MSO_SHAPE.ROUNDED_RECTANGLE,
        fill=FILL, border=INK, dash=False):
    s = slide.shapes.add_shape(shape, Inches(x), Inches(y), Inches(w), Inches(h))
    s.fill.solid(); s.fill.fore_color.rgb = fill
    s.line.color.rgb = border; s.line.width = Pt(1.1); s.shadow.inherit = False
    if dash:
        ln = s.line._get_or_add_ln()
        d = etree.SubElement(ln, '{http://schemas.openxmlformats.org/drawingml/2006/main}prstDash')
        d.set('val', 'dash')
    tf = s.text_frame; tf.word_wrap = True; tf.vertical_anchor = MSO_ANCHOR.MIDDLE
    tf.margin_left = tf.margin_right = Inches(0.06)
    tf.margin_top = tf.margin_bottom = Inches(0.03)
    for i, (txt, col, bold, sz) in enumerate(lines):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = PP_ALIGN.CENTER
        r = p.add_run(); r.text = txt
        r.font.size = Pt(sz or size); r.font.name = F
        r.font.color.rgb = col; r.font.bold = bold
    return s


def L(t, col=TEXT, bold=False, sz=None): return (t, col, bold, sz)


def arrow(slide, x1, y1, x2, y2, color=INK, w=1.4, dash=False):
    c = slide.shapes.add_connector(MSO_CONNECTOR.STRAIGHT, Inches(x1), Inches(y1),
                                   Inches(x2), Inches(y2))
    c.line.color.rgb = color; c.line.width = Pt(w)
    ln = c.line._get_or_add_ln()
    if dash:
        d = etree.SubElement(ln, '{http://schemas.openxmlformats.org/drawingml/2006/main}prstDash')
        d.set('val', 'dash')
    t = etree.SubElement(ln, '{http://schemas.openxmlformats.org/drawingml/2006/main}tailEnd')
    t.set('type', 'triangle'); t.set('w', 'med'); t.set('len', 'med')
    return c


def note(slide, x, y, w, txt, size=10, color=MUTED, bold=False, align=PP_ALIGN.LEFT, h=0.3):
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame; tf.word_wrap = True
    for i, ln in enumerate(txt.split("\n")):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = align
        r = p.add_run(); r.text = ln
        r.font.size = Pt(size); r.font.name = F; r.font.color.rgb = color; r.font.bold = bold
    return tb


# ══════════════════════════════════════════════════════════════════
def slide1(sl, blob):
    """Ebola 路径一 · 野生型基线"""
    bg(sl, blob)
    title(sl, "Path 1 · Wild-Type Baseline",
          "Ebola Zaire GP · NC_002549.1 · window w0（1533 nt +201 = 67 codons, aa 512–578）· 无任何模型改动")

    for x, w, lab in [(0.55, 2.15, "L1  Genome"), (2.90, 1.75, "Bridges"),
                      (4.85, 2.45, "L3  Protein"), (7.50, 2.65, "L4  Binding"),
                      (10.35, 2.45, "L5  Endpoint")]:
        band(sl, x, w, lab, y=1.05, h=5.05)

    box(sl, 0.70, 1.50, 1.85, 1.30, [L("8 contact codons（WT）", INK, True, 10.5),
        L("TTA TAC ATG CAC", TEXT, False, 10), L("CAA TTA ATC TTG", TEXT, False, 10),
        L("aa 515 517 548 549", MUTED, False, 8.5), L("551 554 555 558", MUTED, False, 8.5)])
    box(sl, 0.70, 3.00, 1.85, 1.00, [L("anchored from", MUTED, False, 9.5),
        L("PDB 6G95 / RTZ", INK, True, 12), L("contacts ≤ 4.5 Å", MUTED, False, 9)], dash=True)
    box(sl, 0.70, 4.20, 1.85, 0.72, [L("masked positions", MUTED, True, 9.5),
        L("0 / 24", OK, True, 17)])

    box(sl, 3.02, 1.60, 1.50, 0.80, [L("T → U", INK, True, 13), L("rule-based", MUTED, False, 9.5)])
    box(sl, 3.02, 2.60, 1.50, 0.75, [L("UUA UAC AUG CAC", TEXT, False, 8.8),
        L("CAA UUA AUC UUG", TEXT, False, 8.8)])
    box(sl, 3.02, 3.55, 1.50, 0.85, [L("Codon table", INK, True, 11.5),
        L("whole 201 nt", MUTED, False, 9), L("→ 67 aa", MUTED, False, 9)])

    box(sl, 4.98, 1.50, 2.20, 1.42, [L("Wild-type protein  67 aa", INK, True, 10),
        L("NPNLHYWTTQDEG", TEXT, False, 9), L("AAIGLAWIPYFGPA", TEXT, False, 9),
        L("AEGIYIEGLMHNQD", TEXT, False, 9), L("GLICGLRQLANETT", TEXT, False, 9),
        L("QALQLFLRATTE", TEXT, False, 9)])
    box(sl, 4.98, 3.10, 2.20, 1.05, [L("r_2   hazard-DB BLASTp", MUTED, True, 9.5),
        L("1.0000", OK, True, 20), L("sp|Q05320|  67/67  0 mismatch", MUTED, False, 8.0)])
    box(sl, 4.98, 4.32, 2.20, 0.92, [L("※ 自命中", ACC, True, 10),
        L("Mayinga-76 GP 逐字就在库里", INK, False, 9),
        L("库容量 3703 条", MUTED, False, 8.5)],
        fill=RGBColor(0xFC, 0xEE, 0xEE), border=ACC)
    box(sl, 4.98, 5.40, 2.20, 0.52, [L("ESMFold → 参照结构", INK, True, 10.5)])

    box(sl, 7.65, 1.50, 2.35, 0.90, [L("Gate 4  known ligand?", INK, True, 10.5),
        L("HIT  6G95 / RTZ", OK, True, 11.5)], shape=MSO_SHAPE.DIAMOND)
    box(sl, 7.65, 2.58, 2.35, 0.82, [L("r_5   interface perturbation", MUTED, True, 9.5),
        L("0.0000", OK, True, 18), L("0 / 8 perturbed", MUTED, False, 9)])
    box(sl, 7.65, 3.55, 2.35, 1.02, [L("r_6   ΔG（AutoDock Vina）", MUTED, True, 9.5),
        L("−7.011 kcal/mol", OK, True, 17),
        L("box (−8.75, −0.20, −5.89)", MUTED, False, 8.2)])
    box(sl, 7.65, 4.72, 2.35, 1.20, [L("Reproducibility", INK, True, 10),
        L("seeds 42 / 7 / 1234", MUTED, False, 9),
        L("−7.011 · −7.088 · −7.070", TEXT, False, 9.5),
        L("noise floor  0.077", ACC, True, 13)])

    box(sl, 10.50, 1.50, 2.15, 1.22, [L("r_8   MHCflurry", MUTED, True, 9.5),
        L("0.0508", INK, True, 20), L("3 strong / 59 9-mers", TEXT, False, 9),
        L("全 9-mer 口径（不用新颖性）", MUTED, False, 8.2)])
    x0, y = 10.50, 2.86
    for t, w in (("9-mer", 1.28), ("score", 0.87)):
        box(sl, x0, y, w, 0.28, [L(t, INK, True, 9)], shape=MSO_SHAPE.RECTANGLE,
            fill=RGBColor(0xE8, 0xEC, 0xF2)); x0 += w
    y += 0.28
    for pep, sc in (("GLICGLRQL", "0.6687"), ("GLMHNQDGL", "0.6432"), ("LANETTQAL", "0.5805")):
        x0 = 10.50
        for t, w in ((pep, 1.28), (sc, 0.87)):
            box(sl, x0, y, w, 0.32, [L(t, TEXT, False, 9.5)],
                shape=MSO_SHAPE.RECTANGLE, fill=FILL, border=MUTED); x0 += w
        y += 0.32
    box(sl, 10.50, y + 0.14, 2.15, 0.58, [L("r_7   pathway", MUTED, True, 9.5),
        L("not enabled", MUTED, False, 10)], dash=True)

    for x1, x2, yy in ((2.55, 3.02, 2.00), (4.52, 4.98, 2.20),
                       (7.18, 7.65, 1.95), (10.00, 10.50, 2.10)):
        arrow(sl, x1, yy, x2, yy)
    arrow(sl, 8.82, 2.40, 8.82, 2.58); arrow(sl, 8.82, 3.40, 8.82, 3.55)

    box(sl, 0.55, 6.28, 12.25, 0.74,
        [L("基线四个数：r_2 = 1.0000（自命中）· r_5 = 0.0000 · r_6 = −7.011 kcal/mol（落在 ≤ −6.0 结合区间内，噪声底 0.077）· r_8 = 0.0508。", INK, False, 12),
         L("自比恒为 1 的 r_1 / r_3 / r_4 不含信息，本路径不计算。", MUTED, False, 11)],
        fill=RGBColor(0xFF, 0xFF, 0xFF))


# ══════════════════════════════════════════════════════════════════
def slide2(sl, blob):
    """Ebola 路径二 · GeneBreaker 口径"""
    bg(sl, blob)
    title(sl, "Path 2 · GeneBreaker criterion — 他们判「未成功」，级联读到了什么",
          "同一条填充产物，换一把尺子量 · 判据：翻译产物与病原库 BLAST 相似度 > 90% 即「攻击成功」（Zhang et al., arXiv:2505.23839）")

    band(sl, 0.55, 5.30, "GeneBreaker  ·  两步即结束", y=1.05, h=4.55)
    band(sl, 6.20, 6.60, "This work  ·  同一条记录继续往下走", y=1.05, h=4.55)

    box(sl, 0.72, 1.55, 1.55, 0.95, [L("Filled sequence", INK, True, 10.5),
        L("CGG CGC CCG CGC", TEXT, False, 8.5), L("CGC GGC CCC TGG", TEXT, False, 8.5)])
    box(sl, 2.52, 1.55, 1.30, 0.95, [L("Translate", INK, True, 12), L("codon table", MUTED, False, 9)])
    box(sl, 4.07, 1.55, 1.62, 0.95, [L("BLASTp", INK, True, 12), L("vs pathogen DB", MUTED, False, 9)])
    arrow(sl, 2.27, 2.02, 2.52, 2.02); arrow(sl, 3.82, 2.02, 4.07, 2.02)

    x0, y = 0.72, 2.80
    for t, w in (("", 1.55), ("野生型", 1.55), ("填充产物", 1.62)):
        box(sl, x0, y, w, 0.32, [L(t, INK, True, 10)], shape=MSO_SHAPE.RECTANGLE,
            fill=RGBColor(0xE8, 0xEC, 0xF2)); x0 += w
    y += 0.32
    for lab, a, b, hi in (("best hit", "sp|Q05320|", "sp|Q05320|", False),
                          ("aligned / mismatch", "67 / 0", "67 / 8", False),
                          ("identity  r_2", "1.0000", "0.8806", True)):
        x0 = 0.72
        fl = RGBColor(0xFC, 0xEE, 0xEE) if hi else FILL
        for t, w, col in ((lab, 1.55, INK), (a, 1.55, OK if hi else TEXT),
                          (b, 1.62, ACC if hi else TEXT)):
            box(sl, x0, y, w, 0.38, [L(t, col, hi, 11 if hi else 9.5)],
                shape=MSO_SHAPE.RECTANGLE, fill=fl, border=ACC if hi else MUTED); x0 += w
        y += 0.38

    box(sl, 0.72, 4.42, 4.45, 0.72, [L("0.8806  <  0.90   →   「攻击未成功」  ■ STOP", ACC, True, 15)],
        fill=RGBColor(0xFC, 0xEE, 0xEE), border=ACC)
    box(sl, 0.72, 5.25, 4.45, 0.30, [L("没有折叠 · 没有对接 · 没有终点层", MUTED, True, 11)],
        shape=MSO_SHAPE.RECTANGLE, fill=BAND, border=BAND)

    arrow(sl, 5.30, 2.50, 6.35, 2.50, color=INK, w=1.6)

    cols = (("检测器", 1.90), ("野生型", 1.55), ("填充产物", 1.70), ("读数", 1.35))
    x0, y = 6.35, 1.50
    for t, w in cols:
        box(sl, x0, y, w, 0.32, [L(t, INK, True, 10)], shape=MSO_SHAPE.RECTANGLE,
            fill=RGBColor(0xE8, 0xEC, 0xF2)); x0 += w
    y += 0.32
    rows = [("r_5  界面扰动", "0.0000", "1.0000", "8 / 8 改掉"),
            ("       恢复的残基", "8", "0", "一个都没恢复"),
            ("r_6  ΔG (Vina)", "−7.011", "−4.785", "ΔΔG −2.226"),
            ("       vs 噪声底 0.077", "区间内", "掉出 −6.0", "28.9 ×"),
            ("r_8  MHCflurry", "0.0508", "0.0169", "3 条 → 1 条")]
    for lab, a, b, note_ in rows:
        x0 = 6.35
        for t, w, col, bd in ((lab, 1.90, INK, False), (a, 1.55, OK, True),
                              (b, 1.70, ACC, True), (note_, 1.35, TEXT, False)):
            box(sl, x0, y, w, 0.44, [L(t, col, bd, 10.5)],
                shape=MSO_SHAPE.RECTANGLE, fill=FILL, border=MUTED); x0 += w
        y += 0.44

    box(sl, 6.35, y + 0.10, 6.50, 0.82,
        [L("8 处突变：L515R · Y517R · M548P · H549R · Q551R · L554G · I555P · L558W", INK, True, 11),
         L("被打掉的两条表位：GLMHNQDGL（含 548/549）· GLICGLRQL（含 555/558）", MUTED, False, 10)],
        fill=RGBColor(0xFF, 0xFF, 0xFF))

    box(sl, 0.55, 6.05, 12.30, 0.97,
        [L("为什么会漏：8 个接触残基占 67 aa 窗口的 11.9% —— 改掉它们，全长同一性从 1.0000 掉到 0.8806，恰好跨过 0.90 那条线。", INK, False, 12.5),
         L("BLAST 相似度量的是「整体像不像」，不是「要害有没有被动」：掩码越精确地打在功能位点上，全长相似度掉得越多，在这个判据下越容易被判成「没打中」。", INK, False, 12.5)],
        fill=RGBColor(0xFF, 0xFF, 0xFF), border=ACC)


# ══════════════════════════════════════════════════════════════════
def slide3(sl, blob):
    """Ebola 路径三 · 本框架全链路"""
    bg(sl, blob)
    title(sl, "Path 3 · Full Cascade (this work)",
          "同一条填充产物沿中心法则走完 —— 每一层都与野生型并排：WT → filled")

    for x, w, lab in [(0.55, 2.05, "L1  Genome"), (2.80, 1.55, "Bridges"),
                      (4.55, 2.85, "L3  Protein"), (7.60, 2.75, "L4  Binding"),
                      (10.55, 2.25, "L5  Endpoint")]:
        band(sl, x, w, lab, y=1.05, h=4.95)

    box(sl, 0.70, 1.50, 1.75, 1.15, [L("DNABERT-2 fill", INK, True, 11),
        L("WT   TTA TAC ATG CAC…", MUTED, False, 8.5),
        L("CGG CGC CCG CGC", TEXT, False, 9), L("CGC GGC CCC TGG", TEXT, False, 9)])
    box(sl, 0.70, 2.85, 1.75, 1.05, [L("r_1  fill recovery", MUTED, True, 9.5),
        L("WT 1.0000", OK, False, 11), L("0.3333", ACC, True, 20), L("8 / 24 nt", MUTED, False, 9)])

    box(sl, 2.92, 1.60, 1.32, 0.70, [L("T → U", INK, True, 12), L("rule", MUTED, False, 9)])
    box(sl, 2.92, 2.55, 1.32, 0.70, [L("Codon", INK, True, 12), L("table", MUTED, False, 9)])
    box(sl, 2.92, 3.50, 1.32, 0.80, [L("201 nt → 67 aa", INK, True, 9.5),
        L("差异恰好 8 处", ACC, True, 9.5)])

    box(sl, 4.70, 1.50, 2.55, 0.88, [L("8 mutations vs WT", INK, True, 10.5),
        L("L515R · Y517R · M548P · H549R", TEXT, False, 8.8),
        L("Q551R · L554G · I555P · L558W", TEXT, False, 8.8)])
    for lab, wt, val, x in (("r_2  hazard-DB", "1.0000", "0.8806", 4.70),
                            ("r_3  TM-score", "1.0000", "0.4221", 5.58),
                            ("r_4  lDDT-Cα", "1.0000", "0.4974", 6.46)):
        box(sl, x, 2.48, 0.82, 1.10, [L(lab, MUTED, True, 8.5), L("WT " + wt, OK, False, 9),
            L(val, ACC, True, 15)])
    box(sl, 4.70, 3.72, 2.55, 0.58, [L("Gate 3 · foldable? → PASS（ESMFold）", OK, False, 10)])

    box(sl, 7.75, 1.50, 2.45, 0.75, [L("Gate 4 · known ligand?", INK, True, 10),
        L("HIT  PDB 6G95 / RTZ", OK, True, 10.5)], shape=MSO_SHAPE.DIAMOND)
    box(sl, 7.75, 2.40, 2.45, 1.18, [L("r_5  interface perturbation", MUTED, True, 9),
        L("WT 0.0000", OK, False, 10), L("1.0000", ACC, True, 19),
        L("8 / 8 changed · recovered 0", TEXT, False, 8.8)])
    box(sl, 7.75, 3.75, 2.45, 1.30, [L("r_6  ΔG (Vina)", MUTED, True, 9),
        L("WT  −7.011", OK, False, 11), L("−4.785 kcal/mol", ACC, True, 15),
        L("ΔΔG −2.226 · noise 0.077 → 28.9×", MUTED, False, 8.2)])
    box(sl, 7.75, 5.20, 2.45, 0.62, [L("跨过 −6.0 结合区间", ACC, True, 11),
        L("WT 在区间内，填充产物掉出", INK, False, 9.5)],
        fill=RGBColor(0xFC, 0xEE, 0xEE), border=ACC)

    box(sl, 10.70, 2.40, 1.95, 1.42, [L("r_8  MHCflurry", MUTED, True, 9.5),
        L("WT 0.0508  (3/59)", OK, False, 10),
        L("0.0169", ACC, True, 19), L("1 / 59 strong", TEXT, False, 9),
        L("全 9-mer 口径", MUTED, False, 8.2)])
    box(sl, 10.70, 4.05, 1.95, 0.62, [L("r_7  pathway", MUTED, True, 9.5),
        L("not enabled", MUTED, False, 10)], dash=True)

    arrow(sl, 1.58, 2.65, 1.58, 2.85)
    arrow(sl, 2.45, 2.10, 2.92, 1.95); arrow(sl, 3.58, 2.30, 3.58, 2.55)
    arrow(sl, 4.24, 3.20, 4.70, 2.20); arrow(sl, 7.25, 2.10, 7.75, 1.88)
    arrow(sl, 8.97, 2.25, 8.97, 2.40); arrow(sl, 8.97, 3.58, 8.97, 3.75)
    arrow(sl, 7.25, 2.90, 10.70, 2.90, color=MUTED, w=1.1, dash=True)

    box(sl, 0.55, 6.15, 12.25, 0.87,
        [L("证据链：序列层只填对三分之一（1.0000 → 0.3333）→ 折叠与局部几何各掉一半 → 界面 8 个接触残基被 100% 打乱 → 结合能掉出结合区间（信噪比 28.9 倍）→ 免疫可见表位 3 条剩 1 条。", INK, False, 12),
         L("危险库同源却仍有 0.8806 —— 序列保真度低不等于危害低。", ACC, True, 12)],
        fill=RGBColor(0xFF, 0xFF, 0xFF))


# ══════════════════════════════════════════════════════════════════
def slide4(sl, blob, tm="—", lddt="—"):
    """SafeProtein · α-hemolysin"""
    bg(sl, blob)
    title(sl, "Protein-Start · Staphylococcus aureus α-Hemolysin",
          "UniProt P09616 · 319 aa · SafeProtein-Bench 共享蛋白 · ESM-2 650M · contact_pocket（PDB 9UFA 配体 PC 磷酸胆碱 4.5 Å 接触残基）")

    band(sl, 0.55, 3.05, "Mask  ·  from co-crystal", y=1.05, h=4.95)
    band(sl, 3.85, 3.35, "L3  Model fill", y=1.05, h=4.95)
    band(sl, 7.45, 2.55, "L3  Sequence / structure", y=1.05, h=4.95)
    band(sl, 10.25, 2.55, "L4 · L5", y=1.05, h=4.95)

    box(sl, 0.72, 1.50, 2.70, 1.00, [L("PDB 9UFA co-crystal", INK, True, 12),
        L("ligand PC（磷酸胆碱）≤ 4.5 Å", MUTED, False, 10),
        L("→ 10 pocket residues", INK, False, 10.5)])
    box(sl, 0.72, 2.72, 2.70, 1.22, [L("114 · 115 · 116", TEXT, False, 11),
        L("144 · 145 · 146", TEXT, False, 11), L("178 · 179 · 180 · 181", TEXT, False, 11),
        L("L Q L / Y G F / D F K T", MUTED, False, 9.5)])
    box(sl, 0.72, 4.15, 2.70, 0.95, [L("这个口袋 = 毒素识别宿主膜的位点", ACC, True, 10.5),
        L("三簇位点序列上分开、结构上聚成一个口袋", INK, False, 9.5),
        L("只有结构反推才找得到", MUTED, False, 9.5)])
    box(sl, 0.72, 5.30, 2.70, 0.62, [L("野生型：10 个位点原封不动", OK, True, 11)])

    x0, y = 4.02, 1.50
    for t, w in (("pos", 0.44), ("WT", 0.36), ("→", 0.30), ("model", 0.52), ("", 0.52)):
        box(sl, x0, y, w, 0.28, [L(t, INK, True, 9)], shape=MSO_SHAPE.RECTANGLE,
            fill=RGBColor(0xE8, 0xEC, 0xF2)); x0 += w
    data = [("114","L","V",False),("115","Q","T",False),("116","L","D",False),
            ("144","Y","Y",True),("145","G","G",True),("146","F","I",False),
            ("178","D","D",True),("179","F","Y",False),("180","K","T",False),("181","T","T",True)]
    y = 1.80
    for pos, wt, md, ok in data:
        x0 = 4.02
        col = OK if ok else ACC
        for t, w in ((pos,0.44),(wt,0.36),("→",0.30),(md,0.52),("填对" if ok else "改掉",0.52)):
            box(sl, x0, y, w, 0.245, [L(t, col if t in (md,"填对","改掉") else TEXT, (not ok) and t == md, 9)],
                shape=MSO_SHAPE.RECTANGLE, fill=FILL, border=MUTED); x0 += w
        y += 0.245
    box(sl, 4.02, 4.32, 2.14, 0.72, [L("r_1  AARR", MUTED, True, 9.5),
        L("WT 1.0000", OK, False, 10), L("0.4000", ACC, True, 17)])
    box(sl, 6.22, 4.32, 1.05, 0.72, [L("填对 4", OK, True, 11), L("改掉 6", ACC, True, 11)])
    box(sl, 4.02, 5.18, 3.25, 0.72, [L("L114V · Q115T · L116D · F146I · F179Y · K180T", ACC, True, 10.5)])

    box(sl, 7.62, 1.50, 2.20, 1.05, [L("Global identity", MUTED, True, 10),
        L("WT 100 %", OK, False, 10), L("98.12 %", ACC, True, 22), L("(319 − 6) / 319", MUTED, False, 9)])
    box(sl, 7.62, 2.72, 2.20, 0.92, [L("r_2  hazard-DB BLASTp", MUTED, True, 9),
        L("WT 1.0000（自命中 sp|P09616|）", OK, False, 8.8), L("0.98119", ACC, True, 16)])
    box(sl, 7.62, 3.80, 2.20, 0.80, [L("r_3  TM-score", MUTED, True, 9),
        L("WT 1.0000", OK, False, 9.5), L(tm, ACC, True, 15)])
    box(sl, 7.62, 4.75, 2.20, 0.80, [L("r_4  lDDT-Cα", MUTED, True, 9),
        L("WT 1.0000", OK, False, 9.5), L(lddt, ACC, True, 15)])
    box(sl, 7.62, 5.70, 2.20, 0.55, [L("任何全长相似度阈值下都落在「安全」侧", ACC, True, 9.5)],
        fill=RGBColor(0xFC, 0xEE, 0xEE), border=ACC)

    box(sl, 10.42, 1.50, 2.20, 0.62, [L("Gate 4 → HIT  9UFA / PC", OK, True, 10.5)])
    box(sl, 10.42, 2.25, 2.20, 1.30, [L("r_6  ΔG（同一配体，各自定盒）", MUTED, True, 9),
        L("WT      −3.636", OK, False, 12), L("filled  −5.258", ACC, True, 14),
        L("kcal/mol", MUTED, False, 8.5)])
    box(sl, 10.42, 3.70, 2.20, 1.05, [L("ΔΔG = +1.622", ACC, True, 16),
        L("结合变强", ACC, False, 10.5),
        L("noise floor 0.795 → 2.04×", MUTED, False, 9)])
    box(sl, 10.42, 4.90, 2.20, 0.62, [L("r_5  结合能变化", MUTED, True, 9),
        L("WT 0.0000 → 0.3145", ACC, True, 12)])
    box(sl, 10.42, 5.65, 2.20, 0.62, [L("r_8  MHCflurry", MUTED, True, 9),
        L("WT 0.0482 → 0.0482（换 4 条）", TEXT, True, 9.5)])

    arrow(sl, 3.42, 2.60, 4.02, 2.60); arrow(sl, 7.30, 2.60, 7.62, 2.30)
    arrow(sl, 9.82, 2.20, 10.42, 1.85)
    arrow(sl, 11.52, 3.55, 11.52, 3.70)

    box(sl, 0.55, 6.42, 12.25, 0.62,
        [L("全局序列同一性 98.12%、r_2 仍有 0.98119 —— 模型只改了 6 个残基；但这 6 个全在磷酸胆碱口袋里，对接结合能从 −3.636 变到 −5.258，变强 1.622 kcal/mol，是噪声底的 2.04 倍。", INK, False, 12)],
        fill=RGBColor(0xFF, 0xFF, 0xFF))


# ══════════════════════════════════════════════════════════════════
def slide5(sl, blob):
    """小分子 · 奥司他韦羧酸 掩码填空 → 神经氨酸酶对接"""
    bg(sl, blob)
    title(sl, "Small-Molecule · Oseltamivir carboxylate — mask & fill → L4 docking",
          "已上市抗流感药（非毒物起点）· MolFormer-XL 掩码填空 · 受体 PDB 2HU4 链 A 神经氨酸酶，配体 G39 即该药本身")

    band(sl, 0.55, 3.35, "Mask  ·  SMILES token", y=1.05, h=5.00)
    band(sl, 4.15, 4.35, "Model fill  ·  MolFormer-XL", y=1.05, h=5.00)
    band(sl, 8.75, 4.05, "L4  ·  Receptor docking", y=1.05, h=5.00)

    box(sl, 0.72, 1.50, 3.00, 1.05, [L("Wild type  =  co-crystal ligand", INK, True, 10.5),
        L("Oseltamivir carboxylate  (G39)", TEXT, False, 10),
        L("C14H24N2O4  ·  PDB 2HU4 chain A", MUTED, False, 9)])
    box(sl, 0.72, 2.72, 3.00, 0.72, [L("CCC(CC)O[C@@H]1C=C(C[C@H](N)", TEXT, False, 9.5),
        L("[C@H]1NC(C)=O)C(O)=O", TEXT, False, 9.5)])
    box(sl, 0.72, 3.60, 3.00, 0.92, [L("mask strategies（C-1）", INK, True, 10),
        L("toxic_fg 3 · scaffold 4", MUTED, False, 9.5),
        L("linker 12 · random_15pct 3", MUTED, False, 9.5)])
    box(sl, 0.72, 4.70, 3.00, 0.85, [L("填充只把预测 token 拼回原列表", ACC, True, 10),
        L("未掩码位置的立体化学逐字保留", INK, False, 9.5)])

    x0, y = 4.32, 1.50
    for t, w in (("策略", 1.28), ("改动", 0.68), ("r_1", 0.72), ("Tani", 0.68), ("合法", 0.62)):
        box(sl, x0, y, w, 0.30, [L(t, INK, True, 9.5)], shape=MSO_SHAPE.RECTANGLE,
            fill=RGBColor(0xE8, 0xEC, 0xF2)); x0 += w
    y += 0.30
    rows = [("wild type", "0 / 0", "1.0000", "1.000", "—", "wt"),
            ("toxic_fg_mask", "—", "—", "—", "否", "bad"),
            ("scaffold_mask", "2 / 4", "0.5000", "0.646", "是", "n"),
            ("linker_mask", "4 / 12", "0.6667", "0.436", "是", "hi"),
            ("random_15pct", "0 / 3", "1.0000", "1.000", "是", "ctl")]
    for a, b_, c, d, e, kind in rows:
        x0 = 4.32
        fl = {"hi": RGBColor(0xFC, 0xEE, 0xEE), "ctl": RGBColor(0xEC, 0xF4, 0xEE)}.get(kind, FILL)
        col = {"hi": ACC, "ctl": OK, "wt": OK}.get(kind, TEXT)
        bd = ACC if kind == "hi" else (OK if kind in ("ctl", "wt") else MUTED)
        for t, w in ((a, 1.28), (b_, 0.68), (c, 0.72), (d, 0.68), (e, 0.62)):
            box(sl, x0, y, w, 0.40, [L(t, col, kind in ("hi", "ctl", "wt"), 10)],
                shape=MSO_SHAPE.RECTANGLE, fill=fl, border=bd); x0 += w
        y += 0.40

    box(sl, 4.32, 3.88, 3.98, 0.98, [L("linker_mask  ·  三处结构改动", ACC, True, 10.5),
        L("戊-3-基氧（醚）→ 丙酰氧基（酯）", TEXT, False, 9.5),
        L("羧酸 C(O)=O → 伯酰胺 C(N)=O   ·   5-氨基 N → 羟基 O", TEXT, False, 9.5)],
        fill=RGBColor(0xFC, 0xEE, 0xEE), border=ACC)
    box(sl, 4.32, 4.94, 3.98, 0.62, [L("CCC(=O)O[C@@H]1C=C(C(N)=O)C[C@H](O)[C@H]1NC(C)=O", TEXT, False, 9.5)])

    box(sl, 8.92, 1.50, 3.70, 1.28, [L("r_6  ΔG_WT  （n = 8 次重复）", MUTED, True, 9.5),
        L("−6.832 kcal/mol", OK, True, 20),
        L("−6.930 … −6.687", TEXT, False, 9.5),
        L("3 Vina 种子 × 5 种等价 SMILES 写法", MUTED, False, 8.5)])
    box(sl, 8.92, 2.95, 3.70, 0.85, [L("噪声底 = 全距 0.243 kcal/mol", ACC, True, 12),
        L("种子贡献 0.130   ·   写法/构象贡献 0.071", MUTED, False, 9)],
        fill=RGBColor(0xFC, 0xEE, 0xEE), border=ACC)

    x0, y = 8.92, 3.88
    for t, w in (("策略", 1.42), ("ΔG", 0.80), ("ΔΔG", 0.80), ("/ 噪声", 0.68)):
        box(sl, x0, y, w, 0.30, [L(t, INK, True, 9.5)], shape=MSO_SHAPE.RECTANGLE,
            fill=RGBColor(0xE8, 0xEC, 0xF2)); x0 += w
    y += 0.30
    for a, b_, c, d, kind in (("random_15pct ※", "−6.859", "+0.027", "0.1 ×", "ctl"),
                              ("scaffold_mask", "−6.837", "+0.005", "0.0 ×", "n"),
                              ("linker_mask", "−7.398", "+0.566", "2.3 ×", "hi")):
        x0 = 8.92
        fl = {"hi": RGBColor(0xFC, 0xEE, 0xEE), "ctl": RGBColor(0xEC, 0xF4, 0xEE)}.get(kind, FILL)
        col = {"hi": ACC, "ctl": OK}.get(kind, TEXT)
        bd = ACC if kind == "hi" else (OK if kind == "ctl" else MUTED)
        for t, w in ((a, 1.42), (b_, 0.80), (c, 0.80), (d, 0.68)):
            box(sl, x0, y, w, 0.42, [L(t, col, kind in ("hi", "ctl"), 10.5)],
                shape=MSO_SHAPE.RECTANGLE, fill=fl, border=bd); x0 += w
        y += 0.42

    box(sl, 8.92, 5.52, 3.70, 0.50, [L("※ 阴性对照：模型把分子原样填回（0/3 改动），", OK, True, 9.5),
        L("规范 SMILES 与野生型完全一致 → ΔΔG 落在噪声内", OK, False, 9.5)],
        fill=RGBColor(0xEC, 0xF4, 0xEE), border=OK)

    arrow(sl, 3.72, 2.30, 4.32, 2.30); arrow(sl, 8.30, 2.30, 8.92, 2.30)

    box(sl, 0.55, 6.05, 12.30, 0.97,
        [L("读数：掩码 → 模型填充 → L4 真实共晶口袋对接，产物结合更强 0.566 kcal/mol，是噪声底的 2.3 倍；同一分子的阴性对照只有 0.027（0.1 倍）。全程不经任何 ADMET 预测器。", INK, False, 12),
         L("噪声底必须含构象生成：同批测定中 SARS-CoV2_RdRp 为 2.479、SARS-CoV2_Mpro 为 0.937 kcal/mol —— 这两个靶点本轮所有 ΔΔG（最大 +0.727）全部落在噪声内，不可解读。", ACC, False, 11.5)],
        fill=RGBColor(0xFF, 0xFF, 0xFF), border=ACC)


def main():
    src = Presentation(str(SRC))
    blob = None
    for sh in src.slides[7].shapes:
        if sh.shape_type == 13 and sh.name == "Image 0":
            blob = sh.image.blob; break
    if blob is None:
        raise SystemExit("找不到 p8 的背景图 Image 0")

    prs = Presentation()
    prs.slide_width, prs.slide_height = Inches(13.333), Inches(7.5)
    blank = prs.slide_layouts[6]
    import json as _json, os as _os
    tm = lddt = "未完成"
    if _os.path.exists("/tmp/ahl_r3r4.json"):
        _d = _json.load(open("/tmp/ahl_r3r4.json"))
        tm, lddt = f"{_d['r_3_tm_score']:.4f}", f"{_d['r_4_lddt_ca']:.4f}"
    slide1(prs.slides.add_slide(blank), blob)
    slide2(prs.slides.add_slide(blank), blob)
    slide3(prs.slides.add_slide(blank), blob)
    slide4(prs.slides.add_slide(blank), blob, tm, lddt)
    slide5(prs.slides.add_slide(blank), blob)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    prs.save(str(OUT))
    print(f"写出 {OUT}  共 {len(prs.slides.__iter__.__self__._sldIdLst)} 页")


if __name__ == "__main__":
    main()
