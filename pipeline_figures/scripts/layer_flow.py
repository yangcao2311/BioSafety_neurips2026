#!/usr/bin/env python3
"""层带式流程页渲染器。

版式约定（按用户 2026-09-02 口径）：
  * 顶部 label = **层（layer）**，不是步骤；L1/L2 合并时共用一个颜色、一条带
  * 每个上游 → 产物的转换都必须有**箭头**
  * **只有 bridge / r 指标 / Gate 三类信息进框**，其余一律纯文字
  * 参考 deck 的层义：L5 = Function & Consequence（后果层，含免疫）；**无 L6**；
    L5 里的 pathway 档不启用，不画

图形基元复用 build_four_case_deck_v2 的 S 类。
"""
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE

from build_four_case_deck_v2 import (
    S, rgb, W, H, SANS, MONO,
    C_WT, C_GEN, C_GRAY, C_BLACK, C_WHITE, C_RED, C_GREEN, C_TOOLBG,
    FOOT_Y, FOOT_H,
)

# 层色板：同一层同一色；L1/L2 合并时只取一个
LAYER_COLOR = {
    "L1":  ("1C5A7D", "E8F1F8"),
    "L2":  ("2F8F5B", "E9F6EE"),
    "L12": ("1C5A7D", "E8F1F8"),   # L1/L2 合并 —— 用 L1 的颜色
    "L3":  ("7A4FA3", "F3EDF9"),
    "L4":  ("C9622A", "FCEEE3"),
    "L5":  ("6B5A1A", "FBF7EA"),
    "X":   ("5A5A66", "EFEFEC"),   # 无层号的管线
}

BAND_Y, BAND_H = 1.02, 10.30
LANE1_Y, LANE2_Y = 2.05, 6.05
CHIP_Y = 8.55


def layer_slide(prs, sp):
    sl = S(prs)
    sl.txt(0.4, 0.08, W - 0.8, 0.34, sp["title"], size=17.5, bold=True, color=rgb("12293F"))
    sl.txt(0.4, 0.42, W - 0.8, 0.34, sp["subtitle"], size=9.5, color=rgb("4A4A4A"), ls=1.12)

    layers = sp["layers"]
    tot = sum(l["w"] for l in layers)
    xs, x = [], 0.0
    for l in layers:
        l["w"] = l["w"] * W / tot
        xs.append(x); x += l["w"]

    # 泳道标签（纯色块 + 文字，不算"信息框"）
    sl.rbox(0.10, LANE1_Y - 0.05, 1.05, 1.55, rgb("EEEEEE"), radius=0.2)
    sl.txt(0.10, LANE1_Y - 0.05, 1.05, 1.55, sp["row1_label"], size=12, bold=True,
           color=C_WT, align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE)
    sl.rbox(0.10, LANE2_Y - 0.05, 1.05, 1.75, rgb("FCECEE"), radius=0.2)
    sl.txt(0.10, LANE2_Y - 0.05, 1.05, 1.75, sp["row2_label"], size=12, bold=True,
           color=C_GEN, align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE)

    for i, l in enumerate(layers):
        fg, bg = LAYER_COLOR[l["layer"]]
        accent = rgb(fg)
        # ---- 层带 ----
        sl.rect(xs[i] + 0.05, BAND_Y, l["w"] - 0.10, BAND_H, rgb(bg))
        hw = min(l["w"] - 0.5, 5.6)
        sl.rbox(xs[i] + (l["w"] - hw) / 2, 0.80, hw, 0.40, accent, radius=0.35)
        sl.txt(xs[i] + (l["w"] - hw) / 2, 0.80, hw, 0.40, l["name"], size=12, bold=True,
               color=C_WHITE, align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE)

        cx = xs[i] + (1.30 if i == 0 else 0.28)
        ico = l.get("icon")
        # 文本右边界：为层间 bridge 框与 Gate 菱形让出空间，避免与它们重叠
        gw, gh = 1.30, 1.55
        gx = min(xs[i] + l["w"] - gw - 0.34, W - gw - 1.20) if l.get("gate") else None
        right = xs[i] + l["w"] - (0.98 if i < len(layers) - 1 else 0.28)
        if gx is not None:
            right = min(right, gx - 0.22)
        cw = max(1.4, right - cx)

        # ---- 上泳道：纯文字（不进框） ----
        tx, tw = cx, cw
        if ico:
            sl.pic(ico, cx, LANE1_Y, 0.58, 0.58)
            tx, tw = cx + 0.68, cw - 0.68
        sl.txt(tx, LANE1_Y, tw, 1.45, l.get("row1", ""), size=l.get("row1_size", 10),
               color=C_WT, font=l.get("row1_font", SANS), ls=1.22)

        # ---- 下泳道：纯文字（不进框） ----
        tx, tw = cx, cw
        if ico:
            sl.pic(ico, cx, LANE2_Y, 0.58, 0.58)
            tx, tw = cx + 0.68, cw - 0.68
        sl.txt(tx, LANE2_Y, tw, 1.65, l.get("row2", ""), size=l.get("row2_size", 10),
               color=C_GEN, font=l.get("row2_font", SANS), ls=1.22)

        # ---- r 指标：进框 ----
        yy = CHIP_Y
        for lab, val, hl in l.get("chips", []):
            h = 0.40 + 0.17 * (lab.count("\n") + val.count("\n"))
            sl.chip(cx, yy, min(cw, 4.2), h, lab, val, accent,
                    vcolor=C_RED if hl else None,
                    fill=rgb("FDF1E6") if hl else C_WHITE)
            yy += h + 0.11

        # ---- Gate：菱形 ----
        g = l.get("gate")
        if g:
            gy = (LANE1_Y + LANE2_Y) / 2 - 0.15
            sl.diamond(gx, gy, gw, gh, accent, g["text"], size=8.5)
            sl.tag(gx + gw + 0.22, gy + gh * 0.5, g.get("yes", "YES"), C_GREEN, w=0.54)
            if g.get("els"):
                sl.arrow(gx + gw / 2, gy + gh, gx + gw / 2, gy + gh + 0.40,
                         color=C_RED, w=1.6, dashed=True)
                sl.note_box(gx - 0.30, gy + gh + 0.42, gw + 0.60, g["els"], C_RED, size=8.0)

        # ---- 层间箭头 + bridge 框 ----
        if i < len(layers) - 1:
            ax0 = xs[i] + l["w"] - 0.60
            ax1 = xs[i] + l["w"] + 0.52
            sl.arrow(ax0, LANE1_Y + 0.30, ax1, LANE1_Y + 0.30, color=C_WT, w=2.4)
            sl.arrow(ax0, LANE2_Y + 0.30, ax1, LANE2_Y + 0.30, color=C_GEN, w=2.4)
            if i < len(sp["bridges"]):
                bw = 1.75
                sl.note_box((ax0 + ax1) / 2 - bw / 2, LANE1_Y + 0.60, bw,
                            sp["bridges"][i], rgb("6E6E6E"), size=8.0, fill=C_TOOLBG)

    # ---- 精简页脚：只放读图必需的档位定义 + 一行关键限制 + 出处 ----
    sl.rbox(0.35, FOOT_Y, W - 0.70, FOOT_H, rgb("EFEFEC"),
            line=rgb("ccccc8"), lw=1.0, radius=0.08)
    sl.txt(0.50, FOOT_Y + 0.08, W - 1.0, 0.22, "档位定义（读图必需）", size=10.5,
           bold=True, color=rgb("303030"))
    sl.txt(0.50, FOOT_Y + 0.32, W - 1.0, 0.86, sp["registry"], size=9.5,
           color=rgb("454545"), ls=1.20)
    sl.txt(0.50, FOOT_Y + 1.24, W - 1.0, 0.30, sp["caveat_line"], size=9,
           color=rgb("8A5A00"), ls=1.12)
    sl.txt(0.50, FOOT_Y + 1.58, W - 1.0, 0.20, "数据出处：" + sp["source"], size=7.5,
           color=rgb("606060"), ls=1.05)
    return sl
