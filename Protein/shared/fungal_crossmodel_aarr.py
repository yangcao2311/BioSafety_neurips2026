#!/usr/bin/env python3
"""真菌 panel 的跨模型 AARR 聚合（L3）。

**为什么要做**：真菌 24 蛋白的**推理早就跑完了** —— `results/maskfill*` 四个模型目录下
各有 96 条真菌记录。缺的是聚合：级联文件 `cascade_fungal_real.json` 只有 ESM-2-650M
一个模型，于是"真菌只跑了一个模型"这个说法一直挂着。本脚本从已落盘的推理结果直接聚合，
**零新计算**。

**判据（沿用本层其余各处）**：
* `results/maskfill/` 是**混目录**（ESM-2-650M 与 ProteinDT 混在一起），
  必须按 `model_used` 过滤；**ProteinDT 的文件用 `model` 而非 `model_used`**，两者都要接受。
* 只比较**四个模型都有的掩码档**，否则不同模型的均值来自不同的档组合。
* 真菌与非真菌的对比只用 `random_15pct` —— 那是唯一无锚定、严格可比的档；
  位点档的位点由注释定义，用它比较等于把自变量喂给因变量。
* 取不到的记 `null` 并写原因，不填 0。
"""
from __future__ import annotations

import glob
import json
import statistics as st
from collections import defaultdict
from pathlib import Path

BASE = Path(__file__).resolve().parents[1] / "A-L3.1 Protein Mask & Fill"
PANEL = Path(__file__).resolve().parent / "fungal_virulence_panel.json"
OUT = Path(__file__).resolve().parent / "fungal_crossmodel_aarr.json"

MODEL_DIRS = {
    "ESM-2-650M": ("results/maskfill", "facebook/esm2_t33_650M_UR50D"),
    "ESM-2-150M": ("results/maskfill_esm2_150M", None),
    "ProtBERT": ("results/maskfill_protbert", None),
    "ProteinDT": ("results/maskfill_proteindt", None),
}

# 阴性对照的记录**不与病原体同目录**，且 150M / ProtBERT 还在子目录里 ——
# 用顶层 glob 会得到 n=0，进而把"该模型没有对照"误报成缺数据。
CONTROL_DIRS = {
    "ESM-2-650M": "results/maskfill_negative_control",
    "ESM-2-150M": "results/maskfill_negative_control_extra/ESM2_150M",
    "ProtBERT": "results/maskfill_negative_control_extra/ProtBERT",
    "ProteinDT": "results/maskfill_proteindt",   # 对照与病原体同目录，靠 ctrl_ 前缀区分
}


def model_of(rec: dict) -> str:
    # ProteinDT 的结果文件用 model 而非 model_used
    return str(rec.get("model_used") or rec.get("model") or "")


def load(dirname: str, want_model: str | None, names: set[str]):
    rows = []
    for f in glob.glob(str(BASE / dirname / "*.json")):
        base = Path(f).name
        hit = next((n for n in names if base.startswith(n + "_")), None)
        if hit is None:
            continue
        try:
            r = json.load(open(f))
        except Exception:
            continue
        if want_model and want_model not in model_of(r):
            continue
        aarr = r.get("aarr_total")
        if aarr is None:
            continue
        rows.append({"protein": hit, "strategy": r.get("strategy"),
                     "aarr": float(aarr), "model_field": model_of(r)})
    return rows


def load_group(dirname, want_model, exclude: set[str], prefix_filter=None):
    """按目录取一组记录；exclude 里的蛋白名排除掉。prefix_filter 用于挑对照（ctrl_）。"""
    rows = []
    for f in glob.glob(str(BASE / dirname / "*.json")):
        base = Path(f).name
        if any(base.startswith(n + "_") for n in exclude):
            continue
        if prefix_filter is not None and not base.startswith(prefix_filter):
            continue
        if prefix_filter is None and base.startswith("ctrl_"):
            continue
        try:
            r = json.load(open(f))
        except Exception:
            continue
        if want_model and want_model not in model_of(r):
            continue
        if r.get("aarr_total") is None:
            continue
        rows.append({"strategy": r.get("strategy"), "aarr": float(r["aarr_total"])})
    return rows


def main():
    panel = json.load(open(PANEL))["proteins"]
    names = {p["name"] for p in panel}

    per_model, strat_sets = {}, {}
    for label, (d, want) in MODEL_DIRS.items():
        rows = load(d, want, names)
        per_model[label] = rows
        strat_sets[label] = {r["strategy"] for r in rows}
        print(f"  {label:<12} 真菌记录 {len(rows):>4}  档 {sorted(strat_sets[label])}")

    common = set.intersection(*strat_sets.values()) if strat_sets else set()
    print(f"\n四个模型共有的掩码档：{sorted(common)}")

    out = {
        "analysis": "真菌 panel 跨模型 AARR（从已落盘推理结果聚合，零新计算）",
        "n_proteins": len(names),
        "common_strategies": sorted(common),
        "note": ("results/maskfill 是混目录，按 model_used 过滤；ProteinDT 用 model 字段。"
                 "跨模型比较只用四个模型共有的掩码档。"),
        "by_model": {},
    }
    for label, rows in per_model.items():
        sel = [r for r in rows if r["strategy"] in common]
        by_s = defaultdict(list)
        for r in sel:
            by_s[r["strategy"]].append(r["aarr"])
        out["by_model"][label] = {
            "n_records_common_strata": len(sel),
            "n_proteins": len({r["protein"] for r in sel}),
            "mean_aarr": round(st.mean([r["aarr"] for r in sel]), 4) if sel else None,
            "by_strategy": {s: {"n": len(v), "mean": round(st.mean(v), 4)}
                            for s, v in sorted(by_s.items())},
        }

    print(f"\n{'模型':<12}{'n':>5}{'均值':>9}   逐档均值")
    for label, v in out["by_model"].items():
        per = " · ".join(f"{s} {d['mean']}" for s, d in v["by_strategy"].items())
        print(f"{label:<12}{v['n_records_common_strata']:>5}{v['mean_aarr']:>9}   {per}")

    # ---- 真菌 vs 非真菌病原体 vs 阴性对照（只用 random_15pct）----
    RS = "random_15pct"
    cmp_rows = {}
    print(f"\n{'模型':<12}{'真菌':>9}{'非真菌病原体':>14}{'阴性对照':>10}")
    for label, (d, want) in MODEL_DIRS.items():
        fung = [r["aarr"] for r in per_model[label] if r["strategy"] == RS]
        nonf = [r["aarr"] for r in load_group(d, want, names) if r["strategy"] == RS]
        ctrl = [r["aarr"] for r in load_group(CONTROL_DIRS[label], want, names,
                                              prefix_filter="ctrl_") if r["strategy"] == RS]
        cmp_rows[label] = {
            "fungal": {"n": len(fung), "mean": round(st.mean(fung), 4) if fung else None},
            "non_fungal_pathogen": {"n": len(nonf), "mean": round(st.mean(nonf), 4) if nonf else None},
            "negative_control": {"n": len(ctrl), "mean": round(st.mean(ctrl), 4) if ctrl else None,
                                 "note": None if ctrl else "该模型目录下没有 ctrl_ 前缀记录，如实记 null"},
        }
        f_ = cmp_rows[label]["fungal"]["mean"]
        n_ = cmp_rows[label]["non_fungal_pathogen"]["mean"]
        c_ = cmp_rows[label]["negative_control"]["mean"]
        fmt = lambda x: "n/a" if x is None else f"{x:.4f}"
        print(f"{label:<12}{fmt(f_):>9}{fmt(n_):>14}{fmt(c_):>10}")
    out["random15_comparison"] = {
        "why_only_random_15pct": ("唯一无锚定、严格可比的档。位点档的位点由注释定义，"
                                  "用它比较等于把自变量喂给因变量。"),
        "by_model": cmp_rows,
    }

    OUT.write_text(json.dumps(out, ensure_ascii=False, indent=2))
    print(f"\n-> {OUT}")


if __name__ == "__main__":
    main()
