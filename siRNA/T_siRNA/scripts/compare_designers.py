#!/usr/bin/env python3
"""B-2 三方设计器对比：OligoFormer（专用） vs RNA-FM / RiNALMo（通用 RNA-LM）。

只读两份已产出的级联文件，不重算任何检测器：
  results/rna_fm_designer/cascade_results.json     RNA-FM · RiNALMo 各 30 条
  results/oligoformer_designer/cascade_results.json OligoFormer 30 条

被检验的假设：按沉默效力排序选出的 siRNA，脱靶负担是否更高
（「设计得越好 = 越危险」）。

⚠ 判据：三组走的是**同一个候选池、同一套通道、同一批检测器函数**，
唯一自变量是排序模型。因此组间差异只能来自设计器。
⚠ 不做跨检测器平均：r_1 是计数、r_2 是比例、r_3 可 >1、r_4 是比例，单位不同。
⚠ 样本量小（每组 30 条 / 6 基因），只报中位数、区间与非参数检验，不报效应量结论。
"""
from __future__ import annotations
import json, itertools, statistics as st
from pathlib import Path
from scipy.stats import mannwhitneyu

BASE = Path(__file__).resolve().parent.parent
F_RFD = BASE / "results" / "rna_fm_designer" / "cascade_results.json"
F_OF  = BASE / "results" / "oligoformer_designer" / "cascade_results.json"
OUT   = BASE / "results" / "designer_comparison.json"

DETECTORS = [("r_1", "Channel A 转录组脱靶条数"), ("r_2", "带种子位点基因数占比"),
             ("r_3", "命中基因最高 DepMap 必需性"), ("r_4", "引导链 U 比例")]


def load():
    ents = []
    for f in (F_RFD, F_OF):
        if not f.exists():
            raise SystemExit(f"缺文件：{f}")
        ents += json.load(open(f))["all_results"]
    by = {}
    for e in ents:
        by.setdefault(e["designer"], []).append(e)
    return by


def vals(ents, k):
    return [e[k] for e in ents if isinstance(e[k], (int, float))]


def main():
    by = load()
    labels = list(by)
    print("设计器：", {k: len(v) for k, v in by.items()})

    per = {}
    for lab, ents in by.items():
        d = {"n": len(ents),
             "genes": sorted({e["gene"] for e in ents}),
             "designer_kind": ents[0].get("designer_kind", "通用 RNA 语言模型（MLM NLL 排序）")}
        for k, name in DETECTORS:
            v = vals(ents, k)
            d[k] = {"name": name, "n_measured": len(v),
                    "median": round(st.median(v), 6) if v else None,
                    "min": round(min(v), 6) if v else None,
                    "max": round(max(v), 6) if v else None}
        ess = [e for e in ents if e["L5_details"]["n_common_essential_hit"]]
        d["n_entries_hitting_common_essential"] = len(ess)
        d["most_essential_hits"] = sorted(
            {e["L5_details"]["most_essential_gene"] for e in ents if e["L5_details"]["most_essential_gene"]})
        per[lab] = d

    # 两两非参数检验
    tests = {}
    for a, b in itertools.combinations(labels, 2):
        for k, _ in DETECTORS:
            va, vb = vals(by[a], k), vals(by[b], k)
            if len(va) < 3 or len(vb) < 3:
                continue
            try:
                u, p = mannwhitneyu(va, vb, alternative="two-sided")
            except ValueError:
                continue
            tests[f"{a} vs {b} · {k}"] = {
                "median_a": round(st.median(va), 6), "median_b": round(st.median(vb), 6),
                "U": float(u), "p_mwu": round(float(p), 5),
                "verdict": "有差异" if p < 0.05 else "分不开",
            }

    # 选择重叠：同一候选被不同设计器选中
    sel = {lab: {e["candidate_id"] for e in ents} for lab, ents in by.items()}
    overlap = {f"{a} ∩ {b}": sorted(sel[a] & sel[b]) for a, b in itertools.combinations(labels, 2)}

    out = {
        "analysis": "B-2 siRNA · 三方设计器对比（专用 vs 通用）",
        "design": "同一候选池（每基因 200 条滑窗）、同一脱靶通道、同一批检测器函数；"
                  "唯一自变量 = 用哪个模型排序、取各自 top-5。",
        "hypothesis": "按沉默效力排序（OligoFormer）选出的 siRNA，脱靶负担是否高于按通用 MLM 似然排序（RNA-FM/RiNALMo）。",
        "caveats": [
            "每组 30 条 / 6 基因，样本量小；只报中位数与非参数检验，不出效应量结论。",
            "不做跨检测器平均——四档单位互不相同。",
                "OligoFormer 与两个 RNA-LM 同为判别式 oracle，都不是生成模型。",
        ],
        "per_designer": per,
        "pairwise_tests": tests,
        "selection_overlap": {k: {"n": len(v), "ids": v} for k, v in overlap.items()},
    }
    json.dump(out, open(OUT, "w"), indent=2, ensure_ascii=False)
    print(f"\n写出 {OUT}")

    for lab in labels:
        d = per[lab]
        print(f"\n── {lab}（{d['designer_kind']}）n={d['n']}")
        for k, _ in DETECTORS:
            s = d[k]
            print(f"   {k} 中位 {s['median']}  区间 [{s['min']}, {s['max']}]  测得 {s['n_measured']}/{d['n']}")
        print(f"   命中共同必需基因的记录 {d['n_entries_hitting_common_essential']}/{d['n']}")
    print("\n── 两两检验")
    for k, v in tests.items():
        print(f"   {k}: {v['median_a']} vs {v['median_b']}  p={v['p_mwu']}  → {v['verdict']}")
    print("\n── 选择重叠")
    for k, v in overlap.items():
        print(f"   {k}: {len(v)} 条")


if __name__ == "__main__":
    main()
