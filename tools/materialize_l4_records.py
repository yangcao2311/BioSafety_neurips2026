#!/usr/bin/env python3
"""把两组只有聚合形式的 L4 结果物化成逐记录级联格式，供阈值敏感性表统一读取。

不新增任何测量 —— 只把已落盘的数按**文件自述的口径**换算成 r 档：

1. **RFdiffusion + ProteinMPNN → 官方 AF3**（阴性对照轮，42 条设计）
   源：`results/l4_iptm_control_analysis.json` 的 `by_target`（每靶点 3 条设计的 ipTM）
   口径取自该文件自己的 `metric` 字段：`r = clip(0.2 + 0.8 × ipTM, 0.2, 1.0)`
   → 落到 L4 的 `r_2`（界面置信度）。r_1 / r_3 本轮未算，记 None。

2. **C11 复合物界面掩码补全**（ESM-2 650M / 150M，各 11 靶点）
   源：`A-L4.1/results/interface_maskfill_analysis.json` 的 `per_target`
   → `r_1` = 界面档 AARR，`r_1_control` = 同链非界面对照档 AARR
   ⚠ 这里一条"记录"是一个靶点，不是一条设计。
"""
import json
from pathlib import Path

BASE = Path("/home/ycao95/BioSafety/Task/A_L4_Complex/A-L4.1 PPI Binder Design")
OUT = BASE / "results" / "cascade"


def clip(x, lo, hi): return max(lo, min(hi, x))


def main():
    # ---- 1. 阴性对照轮 ----
    src = json.loads((BASE / "results/l4_iptm_control_analysis.json").read_text())
    recs = []
    for side, targets in src["by_target"].items():
        for tgt, iptms in targets.items():
            for i, v in enumerate(iptms, 1):
                recs.append({
                    "target": tgt, "cohort": side, "design_index": i,
                    "chain_pair_iptm": v,
                    "r_1": None, "r_1_reason": "本轮未算 seq_recovery",
                    "r_2": round(clip(0.2 + 0.8 * v, 0.2, 1.0), 4),
                    "r_2_method": "clip(0.2 + 0.8 × chain_pair_iptm, 0.2, 1.0) —— 取自源文件 metric 字段",
                    "r_3": None, "r_3_reason": "本轮未算 MHCflurry",
                    "model": "RFdiffusion(CPU) + ProteinMPNN v_48_020 → 官方 AlphaFold 3",
                })
    p1 = OUT / "cascade_l4_negative_control_designs.json"
    p1.write_text(json.dumps({
        "task": "A-L4.1 · binder 阴性对照轮（逐设计物化）",
        "model": "RFdiffusion + ProteinMPNN → 官方 AlphaFold 3",
        "source": "results/l4_iptm_control_analysis.json :: by_target",
        "note": "两侧逐字相同流程；靶点级为主分析，本文件是设计级物化。"
                "⚠ 同靶点 3 条设计共享骨架，设计级属伪重复，仅供逐记录统计用。",
        "n_total": len(recs), "all_results": recs,
    }, ensure_ascii=False, indent=2))
    print(f"{p1.name}: {len(recs)} 条")

    # ---- 2. C11 界面掩码补全 ----
    a = json.loads((BASE / "results/interface_maskfill_analysis.json").read_text())
    for tag, m in a["by_model"].items():
        rows = []
        for t in m["per_target"]:
            rows.append({
                "target": t["target"], "cohort": t["cohort"], "n_masked": t["n_masked"],
                "r_1": t["aarr_interface"],
                "r_1_method": "界面残基 AARR（ESMFold 掩码补全，界面 = ≤4.5 Å 于伙伴链）",
                "r_1_control_random_sites": t["aarr_random"],
                "r_2": None, "r_2_reason": "本 task 无 ipTM 档",
                "r_3": None, "r_3_reason": "本 task 未算 MHCflurry",
                "model": tag,
            })
        p = OUT / f"cascade_c11_interface_{tag}.json"
        p.write_text(json.dumps({
            "task": "C11 · A-L4 G1 复合物界面掩码补全（逐靶点）",
            "model": tag, "source": "results/interface_maskfill_analysis.json :: per_target",
            "note": "一条记录 = 一个靶点，不是一条设计。r_1 是界面档 AARR；"
                    "同链非界面对照档在 r_1_control_random_sites。",
            "n_total": len(rows), "all_results": rows,
        }, ensure_ascii=False, indent=2))
        print(f"{p.name}: {len(rows)} 条")


if __name__ == "__main__":
    main()
