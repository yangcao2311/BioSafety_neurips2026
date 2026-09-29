#!/usr/bin/env python3
"""L3 检测器审计补 r_4（结合能）档：阴性对照 vs 病原体面板。

口径与本层其余各档**完全一致**（见 `l3_detector_audit_with_controls.json` 的 method）：

* 病原体侧 = `results/cascade/cascade_results_real_full.json`，只取对照集也有的
  4 个掩码档（active_site / binding_interface / immune_epitope / random_15pct）。
  这与 r_5 / r_6 用的是同一批记录 —— 该文件在 4 档下的 r_5 恰为 673 条，
  与已发布的 `n_pathogen=673` 逐条对上。
* 对照侧 = `results/cascade/cascade_negative_control.json`（27 蛋白 × 4 档 = 108 条）。
* AUROC 0.5 = 分不开，< 0.5 = 方向反，**如实记 REVERSED，不翻符号**。

**r_4 只在登记了共晶配体的蛋白上有值**（r_4 = clamp((−ΔG − 4)/4, 0, 1)，
越高 = 结合越强 = 检测器假定的"越危险"方向）。没有配体的记 None，
不进分母 —— 没有分母就没有数。所以本档的 n 小于其他各档：
病原体 187 / 对照 60（对照 27 蛋白里 15 条登记了本链的共晶配体）。

immune_epitope 档病原体侧 r_4 为 0 条（该档的蛋白都没登记配体），
该档 AUROC 记 null 并写明原因，不拿 0 条去凑一个数。
"""
from __future__ import annotations

import json
import statistics as st
from pathlib import Path

from scipy import stats

TASK = Path(__file__).resolve().parents[1]
L31 = TASK / "A-L3.1 Protein Mask & Fill"
PATHOGEN = L31 / "results" / "cascade" / "cascade_results_real_full.json"
CONTROL = L31 / "results" / "cascade" / "cascade_negative_control.json"
AUDIT = TASK / "shared" / "l3_detector_audit_with_controls.json"
T6 = TASK.parent / "tools" / "t6_detector_audit.json"

STRATA = ["active_site", "binding_interface", "immune_epitope", "random_15pct"]


def load(path: Path) -> list[dict]:
    return json.load(open(path))["all_results"]


def auroc(pathogen: list[float], control: list[float]) -> tuple[float, float]:
    """AUROC = Mann-Whitney U / (n1·n2)，正类 = 病原体。不翻符号。"""
    u = stats.mannwhitneyu(pathogen, control, alternative="two-sided")
    return u.statistic / (len(pathogen) * len(control)), float(u.pvalue)


def verdict(a: float, p: float) -> str:
    if p >= 0.05:
        return "分不开"
    return "病原体更高" if a > 0.5 else "REVERSED · 对照更高"


def block(pathogen: list[float], control: list[float], stratum: str | None = None) -> dict:
    if not pathogen or not control:
        return {
            "auroc": None, "p_mwu": None,
            "n_pathogen": len(pathogen), "n_control": len(control),
            "mean_pathogen": round(st.mean(pathogen), 4) if pathogen else None,
            "mean_control": round(st.mean(control), 4) if control else None,
            "direction": None,
            "note": (f"{stratum} 档在病原体侧有 r_4 的记录为 {len(pathogen)} 条"
                     f"（该档蛋白未登记共晶配体），无法比较，如实记 null。"),
        }
    a, p = auroc(pathogen, control)
    return {
        "auroc": round(a, 4), "p_mwu": float(f"{p:.3g}"),
        "n_pathogen": len(pathogen), "n_control": len(control),
        "mean_pathogen": round(st.mean(pathogen), 4),
        "mean_control": round(st.mean(control), 4),
        "direction": verdict(a, p),
    }


def main():
    pr = [x for x in load(PATHOGEN) if x["strategy"] in STRATA and x.get("r_4") is not None]
    cr = [x for x in load(CONTROL) if x["strategy"] in STRATA and x.get("r_4") is not None]

    by_stratum = {}
    for s in STRATA:
        by_stratum[s] = block([x["r_4"] for x in pr if x["strategy"] == s],
                              [x["r_4"] for x in cr if x["strategy"] == s], s)
    overall = block([x["r_4"] for x in pr], [x["r_4"] for x in cr])

    entry = {
        "detector": "r_4 结合能（AutoDock Vina，filled_seq → ESMFold → 对接）",
        "value_definition": "r_4 = clamp((−ΔG − 4) / 4, 0, 1)，越高 = 结合越强",
        "pathogen_source": "cascade_results_real_full.json，对照集也有的 4 个掩码档",
        "control_source": "cascade_negative_control.json（27 蛋白 × 4 档）",
        "coverage": {
            "n_pathogen_with_r_4": len(pr),
            "n_control_with_r_4": len(cr),
            "n_control_proteins_with_ligand": len({x["protein_name"] for x in cr}),
            "note": ("r_4 只在登记了共晶配体的蛋白上有值，其余记 None 且不进分母。"
                     "对照 27 蛋白中 15 条登记了本链共晶配体（ctrl_TXN 丢弃：3QFA 的 FAD "
                     "属于硫氧还蛋白还原酶那条链；ctrl_HBA1 / ctrl_HBB 按链各取本链血红素；"
                     "ctrl_TUBB 在 8T42 链 B 找不到 GTP，未登记）。"),
        },
        "by_strategy": by_stratum,
        "overall": overall,
    }

    audit = json.load(open(AUDIT))
    audit["full_cascade"]["r_4"] = entry
    json.dump(audit, open(AUDIT, "w"), ensure_ascii=False, indent=1)

    t6 = json.load(open(T6))
    t6["layers"]["L3 蛋白"]["rungs"]["r_4"] = entry
    json.dump(t6, open(T6, "w"), ensure_ascii=False, indent=1)

    print(json.dumps(entry, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
