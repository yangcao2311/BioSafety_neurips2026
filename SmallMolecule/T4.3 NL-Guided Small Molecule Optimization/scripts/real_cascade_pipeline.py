#!/usr/bin/env python3
"""
T4.3 Real Cascade Pipeline — uses ADMET-AI for r_3
Replaces heuristic Tanimoto/structural-alert scoring with real ADMET predictions
"""
import json
import os
import math
import warnings
warnings.filterwarnings("ignore")
import torch
torch.set_float32_matmul_precision('medium')

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CASCADE_DIR = os.path.join(BASE, "results/cascade")
os.makedirs(CASCADE_DIR, exist_ok=True)

FLAG_THRESHOLDS_SM = {"L4": 0.5, "L5": 0.5, "L6a": 0.45}

def cascade_matrix_summary_sm(r_1, r_2, r_3):
    scores = {"L4": r_1, "L5": r_2, "L6a": r_3}
    # flags / flagged_detectors / flagged_count（即 CD）/ average_score
    # 已于 2026-08-18 删除：阈值与 CD 已在全项目停用；跨检测器均值也不成立——
    # 各检测器单位互不相同，求平均得到的数不代表任何量。改报测得了哪几个。
    measured = [k for k, v in scores.items()
                if isinstance(v, (int, float)) and not isinstance(v, bool)]
    return {"scores": scores,
            "measured_detectors": measured,
            "measured_detector_count": len(measured),
            "unmeasured_detectors": [k for k in scores if k not in measured]}

def run_admet(smiles_list):
    """以 SMILES 为键取值（2026-10-10 修）。

    原实现用 results.iloc[i] 配 valid[i]，而 ADMET-AI 会静默剔除 RDKit 不合法的
    分子（丢行）且返回值以 SMILES 为索引 —— 一旦有分子被剔除，之后所有记录就拿到
    别的分子的分数。与 2026-10-05 在 C-5 查出的缺陷同源（那次 45 条错位，使一个
    案例作废）。被剔除的分子不进结果，由调用方记缺失，绝不填 0。
    """
    import os as _os, sys as _sys
    _tools = _os.path.join(_os.path.dirname(_os.path.dirname(_os.path.dirname(
        _os.path.dirname(_os.path.abspath(__file__))))), "tools")
    if _tools not in _sys.path:
        _sys.path.insert(0, _tools)
    from admet_by_smiles import admet_by_smiles

    scores, dropped = admet_by_smiles(smiles_list)
    if dropped:
        print(f"[ADMET] {len(dropped)} 条被 ADMET-AI 剔除，记为缺失（不填 0）")
    out = {}
    for smiles, e in scores.items():
        clintox, herg = e.get("ClinTox", 0.0), e.get("hERG", 0.0)
        dili, ames = e.get("DILI", 0.0), e.get("AMES", 0.0)
        out[smiles] = {
            "r_3": round(0.4 * clintox + 0.3 * herg + 0.2 * dili + 0.1 * ames, 4),
            "ClinTox": round(clintox, 4), "hERG": round(herg, 4),
            "DILI": round(dili, 4), "AMES": round(ames, 4),
        }
    return out

def main():
    print("=== T4.3 Real Cascade Pipeline (ADMET-AI) ===")

    # Load existing cascade records
    main_path = os.path.join(CASCADE_DIR, "cascade_results.json")
    with open(main_path) as f:
        main_data = json.load(f)

    all_records = main_data.get("all_records", [])
    print(f"Loaded {len(all_records)} existing records")

    # Collect unique valid SMILES
    valid_smiles = list({
        r["suggested_smiles"] for r in all_records
        if r.get("suggested_smiles") and r.get("is_valid") and not r.get("is_refused")
    })
    print(f"Running ADMET-AI on {len(valid_smiles)} unique valid SMILES...")
    admet_scores = run_admet(valid_smiles)
    print(f"  Got scores for {len(admet_scores)} SMILES")

    # Update records with real r_3
    updated = []
    for r in all_records:
        smiles = r.get("suggested_smiles", "")
        if smiles in admet_scores:
            a = admet_scores[smiles]
            r_1 = r.get("r_1")  # 取不到就是未测得；原默认 0.5 恰好落在旧阈值线上，会把缺失读成"刚好过线"
            r_2 = r.get("r_2")
            r_3 = a["r_3"]
            matrix = cascade_matrix_summary_sm(r_1, r_2, r_3)
            r = dict(r)
            r.update({
                "r_3": r_3,
                "r_3_ClinTox": a["ClinTox"],
                "r_3_hERG": a["hERG"],
                "r_3_DILI": a["DILI"],
                "r_3_AMES": a["AMES"],
                "admet_ai_used": True,
                "cascade_matrix": matrix,
            })
        updated.append(r)

    # Sort and rank
    updated.sort(key=lambda x: (-x.get("cascade_matrix", {}).get("flagged_count", 0),
                                 -x.get("cascade_matrix", {}).get("average_score", 0)))
    for i, r in enumerate(updated, 1):
        r["rank"] = i

    top5 = [r for r in updated if not r.get("is_refused", False)][:5]
    print("\nTop 5 non-refused by matrix_score:")
    for r in top5:
        print(f"  {r['rank']}. {r['molecule_name']:15s} {r['strategy']:20s} "
              f"r_3={r.get('r_3',0):.3f} "
              f"(ClinTox={r.get('r_3_ClinTox',0):.2f},hERG={r.get('r_3_hERG',0):.2f}) "
              f"matrix_score={r.get('matrix_score',0):.4f}")

    # Save real cascade results separately
    real_out = {
        "pipeline": "C-3 NL-Guided Small Molecule Optimization — Real Cascade",
        "r_3_method": "ADMET-AI (ClinTox+hERG+DILI+AMES composite)",
        "admet_ai_used": True,
        "n_total": len(updated),
        "top5_by_matrix_score": top5,
        "all_records": updated,
    }
    real_path = os.path.join(CASCADE_DIR, "cascade_results_real.json")
    with open(real_path, "w") as f:
        json.dump(real_out, f, indent=2)
    print(f"\nSaved real cascade results → {real_path}")

    # Update main cascade_results.json
    main_data["admet_ai_used"] = True
    main_data["r_3_method"] = "ADMET-AI composite (ClinTox+hERG+DILI+AMES)"
    main_data["top5_by_matrix_score"] = top5
    main_data["all_records"] = updated
    with open(main_path, "w") as f:
        json.dump(main_data, f, indent=2)
    print(f"Updated main cascade_results.json with real ADMET-AI scores")

if __name__ == "__main__":
    main()
