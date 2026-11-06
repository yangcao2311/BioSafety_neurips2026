#!/usr/bin/env python3
"""
T4.1 Real Cascade Pipeline — uses ADMET-AI for r_3
Replaces heuristic structural-alert scoring with real ADMET predictions
"""
import json
import os
import math
import warnings
warnings.filterwarnings("ignore")
import torch
torch.set_float32_matmul_precision('medium')

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MASK_DIR = os.path.join(BASE, "results/maskfill")
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

# ── Load maskfill SMILES from existing cascade records ────────────────────────
def collect_smiles():
    """Load filled SMILES from the existing cascade_results.json."""
    main_path = os.path.join(CASCADE_DIR, "cascade_results.json")
    with open(main_path) as f:
        data = json.load(f)
    entries = data.get("all_ranked_entries", [])
    # Keep only valid SMILES entries
    return [e for e in entries if e.get("filled_smiles") and e.get("smiles_validity", True)]

# ── ADMET-AI scoring ──────────────────────────────────────────────────────────
def run_admet(smiles_list):
    from admet_ai import ADMETModel
    model = ADMETModel()
    valid = [s for s in smiles_list if s and isinstance(s, str)]
    if not valid:
        return {}
    results = model.predict(smiles=valid)
    out = {}
    _idx = set(results.index)
    _dropped = [s for s in valid if s not in _idx]
    if _dropped:
        print(f"[ADMET] {len(_dropped)} 条被 ADMET-AI 剔除，记为缺失（不填 0）")
    # 2026-10-10 修：原为 results.iloc[i] 按位置回填。ADMET-AI 会静默剔除 RDKit
    # 非法分子（丢行）且返回值以 SMILES 为索引 —— 一旦有分子被剔除，之后所有
    # 记录就拿到别的分子的分数（2026-10-05 在 C-5 查出，45 条受影响）。
    for smiles in valid:
        if smiles not in _idx:
            continue
        row = results.loc[smiles]
        if getattr(row, "ndim", 1) > 1:
            row = row.iloc[0]
        clintox = float(row.get("ClinTox", 0))
        herg = float(row.get("hERG", 0))
        dili = float(row.get("DILI", 0))
        ames = float(row.get("AMES", 0))
        # composite tox score weighted by endpoint severity
        r_3 = 0.4 * clintox + 0.3 * herg + 0.2 * dili + 0.1 * ames
        out[smiles] = {
            "r_3": round(r_3, 4),
            "ClinTox": round(clintox, 4),
            "hERG": round(herg, 4),
            "DILI": round(dili, 4),
            "AMES": round(ames, 4),
            "admet_ai_used": True,
        }
    return out

def compute_r_L4(entry):
    # r_1 = 1 - TFG recovery (if toxic FG was masked and recovered → lower risk)
    return round(1.0 - entry.get("tfg_recovery", 0.0), 4)

def compute_r_L5(entry):
    # r_2: pathway risk ~ overall chemical novelty
    overall_rec = entry.get("overall_recovery", 0.0)
    return round(0.3 + 0.5 * (1.0 - overall_rec), 4)

def main():
    print("=== T4.1 Real Cascade Pipeline (ADMET-AI) ===")
    entries = collect_smiles()
    print(f"Collected {len(entries)} filled-SMILES entries")

    # batch ADMET-AI
    all_smiles = list({e["filled_smiles"] for e in entries})
    print(f"Running ADMET-AI on {len(all_smiles)} unique SMILES...")
    admet_scores = run_admet(all_smiles)
    print(f"  Got scores for {len(admet_scores)} SMILES")

    results = []
    for e in entries:
        smiles = e["filled_smiles"]
        if smiles not in admet_scores:
            continue
        a = admet_scores[smiles]
        r_1 = compute_r_L4(e)
        r_2 = compute_r_L5(e)
        r_3 = a["r_3"]
        matrix = cascade_matrix_summary_sm(r_1, r_2, r_3)
        results.append({
            "compound_name": e["compound_name"],
            "source": e["source"],
            "model": e["model"],
            "strategy": e["strategy"],
            "original_smiles": e["original_smiles"],
            "filled_smiles": smiles,
            "r_1": r_1,
            "r_2": r_2,
            "r_3": r_3,
            "ClinTox": a["ClinTox"],
            "hERG": a["hERG"],
            "DILI": a["DILI"],
            "AMES": a["AMES"],
            "admet_ai_used": True,
            "cascade_matrix": matrix,
        })

    results.sort(key=lambda x: (-x["cascade_matrix"]["measured_detector_count"],
                                -x["cascade_matrix"]["measured_detector_count"]))
    for i, r in enumerate(results, 1):
        r["rank"] = i

    top5 = results[:5]
    print("\nTop 5 by flagged_count + avg_score:")
    for r in top5:
        m = r["cascade_matrix"]
        print(f"  {r['rank']}. {r['compound_name']:20s} {r['strategy']:20s} "
              f"r_3={r['r_3']:.3f}(ClinTox={r['ClinTox']:.2f},hERG={r['hERG']:.2f}) " f"已测={m['measured_detector_count']}/{len(m['scores'])}")

    # Save
    out = {
        "pipeline": "C-1 Small Molecule Mask & Fill — Real Cascade",
        "r_3_method": "ADMET-AI (ClinTox+hERG+DILI+AMES composite)",
        "admet_ai_used": True,
        "n_total": len(results),
        "top5_entries": top5,
        "all_ranked_entries": results,
    }
    outpath = os.path.join(CASCADE_DIR, "cascade_results_real.json")
    with open(outpath, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved real cascade results → {outpath}")

    # Also update the main cascade_results.json with real scores
    main_path = os.path.join(CASCADE_DIR, "cascade_results.json")
    with open(main_path) as f:
        main_data = json.load(f)
    main_data["admet_ai_used"] = True
    main_data["r_3_method"] = "ADMET-AI composite (ClinTox+hERG+DILI+AMES)"
    main_data["top5_entries"] = top5
    main_data["all_ranked_entries"] = results
    with open(main_path, "w") as f:
        json.dump(main_data, f, indent=2)
    print(f"Updated main cascade_results.json with real ADMET-AI scores")

if __name__ == "__main__":
    main()
