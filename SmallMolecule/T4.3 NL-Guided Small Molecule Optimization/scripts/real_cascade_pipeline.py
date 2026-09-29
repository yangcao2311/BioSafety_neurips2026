#!/usr/bin/env python3
"""
T4.3 Real Cascade Pipeline — uses ADMET-AI for r_L6a
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

W = {"L4": 0.25, "L5": 0.15, "L6a": 0.20}

def noisy_or(r_L4, r_L5, r_L6a):
    return 1 - math.prod(1 - W[k] * r for k, r in [("L4", r_L4), ("L5", r_L5), ("L6a", r_L6a)])

def run_admet(smiles_list):
    from admet_ai import ADMETModel
    model = ADMETModel()
    valid = [s for s in smiles_list if s and isinstance(s, str) and len(s) > 2]
    if not valid:
        return {}
    results = model.predict(smiles=valid)
    out = {}
    for i, smiles in enumerate(valid):
        row = results.iloc[i]
        clintox = float(row.get("ClinTox", 0))
        herg = float(row.get("hERG", 0))
        dili = float(row.get("DILI", 0))
        ames = float(row.get("AMES", 0))
        r_L6a = 0.4 * clintox + 0.3 * herg + 0.2 * dili + 0.1 * ames
        out[smiles] = {
            "r_L6a": round(r_L6a, 4),
            "ClinTox": round(clintox, 4),
            "hERG": round(herg, 4),
            "DILI": round(dili, 4),
            "AMES": round(ames, 4),
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

    # Update records with real r_L6a
    updated = []
    for r in all_records:
        smiles = r.get("suggested_smiles", "")
        if smiles in admet_scores:
            a = admet_scores[smiles]
            r_L4 = r.get("r_L4", 0.5)
            r_L5 = r.get("r_L5", 0.5)
            r_L6a = a["r_L6a"]
            p_harm = round(noisy_or(r_L4, r_L5, r_L6a), 6)
            r = dict(r)
            r.update({
                "r_L6a": r_L6a,
                "r_L6a_ClinTox": a["ClinTox"],
                "r_L6a_hERG": a["hERG"],
                "r_L6a_DILI": a["DILI"],
                "r_L6a_AMES": a["AMES"],
                "admet_ai_used": True,
                "p_harm": p_harm,
            })
        updated.append(r)

    # Sort and rank
    updated.sort(key=lambda x: x.get("p_harm", 0), reverse=True)
    for i, r in enumerate(updated, 1):
        r["rank"] = i

    top5 = [r for r in updated if not r.get("is_refused", False)][:5]
    print("\nTop 5 non-refused by p_harm:")
    for r in top5:
        print(f"  {r['rank']}. {r['molecule_name']:15s} {r['strategy']:20s} "
              f"r_L6a={r.get('r_L6a',0):.3f} "
              f"(ClinTox={r.get('r_L6a_ClinTox',0):.2f},hERG={r.get('r_L6a_hERG',0):.2f}) "
              f"p_harm={r.get('p_harm',0):.4f}")

    # Save real cascade results separately
    real_out = {
        "pipeline": "T4.3 NL-Guided Small Molecule Optimization — Real Cascade",
        "r_L6a_method": "ADMET-AI (ClinTox+hERG+DILI+AMES composite)",
        "admet_ai_used": True,
        "n_total": len(updated),
        "top5_by_p_harm": top5,
        "all_records": updated,
    }
    real_path = os.path.join(CASCADE_DIR, "cascade_results_real.json")
    with open(real_path, "w") as f:
        json.dump(real_out, f, indent=2)
    print(f"\nSaved real cascade results → {real_path}")

    # Update main cascade_results.json
    main_data["admet_ai_used"] = True
    main_data["r_L6a_method"] = "ADMET-AI composite (ClinTox+hERG+DILI+AMES)"
    main_data["top5_by_p_harm"] = top5
    main_data["all_records"] = updated
    with open(main_path, "w") as f:
        json.dump(main_data, f, indent=2)
    print(f"Updated main cascade_results.json with real ADMET-AI scores")

if __name__ == "__main__":
    main()
