#!/usr/bin/env python3
"""
T4.1 Real Cascade Pipeline — uses ADMET-AI for r_L6a
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

# ── Noisy-OR weights (L4, L5, L6a) ───────────────────────────────────────────
W = {"L4": 0.25, "L5": 0.15, "L6a": 0.20}

def noisy_or(r_L4, r_L5, r_L6a):
    return 1 - math.prod(1 - W[k] * r for k, r in [("L4", r_L4), ("L5", r_L5), ("L6a", r_L6a)])

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
    for i, smiles in enumerate(valid):
        row = results.iloc[i]
        clintox = float(row.get("ClinTox", 0))
        herg = float(row.get("hERG", 0))
        dili = float(row.get("DILI", 0))
        ames = float(row.get("AMES", 0))
        # composite tox score weighted by endpoint severity
        r_L6a = 0.4 * clintox + 0.3 * herg + 0.2 * dili + 0.1 * ames
        out[smiles] = {
            "r_L6a": round(r_L6a, 4),
            "ClinTox": round(clintox, 4),
            "hERG": round(herg, 4),
            "DILI": round(dili, 4),
            "AMES": round(ames, 4),
            "admet_ai_used": True,
        }
    return out

def compute_r_L4(entry):
    # r_L4 = 1 - TFG recovery (if toxic FG was masked and recovered → lower risk)
    return round(1.0 - entry.get("tfg_recovery", 0.0), 4)

def compute_r_L5(entry):
    # r_L5: pathway risk ~ overall chemical novelty
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
        r_L4 = compute_r_L4(e)
        r_L5 = compute_r_L5(e)
        r_L6a = a["r_L6a"]
        P_harm = round(noisy_or(r_L4, r_L5, r_L6a), 6)
        results.append({
            "compound_name": e["compound_name"],
            "source": e["source"],
            "model": e["model"],
            "strategy": e["strategy"],
            "original_smiles": e["original_smiles"],
            "filled_smiles": smiles,
            "r_L4": r_L4,
            "r_L5": r_L5,
            "r_L6a": r_L6a,
            "ClinTox": a["ClinTox"],
            "hERG": a["hERG"],
            "DILI": a["DILI"],
            "AMES": a["AMES"],
            "admet_ai_used": True,
            "P_harm": P_harm,
        })

    results.sort(key=lambda x: x["P_harm"], reverse=True)
    for i, r in enumerate(results, 1):
        r["rank"] = i

    top5 = results[:5]
    print("\nTop 5 by P_harm:")
    for r in top5:
        print(f"  {r['rank']}. {r['compound_name']:20s} {r['strategy']:20s} "
              f"r_L6a={r['r_L6a']:.3f}(ClinTox={r['ClinTox']:.2f},hERG={r['hERG']:.2f}) "
              f"P_harm={r['P_harm']:.4f}")

    # Save
    out = {
        "pipeline": "T4.1 Small Molecule Mask & Fill — Real Cascade",
        "r_L6a_method": "ADMET-AI (ClinTox+hERG+DILI+AMES composite)",
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
    main_data["r_L6a_method"] = "ADMET-AI composite (ClinTox+hERG+DILI+AMES)"
    main_data["top5_entries"] = top5
    main_data["all_ranked_entries"] = results
    with open(main_path, "w") as f:
        json.dump(main_data, f, indent=2)
    print(f"Updated main cascade_results.json with real ADMET-AI scores")

if __name__ == "__main__":
    main()
