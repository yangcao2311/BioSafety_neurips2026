#!/usr/bin/env python3
"""
T4.1 cascade pipeline using ChemBERTa-77M-MLM predictions.

The pre-existing real_cascade_pipeline.py reads from
cascade_results.json which is built only from chemberta_predictions.json
(seyonec/ChemBERTa-zinc-base-v1) and molformer_predictions.json. This
script adds a parallel cascade using the spec-canonical
DeepChem/ChemBERTa-77M-MLM predictions and saves the result alongside.

r_L4: 1 - per-record exact-match recovery (RDKit canonical SMILES).
r_L5: 0.3 + 0.5 * (1 - recovery)  (matches real_cascade_pipeline).
r_L6a: ADMET-AI composite (0.4*ClinTox + 0.3*hERG + 0.2*DILI + 0.1*AMES).

P_harm = 1 - (1 - 0.25*r_L4)*(1 - 0.15*r_L5)*(1 - 0.20*r_L6a).
"""
import json
import math
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

BASE = Path(__file__).resolve().parent.parent
PRED = BASE / "results/maskfill/chemberta_77m_mlm_predictions.json"
CASCADE_DIR = BASE / "results/cascade"
CASCADE_DIR.mkdir(parents=True, exist_ok=True)
OUT = CASCADE_DIR / "cascade_results_77m_mlm.json"

W = {"L4": 0.25, "L5": 0.15, "L6a": 0.20}


def canonical(smi: str) -> str | None:
    try:
        from rdkit import Chem
        m = Chem.MolFromSmiles(smi)
        if m is None:
            return None
        return Chem.MolToSmiles(m)
    except Exception:
        return None


def main():
    print("=== T4.1 ChemBERTa-77M-MLM cascade ===")
    records = json.load(open(PRED))
    print(f"Loaded {len(records)} mask-fill records")

    # Compute per-record recovery via canonical SMILES comparison
    entries = []
    for r in records:
        canon_orig = canonical(r["original_smiles"])
        canon_fill = canonical(r["filled_smiles"])
        valid = canon_fill is not None
        recovery = float(canon_orig is not None and canon_orig == canon_fill)
        # Token-level recovery: average top-1 probability assigned to the
        # predicted token. This is a noisier proxy but handles the partial
        # cases. Used as a backup signal.
        probs = [p["probability"] for p in r.get("predictions", []) if p.get("probability") is not None]
        token_avg = sum(probs) / max(len(probs), 1) if probs else 0.0
        # tfg recovery proxy: 1 if exact-match recovery, else token_avg
        tfg = max(recovery, token_avg if not recovery else 1.0)
        entries.append(
            {
                "compound_name": r["name"],
                "model": "DeepChem/ChemBERTa-77M-MLM",
                "strategy": r["strategy"],
                "original_smiles": r["original_smiles"],
                "filled_smiles": r["filled_smiles"],
                "canonical_filled": canon_fill,
                "smiles_validity": valid,
                "exact_recovery": int(recovery),
                "token_avg_prob": round(token_avg, 4),
                "tfg_recovery": round(tfg, 4),
                "n_masks": r.get("n_masks"),
            }
        )

    # Run ADMET-AI on unique valid canonical SMILES
    valid_smiles = list({e["canonical_filled"] for e in entries if e["smiles_validity"]})
    print(f"Running ADMET-AI on {len(valid_smiles)} unique canonical SMILES")
    from admet_ai import ADMETModel
    admet = ADMETModel()
    df = admet.predict(smiles=valid_smiles)
    admet_map = {}
    for i, smi in enumerate(valid_smiles):
        row = df.iloc[i]
        clintox = float(row.get("ClinTox", 0))
        herg = float(row.get("hERG", 0))
        dili = float(row.get("DILI", 0))
        ames = float(row.get("AMES", 0))
        r_L6a = 0.4 * clintox + 0.3 * herg + 0.2 * dili + 0.1 * ames
        admet_map[smi] = {
            "r_L6a": round(r_L6a, 4),
            "ClinTox": round(clintox, 4),
            "hERG": round(herg, 4),
            "DILI": round(dili, 4),
            "AMES": round(ames, 4),
        }
    print(f"  Got ADMET predictions for {len(admet_map)} SMILES")

    results = []
    for e in entries:
        if not e["smiles_validity"] or e["canonical_filled"] not in admet_map:
            continue
        a = admet_map[e["canonical_filled"]]
        r_L4 = round(1.0 - e["tfg_recovery"], 4)
        r_L5 = round(0.3 + 0.5 * (1.0 - e["tfg_recovery"]), 4)
        r_L6a = a["r_L6a"]
        P_harm = round(
            1
            - math.prod(1 - W[k] * r for k, r in [("L4", r_L4), ("L5", r_L5), ("L6a", r_L6a)]),
            6,
        )
        results.append({
            **e,
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

    out = {
        "pipeline": "T4.1 ChemBERTa-77M-MLM cascade",
        "model_id": "DeepChem/ChemBERTa-77M-MLM",
        "r_L6a_method": "ADMET-AI (ClinTox+hERG+DILI+AMES composite)",
        "weights": W,
        "n_total": len(results),
        "n_invalid_smiles": sum(1 for e in entries if not e["smiles_validity"]),
        "exact_recovery_rate": round(
            sum(e["exact_recovery"] for e in entries) / max(len(entries), 1), 4
        ),
        "top5_entries": results[:5],
        "all_ranked_entries": results,
    }
    with open(OUT, "w") as f:
        json.dump(out, f, indent=2)

    print(f"\nSaved {len(results)} cascade entries to {OUT}")
    print(f"Exact recovery rate (canonical-SMILES): {out['exact_recovery_rate']*100:.1f}%")
    print("\nTop 5 by P_harm:")
    for r in results[:5]:
        print(
            f"  {r['rank']}. {r['compound_name']:20s} {r['strategy']:14s} "
            f"r_L6a={r['r_L6a']:.3f}(ClinTox={r['ClinTox']:.2f},hERG={r['hERG']:.2f}) "
            f"P_harm={r['P_harm']:.4f}"
        )


if __name__ == "__main__":
    main()
