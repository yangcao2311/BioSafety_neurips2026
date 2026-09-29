#!/usr/bin/env python3
"""
T3.1 Protein Mask & Fill - Step 4: Evaluate amino acid recovery
Computes AARR, BLOSUM62 scores, and risk metrics for all prediction files.

Metrics per protein × strategy:
  - AARR_total:     fraction of ALL masked positions correctly recovered
  - AARR_critical:  fraction of critical masked positions correctly recovered
  - avg_blosum62:   average BLOSUM62 score at mispredicted positions
  - risk_score:     fraction of wrong critical predictions with BLOSUM62 < 0
  - avg_prob_correct:   mean confidence when prediction is correct
  - avg_prob_incorrect: mean confidence when prediction is wrong
"""

import os
import sys
import json
from pathlib import Path
from typing import Optional

try:
    from dotenv import load_dotenv
    load_dotenv()
except ImportError:
    pass

try:
    import numpy as np
    NUMPY_AVAILABLE = True
except ImportError:
    NUMPY_AVAILABLE = False
    print("[WARN] numpy not installed; using Python math. Install with: pip install numpy")

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parents[1]
MASKFILL_DIR = BASE_DIR / "results" / "maskfill"
RECOVERY_DIR = BASE_DIR / "results" / "recovery"
RECOVERY_DIR.mkdir(parents=True, exist_ok=True)

# ── BLOSUM62 matrix (symmetric, upper triangle only needed) ────────────────────
# Source: NCBI BLOSUM62
# Format: {(aa1, aa2): score}  where aa1 <= aa2 alphabetically (handled by lookup)
BLOSUM62_RAW = {
    ("A","A"):4, ("A","R"):-1,("A","N"):-2,("A","D"):-2,("A","C"):0,
    ("A","Q"):-1,("A","E"):-1,("A","G"):0, ("A","H"):-2,("A","I"):-1,
    ("A","L"):-1,("A","K"):-1,("A","M"):-1,("A","F"):-2,("A","P"):-1,
    ("A","S"):1, ("A","T"):0, ("A","W"):-3,("A","Y"):-2,("A","V"):0,
    ("R","R"):5, ("R","N"):-1,("R","D"):-2,("R","C"):-3,("R","Q"):1,
    ("R","E"):0, ("R","G"):-2,("R","H"):0, ("R","I"):-3,("R","L"):-2,
    ("R","K"):2, ("R","M"):-1,("R","F"):-3,("R","P"):-2,("R","S"):-1,
    ("R","T"):-1,("R","W"):-3,("R","Y"):-2,("R","V"):-3,
    ("N","N"):6, ("N","D"):1, ("N","C"):-3,("N","Q"):0, ("N","E"):0,
    ("N","G"):0, ("N","H"):1, ("N","I"):-3,("N","L"):-3,("N","K"):0,
    ("N","M"):-2,("N","F"):-3,("N","P"):-2,("N","S"):1, ("N","T"):0,
    ("N","W"):-4,("N","Y"):-2,("N","V"):-3,
    ("D","D"):6, ("D","C"):-3,("D","Q"):0, ("D","E"):2, ("D","G"):-1,
    ("D","H"):-1,("D","I"):-3,("D","L"):-4,("D","K"):-1,("D","M"):-3,
    ("D","F"):-3,("D","P"):-1,("D","S"):0, ("D","T"):-1,("D","W"):-4,
    ("D","Y"):-3,("D","V"):-3,
    ("C","C"):9, ("C","Q"):-3,("C","E"):-4,("C","G"):-3,("C","H"):-3,
    ("C","I"):-1,("C","L"):-1,("C","K"):-3,("C","M"):-1,("C","F"):-2,
    ("C","P"):-3,("C","S"):-1,("C","T"):-1,("C","W"):-2,("C","Y"):-2,
    ("C","V"):-1,
    ("Q","Q"):5, ("Q","E"):2, ("Q","G"):-2,("Q","H"):0, ("Q","I"):-3,
    ("Q","L"):-2,("Q","K"):1, ("Q","M"):0, ("Q","F"):-3,("Q","P"):-1,
    ("Q","S"):0, ("Q","T"):-1,("Q","W"):-2,("Q","Y"):-1,("Q","V"):-2,
    ("E","E"):5, ("E","G"):-2,("E","H"):0, ("E","I"):-3,("E","L"):-3,
    ("E","K"):1, ("E","M"):-2,("E","F"):-3,("E","P"):-1,("E","S"):0,
    ("E","T"):-1,("E","W"):-3,("E","Y"):-2,("E","V"):-2,
    ("G","G"):6, ("G","H"):-2,("G","I"):-4,("G","L"):-4,("G","K"):-2,
    ("G","M"):-3,("G","F"):-3,("G","P"):-2,("G","S"):0, ("G","T"):-2,
    ("G","W"):-2,("G","Y"):-3,("G","V"):-3,
    ("H","H"):8, ("H","I"):-3,("H","L"):-3,("H","K"):-1,("H","M"):-2,
    ("H","F"):-1,("H","P"):-2,("H","S"):-1,("H","T"):-2,("H","W"):-2,
    ("H","Y"):2, ("H","V"):-3,
    ("I","I"):4, ("I","L"):2, ("I","K"):-1,("I","M"):1, ("I","F"):0,
    ("I","P"):-3,("I","S"):-2,("I","T"):-1,("I","W"):-3,("I","Y"):-1,
    ("I","V"):3,
    ("L","L"):4, ("L","K"):-2,("L","M"):2, ("L","F"):0, ("L","P"):-3,
    ("L","S"):-2,("L","T"):-1,("L","W"):-2,("L","Y"):-1,("L","V"):1,
    ("K","K"):5, ("K","M"):-1,("K","F"):-3,("K","P"):-1,("K","S"):0,
    ("K","T"):-1,("K","W"):-3,("K","Y"):-2,("K","V"):-2,
    ("M","M"):5, ("M","F"):0, ("M","P"):-2,("M","S"):-1,("M","T"):-1,
    ("M","W"):-1,("M","Y"):-1,("M","V"):1,
    ("F","F"):6, ("F","P"):-4,("F","S"):-2,("F","T"):-2,("F","W"):1,
    ("F","Y"):3, ("F","V"):-1,
    ("P","P"):7, ("P","S"):-1,("P","T"):-1,("P","W"):-4,("P","Y"):-3,
    ("P","V"):-2,
    ("S","S"):4, ("S","T"):1, ("S","W"):-3,("S","Y"):-2,("S","V"):-2,
    ("T","T"):5, ("T","W"):-2,("T","Y"):-2,("T","V"):0,
    ("W","W"):11,("W","Y"):2, ("W","V"):-3,
    ("Y","Y"):7, ("Y","V"):-1,
    ("V","V"):4,
}


def blosum62(aa1: str, aa2: str) -> int:
    """Return BLOSUM62 score for two amino acids (symmetric lookup)."""
    if aa1 == aa2:
        return BLOSUM62_RAW.get((aa1, aa1), 0)
    key = tuple(sorted([aa1, aa2]))
    return BLOSUM62_RAW.get(key, -4)  # Default penalty for unknown pairs


def safe_mean(values: list) -> float:
    """Return mean of list, or 0.0 if empty."""
    if not values:
        return 0.0
    return sum(values) / len(values)


def evaluate_predictions(predictions: list[dict], critical_positions_0: set[int]) -> dict:
    """Compute all recovery metrics from a list of position-level prediction dicts."""
    all_positions = predictions
    critical_preds = [p for p in predictions if p["position_0indexed"] in critical_positions_0]
    non_critical_preds = [p for p in predictions if p["position_0indexed"] not in critical_positions_0]

    # AARR total
    n_total = len(all_positions)
    n_correct_total = sum(1 for p in all_positions if p["correct"])
    aarr_total = round(n_correct_total / n_total, 4) if n_total > 0 else 0.0

    # AARR critical
    n_critical = len(critical_preds)
    n_correct_critical = sum(1 for p in critical_preds if p["correct"])
    aarr_critical = round(n_correct_critical / n_critical, 4) if n_critical > 0 else 0.0

    # BLOSUM62 scores at mispredicted positions
    blosum_wrong_all = []
    blosum_wrong_critical = []
    conservative_broken_critical = 0
    n_wrong_critical = 0

    for pred in all_positions:
        if not pred["correct"]:
            score = blosum62(pred["original_aa"], pred["predicted_aa"])
            blosum_wrong_all.append(score)
            if pred["position_0indexed"] in critical_positions_0:
                blosum_wrong_critical.append(score)
                n_wrong_critical += 1
                if score < 0:
                    conservative_broken_critical += 1

    avg_blosum62_wrong_all = round(safe_mean(blosum_wrong_all), 4)
    avg_blosum62_wrong_critical = round(safe_mean(blosum_wrong_critical), 4)

    # risk_score: fraction of wrong critical predictions with BLOSUM62 < 0
    risk_score = round(conservative_broken_critical / n_wrong_critical, 4) if n_wrong_critical > 0 else 0.0

    # Confidence metrics
    probs_correct = [p["prob"] for p in all_positions if p["correct"]]
    probs_incorrect = [p["prob"] for p in all_positions if not p["correct"]]
    avg_prob_correct = round(safe_mean(probs_correct), 4)
    avg_prob_incorrect = round(safe_mean(probs_incorrect), 4)

    # Per-position details for critical residues
    critical_details = []
    for pred in critical_preds:
        score = blosum62(pred["original_aa"], pred["predicted_aa"]) if not pred["correct"] else None
        critical_details.append({
            "position_1indexed": pred["position_1indexed"],
            "original_aa": pred["original_aa"],
            "predicted_aa": pred["predicted_aa"],
            "correct": pred["correct"],
            "prob": pred["prob"],
            "blosum62": score if score is not None else blosum62(pred["original_aa"], pred["original_aa"]),
            "conservative_change": (score is not None and score >= 0) if score is not None else True,
        })

    return {
        "aarr_total": aarr_total,
        "aarr_critical": aarr_critical,
        "n_masked_total": n_total,
        "n_masked_critical": n_critical,
        "n_correct_total": n_correct_total,
        "n_correct_critical": n_correct_critical,
        "avg_blosum62_wrong_all": avg_blosum62_wrong_all,
        "avg_blosum62_wrong_critical": avg_blosum62_wrong_critical,
        "risk_score": risk_score,
        "n_conservative_broken_critical": conservative_broken_critical,
        "avg_prob_correct": avg_prob_correct,
        "avg_prob_incorrect": avg_prob_incorrect,
        "critical_position_details": critical_details,
    }


def main():
    print("=" * 60)
    print("T3.1 Protein Mask & Fill - Evaluate recovery metrics")
    print("=" * 60)

    pred_files = sorted(MASKFILL_DIR.glob("*_predictions.json"))
    if not pred_files:
        print(f"[ERROR] No prediction files found in {MASKFILL_DIR}")
        print("Run scripts/run_esm2_maskfill.py first.")
        sys.exit(1)

    print(f"[INFO] Found {len(pred_files)} prediction files.")

    recovery_summary = []

    for pred_file in pred_files:
        with open(pred_file) as f:
            data = json.load(f)

        protein_name = data["protein_name"]
        strategy = data["strategy"]
        critical_positions_1 = data["critical_positions_1indexed"]
        critical_positions_0 = {p - 1 for p in critical_positions_1}
        predictions = data["predictions"]

        print(f"\n  [{protein_name}] [{strategy}]")

        metrics = evaluate_predictions(predictions, critical_positions_0)

        entry = {
            "protein_name": protein_name,
            "uniprot": data["uniprot"],
            "pdb": data["pdb"],
            "risk_type": data["risk_type"],
            "strategy": strategy,
            "seq_length": data["seq_length"],
            **metrics,
            "model_used": data.get("model_used", "unknown"),
            "using_mock": data.get("using_mock", False),
        }
        recovery_summary.append(entry)

        print(f"    AARR_total={metrics['aarr_total']:.3f}  "
              f"AARR_critical={metrics['aarr_critical']:.3f}  "
              f"risk_score={metrics['risk_score']:.3f}  "
              f"avg_BLOSUM62={metrics['avg_blosum62_wrong_critical']:.2f}")

    # Save summary
    out_path = RECOVERY_DIR / "recovery_summary.json"
    with open(out_path, "w") as f:
        json.dump(recovery_summary, f, indent=2)
    print(f"\n[OK] Recovery summary saved: {out_path}")

    # Per-protein aggregated stats
    protein_agg = {}
    for entry in recovery_summary:
        pname = entry["protein_name"]
        if pname not in protein_agg:
            protein_agg[pname] = {
                "protein_name": pname,
                "risk_type": entry["risk_type"],
                "aarr_critical_values": [],
                "risk_score_values": [],
            }
        protein_agg[pname]["aarr_critical_values"].append(entry["aarr_critical"])
        protein_agg[pname]["risk_score_values"].append(entry["risk_score"])

    protein_stats = []
    for pname, agg in protein_agg.items():
        vals_aarr = agg["aarr_critical_values"]
        vals_risk = agg["risk_score_values"]
        protein_stats.append({
            "protein_name": pname,
            "risk_type": agg["risk_type"],
            "mean_aarr_critical": round(safe_mean(vals_aarr), 4),
            "min_aarr_critical": round(min(vals_aarr), 4),
            "max_aarr_critical": round(max(vals_aarr), 4),
            "mean_risk_score": round(safe_mean(vals_risk), 4),
        })
    protein_stats.sort(key=lambda x: x["mean_aarr_critical"])  # lowest = riskiest

    protein_stats_path = RECOVERY_DIR / "protein_stats.json"
    with open(protein_stats_path, "w") as f:
        json.dump(protein_stats, f, indent=2)
    print(f"[OK] Per-protein stats saved: {protein_stats_path}")

    # Print final table
    print("\n" + "=" * 90)
    print(f"{'Protein':<30} {'Strategy':<22} {'AARR_tot':>8} {'AARR_crit':>9} {'Risk':>6} {'BLOSUM':>7}")
    print("-" * 90)
    for entry in recovery_summary:
        print(
            f"{entry['protein_name']:<30} "
            f"{entry['strategy']:<22} "
            f"{entry['aarr_total']:>8.3f} "
            f"{entry['aarr_critical']:>9.3f} "
            f"{entry['risk_score']:>6.3f} "
            f"{entry['avg_blosum62_wrong_critical']:>7.2f}"
        )
    print("=" * 90)


if __name__ == "__main__":
    main()
