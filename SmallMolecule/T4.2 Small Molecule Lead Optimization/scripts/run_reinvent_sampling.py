#!/usr/bin/env python3
"""
T4.2 REINVENT 4 De Novo Sampling baseline.

Sample N molecules from the bare reinvent.prior (ChEMBL-25-trained RNN),
no scaffold conditioning, as Lead Optimization baseline.
"""
import json
import subprocess
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
RESULTS = BASE / "results" / "sampling"
RESULTS.mkdir(parents=True, exist_ok=True)

PRIOR = Path("/home/xliu316/bio/REINVENT4/priors/reinvent.prior")
NUM_SMILES = 600


def main():
    print("=== T4.2 REINVENT 4 De Novo Sampling (baseline) ===")
    out_csv = RESULTS / "denovo_sampling.csv"
    toml_path = RESULTS / "sampling.toml"
    toml_path.write_text(f"""run_type = "sampling"
device = "cuda:0"

[parameters]
model_file = "{PRIOR}"
output_file = "{out_csv}"
num_smiles = {NUM_SMILES}
unique_molecules = true
randomize_smiles = true
""")
    cmd = ["reinvent", str(toml_path)]
    res = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    if res.returncode != 0:
        print(f"REINVENT error: {res.stderr[-500:]}")

    if not out_csv.exists():
        raise SystemExit("No output CSV produced")

    import csv
    records = []
    with open(out_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            smi = row.get("SMILES", "").strip()
            if not smi:
                continue
            records.append({
                "generated_smiles": smi,
                "nll": float(row.get("NLL", 0)) if row.get("NLL") else 0.0,
                "model": "reinvent_denovo",
            })

    out_json = RESULTS / "denovo_generated.json"
    json.dump({
        "method": "REINVENT 4 de novo sampling (reinvent.prior, ChEMBL-25)",
        "n_total_generated": len(records),
        "records": records,
    }, open(out_json, "w"), indent=2)
    print(f"Saved {len(records)} de novo SMILES → {out_json}")


if __name__ == "__main__":
    main()
