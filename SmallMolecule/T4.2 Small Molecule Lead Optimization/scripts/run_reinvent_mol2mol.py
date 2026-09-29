#!/usr/bin/env python3
"""
T4.2 REINVENT 4 Mol2Mol Lead Optimization.

For each high-tox seed from T4.1/T4.3, run REINVENT mol2mol (medium similarity)
to generate analog molecules. Outputs CSV per seed and a unified JSON.

Mol2Mol mode is generative chemical FM lead-optimization: given a parent SMILES,
sample structurally similar variants using the pretrained transformer prior.
"""
import json
import os
import subprocess
import tempfile
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
DATA = BASE / "data"
RESULTS = BASE / "results" / "mol2mol"
RESULTS.mkdir(parents=True, exist_ok=True)

PRIOR = Path("/home/xliu316/bio/REINVENT4/priors/mol2mol_medium_similarity.prior")
SEEDS_FILE = DATA / "lead_opt_seeds.json"

NUM_PER_SEED = 30           # number of analogs per seed
SAMPLE_STRATEGY = "multinomial"
TEMPERATURE = 1.0


def make_toml(seed_smi: str, output_csv: Path) -> str:
    """Generate a REINVENT TOML config for mol2mol sampling on one seed."""
    smi_path = output_csv.with_suffix(".smi")
    smi_path.write_text(seed_smi.strip() + "\n")
    return f"""run_type = "sampling"
device = "cuda:0"

[parameters]
model_file = "{PRIOR}"
smiles_file = "{smi_path}"
sample_strategy = "{SAMPLE_STRATEGY}"
temperature = {TEMPERATURE}
output_file = "{output_csv}"
num_smiles = {NUM_PER_SEED}
unique_molecules = true
randomize_smiles = true
"""


def run_reinvent_one(seed_smi: str, seed_name: str) -> list:
    """Run REINVENT on one seed; return list of generated dict records."""
    out_csv = RESULTS / f"{seed_name}.csv"
    toml_path = out_csv.with_suffix(".toml")
    toml_path.write_text(make_toml(seed_smi, out_csv))

    cmd = ["reinvent", str(toml_path)]
    res = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
    if res.returncode != 0 and not out_csv.exists():
        print(f"  REINVENT error on {seed_name}: {res.stderr[-500:]}", flush=True)
        return []

    records = []
    if out_csv.exists():
        import csv
        with open(out_csv) as f:
            reader = csv.DictReader(f)
            for row in reader:
                smi = row.get("SMILES", "").strip()
                # mol2mol output column may be 'SMILES' or 'Output_SMILES'
                if not smi:
                    smi = row.get("Output_SMILES", "").strip()
                if not smi:
                    continue
                rec = {
                    "seed_name": seed_name,
                    "seed_smiles": seed_smi,
                    "generated_smiles": smi,
                    "input_smiles": row.get("Input_SMILES", "").strip() or seed_smi,
                    "nll": float(row.get("NLL", 0)) if row.get("NLL") else 0.0,
                    "tanimoto": float(row.get("Tanimoto", 0)) if row.get("Tanimoto") else 0.0,
                    "model": "mol2mol_medium_similarity",
                }
                records.append(rec)
    return records


def main():
    print("=== T4.2 REINVENT 4 Mol2Mol Lead Optimization ===")
    print(f"Prior: {PRIOR.name}")
    print(f"Seeds: {SEEDS_FILE}")

    if not PRIOR.exists():
        raise SystemExit(f"Prior not found: {PRIOR}")
    if not SEEDS_FILE.exists():
        raise SystemExit(f"Seeds JSON not found: {SEEDS_FILE}. Run prepare_seeds.py first.")

    data = json.load(open(SEEDS_FILE))
    seeds = data["seeds"]
    print(f"\nProcessing {len(seeds)} seeds with {NUM_PER_SEED} analogs each...")

    all_records = []
    for i, s in enumerate(seeds, 1):
        seed_smi = s["seed_smiles"]
        # safe identifier for filename
        seed_name = f"{s['source']}_{s['compound_name']}_{i:02d}".replace("/", "_").replace(" ", "_")
        seed_name = "".join(c for c in seed_name if c.isalnum() or c in "_-")[:60]
        print(f"  [{i:2d}/{len(seeds)}] {seed_name} | seed P_harm={s['P_harm_origin']:.4f}")
        recs = run_reinvent_one(seed_smi, seed_name)
        # tag with origin P_harm
        for r in recs:
            r["origin_P_harm"] = s["P_harm_origin"]
            r["origin_compound"] = s["compound_name"]
            r["origin_source"] = s["source"]
        all_records.extend(recs)
        print(f"     → {len(recs)} analogs")

    out_json = BASE / "results" / "mol2mol" / "mol2mol_generated.json"
    json.dump({
        "method": "REINVENT 4 mol2mol_medium_similarity prior",
        "n_seeds": len(seeds),
        "num_per_seed": NUM_PER_SEED,
        "n_total_generated": len(all_records),
        "records": all_records,
    }, open(out_json, "w"), indent=2)
    print(f"\nSaved {len(all_records)} analogs → {out_json}")


if __name__ == "__main__":
    main()
