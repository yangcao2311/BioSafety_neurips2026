#!/usr/bin/env python3
"""
T3.2 ESM-IF1 inverse-folding design.

ESM-IF1 (esm_if1_gvp4_t16_142M_UR50, 141.7M params) is a graph neural network
based inverse-folding model. Given a protein backbone (PDB), it samples
candidate amino-acid sequences that fit the structure.

This script runs ESM-IF1 on the 9 biosafety proteins, mirroring the
ProteinMPNN runner (run_proteinmpnn_design.py) so the cascade can score
both models with the same downstream pipeline.

Strategies:
  S1_low_temp_recovery     temperature 0.2  conservative
  S2_medium_temp_diversity temperature 0.5  moderate
  S3_high_temp_diversity   temperature 1.0  aggressive

Output: results/esm_if1/esm_if1_designs.json
"""
import json
import warnings
import numpy as np
import torch
from pathlib import Path
warnings.filterwarnings("ignore")

BASE = Path(__file__).resolve().parent.parent
RESULTS = BASE / "results" / "esm_if1"
RESULTS.mkdir(parents=True, exist_ok=True)

PROTEIN_INFO = "/home/ycao95/BioSafety/code/Task_2.2/data/tier_a_proteins.json"
PDB_DIR = "/home/ycao95/BioSafety/code/Task_2.2/data/pdb"

STRATEGIES = [
    ("S1_low_temp_recovery", 0.2),
    ("S2_medium_temp_diversity", 0.5),
    ("S3_high_temp_diversity", 1.0),
]
NUM_SAMPLES = 1  # CPU sampling is slow; one sample per (protein, strategy) gives 27 designs total
DEVICE = "cpu"  # ESM-IF1 sampling has CPU-bound autoregressive loop; running fully on CPU avoids device mismatch


def load_pdb(pdb_path, chain_id="A"):
    """Load CA-N-C backbone coords from a single chain of a PDB file."""
    import esm.inverse_folding as eif
    structure = eif.util.load_structure(pdb_path, chain_id)
    coords, seq = eif.util.extract_coords_from_structure(structure)
    return coords, seq


def sample_one(model, alphabet, coords, temperature, seed=42):
    sampled = model.sample(coords, temperature=temperature)
    return sampled


def main():
    print("=== T3.2 ESM-IF1 Inverse-Folding Design ===")
    print(f"Loading ESM-IF1 ...")
    import esm
    model, alphabet = esm.pretrained.esm_if1_gvp4_t16_142M_UR50()
    model = model.eval().to(DEVICE)
    print(f"  ESM-IF1 loaded ({sum(p.numel() for p in model.parameters())/1e6:.1f} M params)")

    proteins = json.load(open(PROTEIN_INFO))
    all_records = []

    for p in proteins:
        short = p["short_name"]
        pdb_id = p["pdb_id"]
        pdb_lower = Path(PDB_DIR) / f"{pdb_id.lower()}.pdb"
        pdb_upper = Path(PDB_DIR) / f"{pdb_id.upper()}.pdb"
        pdb_path = pdb_lower if pdb_lower.exists() else pdb_upper
        if not pdb_path.exists():
            print(f"  Skip {short}: PDB {pdb_id} missing")
            continue

        try:
            coords, native_seq = load_pdb(str(pdb_path), chain_id="A")
        except Exception as e:
            print(f"  {short}: load failed: {e}")
            continue

        print(f"\n[{short}] PDB {pdb_id} chain A loaded ({len(native_seq)} residues)")

        for strat_name, temperature in STRATEGIES:
            print(f"  {strat_name} (T={temperature})")
            for sample_idx in range(NUM_SAMPLES):
                torch.manual_seed(42 + sample_idx)
                np.random.seed(42 + sample_idx)
                try:
                    sampled = sample_one(model, alphabet, coords, temperature)
                except Exception as e:
                    print(f"    sample {sample_idx}: error {e}")
                    continue

                n = min(len(sampled), len(native_seq))
                matches = sum(1 for i in range(n) if sampled[i] == native_seq[i])
                recovery = matches / n if n else 0.0
                print(f"    sample {sample_idx + 1}: recovery {recovery:.2%}")

                rec = {
                    "protein_short_name": short,
                    "uniprot_id": p["uniprot_id"],
                    "pdb_id": pdb_id,
                    "wt_sequence": native_seq,
                    "wt_length": len(native_seq),
                    "functional_positions": p.get("functional_positions", []),
                    "strategy": strat_name,
                    "sampling_temp": float(temperature),
                    "sample_index": sample_idx + 1,
                    "designed_sequence": sampled,
                    "designed_length": len(sampled),
                    "seq_recovery": round(recovery, 4),
                    "model": "ESM-IF1 esm_if1_gvp4_t16_142M_UR50 (fair-esm)",
                }
                all_records.append(rec)

    out_json = RESULTS / "esm_if1_designs.json"
    json.dump({
        "method": "ESM-IF1 graph neural network inverse folding (esm_if1_gvp4_t16_142M_UR50)",
        "n_proteins": len(set(r["protein_short_name"] for r in all_records)),
        "n_strategies": len(STRATEGIES),
        "num_seq_per_strategy": NUM_SAMPLES,
        "n_total_designs": len(all_records),
        "records": all_records,
    }, open(out_json, "w"), indent=2)
    print(f"\nSaved {len(all_records)} ESM-IF1 designs to {out_json}")


if __name__ == "__main__":
    main()
