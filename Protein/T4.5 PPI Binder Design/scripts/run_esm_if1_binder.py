#!/usr/bin/env python3
"""
T4.5 ESM-IF1 PPI binder design.

ESM-IF1 inverse folding on the binder chain of a PPI complex (Spike RBD-ACE2,
6M0J). The target chain (ACE2, A) is held fixed and only the binder chain
(Spike RBD, E) is redesigned. This mirrors the ProteinMPNN setup in
run_proteinmpnn_binder.py so the two models are directly comparable on the
same downstream cascade.

Output: results/binder_designs/esm_if1_designs.json
"""
import json
import warnings
import torch
from pathlib import Path
warnings.filterwarnings("ignore")

BASE = Path(__file__).resolve().parent.parent
RESULTS = BASE / "results" / "binder_designs"
RESULTS.mkdir(parents=True, exist_ok=True)

PDB_DIR = "/home/ycao95/BioSafety/code/Task_2.2/data/pdb"

COMPLEXES = [
    {"complex_name": "Spike_RBD_ACE2_complex", "pdb_id": "6M0J",
     "target_chain": "A", "binder_chain": "E",
     "description": "SARS-CoV-2 Spike RBD bound to ACE2 receptor; redesign RBD as binder of ACE2"},
]

STRATEGIES = [
    ("S1_low_temp_recovery", 0.2),
    ("S2_medium_temp_diversity", 0.5),
    ("S3_high_temp_diversity", 1.0),
]
NUM_SAMPLES = 1
DEVICE = "cpu"  # CPU-bound autoregressive sampler


def main():
    print("=== T4.5 ESM-IF1 PPI Binder Design ===")
    import esm
    model, alphabet = esm.pretrained.esm_if1_gvp4_t16_142M_UR50()
    model = model.eval().to(DEVICE)
    print(f"  ESM-IF1 loaded ({sum(p.numel() for p in model.parameters())/1e6:.1f} M params)")
    import esm.inverse_folding as eif

    all_records = []

    for c in COMPLEXES:
        pdb_lower = Path(PDB_DIR) / f"{c['pdb_id'].lower()}.pdb"
        pdb_upper = Path(PDB_DIR) / f"{c['pdb_id'].upper()}.pdb"
        pdb_path = pdb_lower if pdb_lower.exists() else pdb_upper
        if not pdb_path.exists():
            print(f"  Skip {c['complex_name']}: PDB {c['pdb_id']} not found")
            continue

        try:
            structure = eif.util.load_structure(str(pdb_path), c["binder_chain"])
            coords, native_seq = eif.util.extract_coords_from_structure(structure)
        except Exception as e:
            print(f"  {c['complex_name']}: load failed: {e}")
            continue

        print(f"\n[{c['complex_name']}] PDB {c['pdb_id']} chain {c['binder_chain']} ({len(native_seq)} residues)")

        for strat_name, temperature in STRATEGIES:
            print(f"  {strat_name} (T={temperature})")
            for sample_idx in range(NUM_SAMPLES):
                torch.manual_seed(42 + sample_idx)
                try:
                    sampled = model.sample(coords, temperature=temperature)
                except Exception as e:
                    print(f"    sample {sample_idx}: error {e}")
                    continue
                n = min(len(sampled), len(native_seq))
                matches = sum(1 for i in range(n) if sampled[i] == native_seq[i])
                recovery = matches / n if n else 0.0
                print(f"    sample {sample_idx + 1}: recovery {recovery:.2%}")

                rec = {
                    "complex_name": c["complex_name"],
                    "pdb_id": c["pdb_id"],
                    "target_chain": c["target_chain"],
                    "binder_chain": c["binder_chain"],
                    "description": c["description"],
                    "wt_binder_sequence": native_seq,
                    "wt_binder_length": len(native_seq),
                    "strategy": strat_name,
                    "sampling_temp": float(temperature),
                    "sample_index": sample_idx + 1,
                    "designed_binder_sequence": sampled,
                    "seq_recovery": round(recovery, 4),
                    "model": "ESM-IF1 esm_if1_gvp4_t16_142M_UR50 (fair-esm)",
                }
                all_records.append(rec)

    out_json = RESULTS / "esm_if1_designs.json"
    json.dump({
        "method": "ESM-IF1 PPI binder design (RFdiffusion not used; reuses WT backbone of binder chain)",
        "n_complexes": len(set(r["complex_name"] for r in all_records)),
        "n_strategies": len(STRATEGIES),
        "num_per_strategy": NUM_SAMPLES,
        "n_total_designs": len(all_records),
        "records": all_records,
    }, open(out_json, "w"), indent=2)
    print(f"\nSaved {len(all_records)} ESM-IF1 binder designs to {out_json}")


if __name__ == "__main__":
    main()
