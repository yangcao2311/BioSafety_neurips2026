#!/usr/bin/env python3
"""
T3.2 ProteinMPNN Inverse-Folding Lead Optimization.

Run ProteinMPNN on the 9 biosecurity PDB structures with three sampling
temperatures (0.2 = high recovery, 0.4 = moderate, 0.6 = high diversity)
to mimic three lead-optimization strategies:

  - S1_low_temp_recovery     (T=0.2): conservative inverse-folding, ~50% AARR
  - S2_medium_temp_diversity (T=0.4): moderate variation
  - S3_high_temp_diversity   (T=0.6): aggressive diversification (more risk)

Output: results/proteinmpnn/proteinmpnn_designs.json — list of design records.
"""
import json
import re
import subprocess
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
DATA_DIR = BASE / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)
RESULTS = BASE / "results" / "proteinmpnn"
RESULTS.mkdir(parents=True, exist_ok=True)

PROTEIN_MPNN_RUN = "/home/ycao95/BioSafety/code/Task_2.2/models/ProteinMPNN/protein_mpnn_run.py"
PDB_DIR = "/home/ycao95/BioSafety/code/Task_2.2/data/pdb"
PROTEIN_INFO = "/home/ycao95/BioSafety/code/Task_2.2/data/tier_a_proteins.json"

STRATEGIES = [
    ("S1_low_temp_recovery", "0.2"),
    ("S2_medium_temp_diversity", "0.4"),
    ("S3_high_temp_diversity", "0.6"),
]

NUM_SEQ_PER_STRAT = 3


def load_protein_info():
    return json.load(open(PROTEIN_INFO))


def parse_proteinmpnn_fasta(fa_path: Path):
    """Parse ProteinMPNN .fa output and return (wt_seq, sampled_seqs).
    The file is:
        >7BQY, score=..., fixed_chains=..., designed_chains=...
        WT_SEQ
        >T=0.2, sample=1, score=..., seq_recovery=0.50
        SAMPLED_SEQ
        ...
    """
    text = fa_path.read_text()
    blocks = text.strip().split("\n>")
    blocks = [b if b.startswith(">") else ">" + b for b in blocks]
    wt_seq = None
    samples = []
    for b in blocks:
        lines = b.strip().splitlines()
        if not lines:
            continue
        header = lines[0].lstrip(">")
        seq = "".join(l for l in lines[1:] if not l.startswith(">"))
        if "T=" not in header:
            wt_seq = seq
        else:
            m_score = re.search(r"score=([\d.]+)", header)
            m_rec = re.search(r"seq_recovery=([\d.]+)", header)
            samples.append({
                "header": header,
                "sequence": seq,
                "score": float(m_score.group(1)) if m_score else None,
                "seq_recovery": float(m_rec.group(1)) if m_rec else None,
            })
    return wt_seq, samples


def run_one(pdb_id: str, pdb_file: Path, strategy_name: str, temp: str, out_dir: Path):
    """Run ProteinMPNN on one PDB at given temperature."""
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "python", PROTEIN_MPNN_RUN,
        "--pdb_path", str(pdb_file),
        "--pdb_path_chains", "A",
        "--out_folder", str(out_dir),
        "--num_seq_per_target", str(NUM_SEQ_PER_STRAT),
        "--sampling_temp", temp,
        "--seed", "42",
        "--batch_size", "1",
        "--suppress_print", "1",
    ]
    res = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    if res.returncode != 0:
        print(f"    ERROR: {res.stderr[-300:]}")
        return None
    fa_dir = out_dir / "seqs"
    fa_files = list(fa_dir.glob("*.fa")) if fa_dir.exists() else []
    if not fa_files:
        return None
    wt_seq, samples = parse_proteinmpnn_fasta(fa_files[0])
    return {"wt_seq": wt_seq, "samples": samples, "fa_file": str(fa_files[0])}


def main():
    info = load_protein_info()
    proteins = []
    for p in info:
        short = p["short_name"]
        pdb_id = p["pdb_id"]
        # find local PDB file (case-insensitive)
        pdb_file_l = Path(PDB_DIR) / f"{pdb_id.lower()}.pdb"
        pdb_file_u = Path(PDB_DIR) / f"{pdb_id.upper()}.pdb"
        if pdb_file_l.exists():
            pdb_file = pdb_file_l
        elif pdb_file_u.exists():
            pdb_file = pdb_file_u
        else:
            print(f"  Skip {short}: PDB {pdb_id} not found")
            continue
        proteins.append({**p, "pdb_file": str(pdb_file)})

    print(f"=== T3.2 ProteinMPNN Lead Optimization ===")
    print(f"Designing {len(proteins)} proteins × {len(STRATEGIES)} strategies "
          f"× {NUM_SEQ_PER_STRAT} sequences = {len(proteins)*len(STRATEGIES)*NUM_SEQ_PER_STRAT} total")

    all_records = []
    for p in proteins:
        short = p["short_name"]
        for strat, temp in STRATEGIES:
            print(f"  [{short}] {strat} (T={temp})")
            out_dir = RESULTS / "raw" / f"{p['pdb_id']}_{strat}"
            r = run_one(p["pdb_id"], Path(p["pdb_file"]), strat, temp, out_dir)
            if r is None:
                print(f"    failed")
                continue
            for i, samp in enumerate(r["samples"]):
                rec = {
                    "protein_short_name": short,
                    "uniprot_id": p["uniprot_id"],
                    "pdb_id": p["pdb_id"],
                    "wt_sequence": p["sequence"],
                    "wt_length": p["seq_len"],
                    "functional_positions": p.get("functional_positions", []),
                    "strategy": strat,
                    "sampling_temp": float(temp),
                    "sample_index": i + 1,
                    "designed_sequence": samp["sequence"],
                    "designed_length": len(samp["sequence"]),
                    "proteinmpnn_score": samp["score"],
                    "seq_recovery": samp["seq_recovery"],
                    "model": "ProteinMPNN v_48_020 (vanilla, full-backbone)",
                }
                all_records.append(rec)
            print(f"    → {len(r['samples'])} designs (avg recovery: "
                  f"{sum(s['seq_recovery'] or 0 for s in r['samples'])/max(1,len(r['samples'])):.2%})")

    out_json = RESULTS / "proteinmpnn_designs.json"
    json.dump({
        "method": "ProteinMPNN v_48_020 (vanilla full-backbone, temp ∈ {0.2, 0.4, 0.6})",
        "n_proteins": len(proteins),
        "n_strategies": len(STRATEGIES),
        "num_seq_per_strategy": NUM_SEQ_PER_STRAT,
        "n_total_designs": len(all_records),
        "records": all_records,
    }, open(out_json, "w"), indent=2)
    print(f"\nSaved {len(all_records)} ProteinMPNN designs → {out_json}")


if __name__ == "__main__":
    main()
