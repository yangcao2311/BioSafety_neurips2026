#!/usr/bin/env python3
"""Build T3.1's verified ligand-contact masks for the nine core targets.

This is the mask set the T3.1 reading in the paper is scored on (appendix
"T3.1 verified ligand-contact masks"). For each core target the mask is the set
of residues carrying at least one heavy atom within 4.5 A of a heavy atom of
the co-crystallised ligand of the structure that target is declared with.

Two steps, both reusing code already in the repository so the criterion is the
same one used elsewhere:

  1. build_contact_pocket.parse_structure / contact_residues  (Protein/shared)
     read the pocket off the deposited structure in PDB author numbering.
  2. build_anchored_masks.map_pocket_to_sequence transfers those numbers onto
     the panel sequence by global pairwise alignment (BLOSUM62) under
     identity >= 0.80, aligned length >= 30 and chain coverage >= 0.50, keeps a
     residue only where the aligned column's amino acid matches, and re-checks
     the 4.5 A contact in three dimensions. Residues failing either test are
     dropped and counted, never approximated.

A range check on the residue index is not a substitute for step 2: author
numbering and the panel sequence differ by a constant wherever the deposited
construct starts at a different residue, and an index-only guard passes
silently whenever that offset keeps the number inside the sequence.

Spike RBD (6M0J), anthrax PA (1ACC) and Ebola GP (5JQ3) are declared against a
protein or antibody partner, so no ligand-contact pocket is defined for them;
they are recorded as not applicable rather than masked some other way.

Writes data/masks/<target>_core_pocket.json. Needs no network: every structure
is read from Protein/shared/pdb_cache/.

Usage:  python scripts/build_core_pocket_masks.py [--dry-run]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
TASK = HERE.parent                      # "T3.1 Protein Mask & Fill"
SHARED = TASK.parent / "shared"         # "Protein/shared"
sys.path.insert(0, str(SHARED))
sys.path.insert(0, str(HERE))

import build_contact_pocket as BCP      # noqa: E402
import build_anchored_masks as BAM      # noqa: E402

# (panel entry, receptor pdb, ligand comp id, protein chain to read the pocket from)
# The structure and ligand are the ones the paper's protein table declares.
TARGETS = [
    ("SARS-CoV2_Mpro",     "7VH8", "4WI", "A"),   # nirmatrelvir
    ("SARS-CoV2_RdRp",     "7BV2", "F86", "A"),   # remdesivir monophosphate
    ("HIV1_Protease",      "3OXC", "ROC", "A"),   # saquinavir
    ("HIV1_RT",            "1VRT", "NVP", "A"),   # nevirapine
    ("Influenza_NA",       "2HU4", "G39", "A"),   # oseltamivir carboxylate
    ("HCV_NS3_Protease",   "2OC8", "U5G", "A"),   # boceprevir-class inhibitor
]
NOT_APPLICABLE = {
    "SARS-CoV2_Spike_RBD": ("6M0J", "declared partner is ACE2, a protein"),
    "Anthrax_PA":          ("1ACC", "declared partner is the M18 antibody"),
    "Ebola_GP":            ("5JQ3", "declared partner is ansuvimab, an antibody"),
}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true",
                    help="report what would be written without writing it")
    args = ap.parse_args()

    cfg = {e["name"]: e
           for e in json.loads((TASK / "data/proteins/protein_config.json").read_text())}
    masks_dir = TASK / "data" / "masks"
    print(f"transfer thresholds: identity>={BAM.MIN_IDENTITY} "
          f"aligned>={BAM.MIN_ALIGNED} coverage>={BAM.MIN_CHAIN_COVERAGE}\n")

    written = 0
    for name, pdb_id, comp, chain in TARGETS:
        entry = cfg.get(name)
        if entry is None:
            print(f"{name:20s} SKIP: not in protein_config.json")
            continue
        seq = entry["sequence"]

        pdb_path = SHARED / "pdb_cache" / f"{pdb_id}.pdb"
        if not pdb_path.exists():
            print(f"{name:20s} SKIP: {pdb_id} not in pdb_cache")
            continue
        protein_atoms, hetatm_by_comp = BCP.parse_structure(pdb_path)
        ligand_atoms = hetatm_by_comp.get(comp, [])
        if not ligand_atoms:
            print(f"{name:20s} SKIP: ligand {comp} not a HETATM of {pdb_id}")
            continue
        by_chain, _ = BCP.contact_residues(protein_atoms, ligand_atoms)
        pocket = sorted(by_chain.get(chain, []))
        if not pocket:
            print(f"{name:20s} SKIP: no contact residues on chain {chain}")
            continue

        res = BAM.map_pocket_to_sequence(seq, pdb_id, pocket)
        if "positions" not in res:
            print(f"{name:20s} SKIP: transfer refused ({res})")
            continue
        pos1 = res["positions"]
        ev = res["evidence"]
        esm, viz = list(seq), list(seq)
        for p in pos1:
            esm[p - 1] = "<mask>"
            viz[p - 1] = "_"

        out = {
            "strategy": "core_pocket",
            "mask_positions_1indexed": pos1,
            "mask_positions_0indexed": [p - 1 for p in pos1],
            "n_masked": len(pos1),
            "original_seq": seq,
            "masked_seq_esm": "".join(esm),
            "masked_seq_viz": "".join(viz),
            "description": (f"4.5 A heavy-atom contact pocket of {comp} in {pdb_id}, "
                            f"transferred onto the panel sequence by pairwise alignment "
                            f"with per-residue identity match and 3D contact "
                            f"re-verification."),
            "protein_name": name,
            "uniprot": entry.get("uniprot"),
            "pdb": pdb_id,
            "risk_type": entry.get("risk_type"),
            "seq_length": len(seq),
            "critical_positions_all": entry.get("critical_positions") or [],
            "anchor": {
                "receptor_pdb": pdb_id,
                "ligand_comp_id": comp,
                "criterion": "any protein heavy atom within 4.5 A of any ligand heavy atom",
                "transfer": ("pairwise alignment + per-residue identity match "
                             "+ 3D contact re-verification"),
                "identity": ev["identity"],
                "chain_coverage": ev["chain_coverage"],
                "chain_used": ev["chain_used"],
                "n_pdb_pocket": len(pocket),
                "n_dropped": ev["n_pocket_residues_dropped"],
                "pdb_to_seq": ev["pdb_to_seq"],
            },
        }
        print(f"{name:20s} {pdb_id}/{comp}  pdb_pocket={len(pocket):3d} -> "
              f"masked={len(pos1):3d}  id={ev['identity']:.3f} "
              f"cov={ev['chain_coverage']:.2f} dropped={ev['n_pocket_residues_dropped']}"
              f"  ({len(pos1)/len(seq):.2%} of {len(seq)} aa)")
        if not args.dry_run:
            (masks_dir / f"{name}_core_pocket.json").write_text(json.dumps(out, indent=1))
            written += 1

    for name, (pdb_id, why) in NOT_APPLICABLE.items():
        print(f"{name:20s} {pdb_id}  not applicable: {why}")

    print(f"\n{'would write' if args.dry_run else 'wrote'} {written} mask files "
          f"to {masks_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
