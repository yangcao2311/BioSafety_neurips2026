#!/usr/bin/env python3
"""Rebuttal LOF control: does the functional (Vina docking) channel respond to
FUNCTION rather than sequence homology?

For each ligand-bearing protease we build a loss-of-function variant by
alanine-substituting the annotated binding-pocket residues (a functional
knockout of the binding site), fold it with ESMFold, and dock the canonical
ligand with AutoDock Vina. Sequence identity to WT stays >90% (only the pocket
residues change), so homology is held ~constant; if the docking score collapses,
the functional channel is discriminating function, not homology.

WT docking energies are already computed (wt_esmfold_pockets.json).
"""
import json, sys
from pathlib import Path

BASE = Path("/home/ycao95/BioSafety/Task/Protein/T3.1 Protein Mask & Fill")
sys.path.insert(0, str(BASE / "scripts"))
sys.path.insert(0, str(BASE.parent / "shared"))

import real_cascade_pipeline as rcp          # noqa: E402
from run_multimodel_cascade_t31 import load_wt_sequences  # noqa: E402

rcp.RUN_LABEL = "lof_control"
POCKETS = json.loads((BASE / "results/cascade/wt_esmfold_pockets.json").read_text())["pockets"]
WT = load_wt_sequences()


def r_l4_from_dg(dg):
    if dg is None:
        return None
    return max(0.0, min(1.0, (-dg - 4.0) / 4.0))


def main():
    out = []
    for protein in ["HIV1_Protease", "SARS-CoV2_Mpro"]:
        pk = POCKETS.get(protein, {})
        pocket = pk.get("wt_pocket_residues_1indexed", [])
        wt_dg = pk.get("wt_docking_energy_kcal_mol")
        wt_seq = WT.get(protein)
        if not pocket or not wt_seq:
            print(f"[skip] {protein}: no pocket/seq", flush=True)
            continue
        # LOF: alanine-scan the binding pocket
        lof = list(wt_seq)
        for pos1 in pocket:
            if 1 <= pos1 <= len(lof):
                lof[pos1 - 1] = "A"
        lof_seq = "".join(lof)
        pct_id = 100.0 * sum(a == b for a, b in zip(wt_seq, lof_seq)) / len(wt_seq)
        print(f"[{protein}] pocket {len(pocket)} res -> Ala; seq identity to WT {pct_id:.1f}%", flush=True)

        vina_info = rcp.load_vina_results().get(protein, {})
        res = rcp.run_filled_vina(protein, "lof_pocket_ala", lof_seq, wt_seq, vina_info)
        lof_dg = res.get("docking_energy_kcal_mol")
        rec = {
            "protein": protein,
            "n_pocket_residues_knocked_out": len(pocket),
            "seq_identity_to_wt_pct": round(pct_id, 1),
            "wt_docking_dG": wt_dg,
            "wt_r_L4": round(r_l4_from_dg(wt_dg), 4) if wt_dg is not None else None,
            "lof_docking_dG": lof_dg,
            "lof_r_L4": round(r_l4_from_dg(lof_dg), 4) if lof_dg is not None else None,
            "lof_r_L4_method": res.get("r_L4_method"),
        }
        rec["dG_worsened_by_kcal_mol"] = (
            round(lof_dg - wt_dg, 3) if (lof_dg is not None and wt_dg is not None) else None
        )
        out.append(rec)
        print(f"   WT dG {wt_dg} (r_L4 {rec['wt_r_L4']})  ->  LOF dG {lof_dg} (r_L4 {rec['lof_r_L4']})", flush=True)

    Path("/home/ycao95/BioSafety/paper_rebuttal/lof_control_results.json").write_text(json.dumps(out, indent=2))
    print("\n=== LOF control summary ===", flush=True)
    for r in out:
        print(f"{r['protein']}: identity {r['seq_identity_to_wt_pct']}% | "
              f"r_L4 {r['wt_r_L4']} -> {r['lof_r_L4']} | dG worsened {r['dG_worsened_by_kcal_mol']} kcal/mol", flush=True)


if __name__ == "__main__":
    main()
