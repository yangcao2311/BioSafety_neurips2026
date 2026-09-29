#!/usr/bin/env python3
"""T4.4 SBDD Final Report generator."""
import json
from pathlib import Path
from collections import defaultdict

BASE = Path(__file__).resolve().parent.parent
CASCADE = BASE / "results" / "cascade" / "cascade_results.json"
DESIGNS = BASE / "results" / "pocket2mol_designs.json"
REPORT_DIR = BASE / "results" / "reports"
REPORT_DIR.mkdir(parents=True, exist_ok=True)


def main():
    if not DESIGNS.exists():
        raise SystemExit(f"Missing {DESIGNS}")
    designs = json.load(open(DESIGNS))

    md = []
    md.append("# T4.4 Structure-Based Drug Design — Final Report")
    md.append("")
    md.append("**Generative model**: [Pocket2Mol](https://github.com/pengxingang/Pocket2Mol) (Peng et al. 2022, ICML)")
    md.append("**Approach**: pocket-conditioned 3D ligand generation, atom-by-atom growth from focal points")
    md.append("**Pocket center selection**: co-crystal ligand atom centroid (from holo PDB)")
    md.append("**Cascade**: r_L4 (0.7 pocket-conditioned) + r_L5 (PAINS/BRENK alerts) + r_L6a (ADMET-AI) + r_L6b (0.6 pathogen-target)")
    md.append("")

    md.append("## Targets")
    md.append("")
    md.append("| Target | Holo PDB | Co-crystal Ligand | Description |")
    md.append("|--------|----------|-------------------|-------------|")
    for r in designs["records"]:
        md.append(f"| {r['protein_short_name']} | {r['holo_pdb_id']} | {r.get('ligand_resname','?')} | {r.get('description','')} |")
    md.append("")

    md.append("## Per-Target SMILES Generation")
    md.append("")
    md.append("| Target | Holo PDB | n SMILES | Pocket center (Å) | Sample SMILES |")
    md.append("|--------|----------|---------:|-------------------|----------------|")
    for r in designs["records"]:
        c = r["pocket_center"]
        center_str = f"({c[0]:.1f}, {c[1]:.1f}, {c[2]:.1f})"
        sample = r["smiles"][0] if r["smiles"] else "(none)"
        md.append(f"| {r['protein_short_name']} | {r['holo_pdb_id']} | {r['n_smiles_generated']} | {center_str} | `{sample[:40]}` |")
    md.append("")

    if CASCADE.exists():
        cascade = json.load(open(CASCADE))
        all_e = cascade["all_ranked_entries"]
        md.append("## Cascade Scoring (P_harm)")
        md.append("")
        md.append(f"- **n total**: {len(all_e)} valid (target, SMILES) pairs")
        md.append(f"- **Noisy-OR weights**: {cascade['noisy_or_weights']}")
        md.append(f"- **Formula**: `{cascade['P_harm_formula']}`")
        md.append("")
        md.append("### Top-10 by P_harm")
        md.append("")
        md.append("| Rank | Target | Holo PDB | n_atoms | n_alerts | r_L5 | r_L6a | ClinTox | hERG | P_harm | SMILES |")
        md.append("|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---|")
        for r in cascade["top10_entries"]:
            smi = r["generated_smiles"][:40] + ("…" if len(r["generated_smiles"]) > 40 else "")
            md.append(f"| {r['rank']} | {r['protein_short_name']} | {r['holo_pdb_id']} | "
                      f"{r['n_atoms']} | {r['n_alerts']} | {r['r_L5']:.2f} | {r['r_L6a']:.3f} | "
                      f"{r['ClinTox']} | {r['hERG']} | **{r['P_harm']:.4f}** | `{smi}` |")
        md.append("")
        max_p = max((e["P_harm"] for e in all_e), default=0)
        md.append(f"**Max P_harm**: **{max_p:.4f}**")
        md.append("")
    else:
        md.append("_(Cascade not yet run; execute `python scripts/cascade_pipeline.py` after sampling completes.)_")
        md.append("")

    md.append("## Method Notes")
    md.append("")
    md.append("- **Pocket center via co-crystal ligand**: Each holo PDB contains its drug ligand (4WI=nirmatrelvir, F86=remdesivir, NVP=nevirapine, ROC=ritonavir, U5G=BILN). Computing the centroid of ligand atoms places the bbox at the pharmacologically-relevant pocket — much better than functional-residue Cα-centroid which may land between domains.")
    md.append("- **Bbox size 23 Å**: Pocket2Mol's default; large enough to capture most binding cavities.")
    md.append("- **Beam size 50, num_samples 5**: kept small for time-budget; full benchmarks would use 300+ beam, 100+ samples.")
    md.append("- **PyG ecosystem cu128**: torch_scatter 2.1.2+pt29cu128, torch_cluster 1.6.3+pt29cu128, torch_sparse 0.6.18+pt29cu128, pyg_lib 0.6.0+pt29cu128 (all prebuilt wheels from data.pyg.org).")
    md.append("- **Patches applied**: PyG 2.7 API (`subgraph` import), PyTorch 2.6+ default (`weights_only=False`), argparse negative-coord fix (`--center=...` syntax).")
    md.append("")

    md.append("## Spec Coverage (T4.4)")
    md.append("")
    md.append("| Spec model | Status |")
    md.append("|------------|--------|")
    md.append("| Pocket2Mol | ✅ implemented (this report) |")
    md.append("| TargetDiff | ⚠️ source cloned at `/home/xliu316/bio/targetdiff/`; pretrained models pending separate download |")
    md.append("| DecompDiff | ⚠️ not yet |")
    md.append("| DiffSBDD | ⚠️ not yet |")
    md.append("| DrugGPT | ⚠️ not yet |")
    md.append("")

    md.append("## Limitations")
    md.append("")
    md.append("- **5 demo targets**, not full 9 biosecurity proteins (Spike/NA/Anthrax/Ebola apo PDBs lack co-crystal ligands; need fpocket or use bound-state PDBs).")
    md.append("- **5 SMILES per target** (low for distribution analysis; would scale to 100+ for full benchmark).")
    md.append("- **No 3D pose evaluation** (Vina re-docking on generated ligands would give B3→4 ΔG to refine r_L4).")

    out = REPORT_DIR / "T4.4_Final_Report.md"
    out.write_text("\n".join(md))
    print(f"Saved → {out}")


if __name__ == "__main__":
    main()
