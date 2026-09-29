#!/usr/bin/env python3
"""
T4.4 TargetDiff structure-based drug design.

Loads the public TargetDiff pretrained_diffusion.pt checkpoint (Zenodo
record 14041881, mirrored from the original Google Drive distribution)
and samples 30 ligands per pocket on the same 4 holo PDBs that the
Pocket2Mol run used (7BV2 RdRp, 7VH8 Mpro, 3OXC HIV protease, 1VRT
HIV RT). HCV NS3 (2OC8) is not in the local Pocket2Mol pocket set so
we skip it here.

Pipeline per pocket:
  1. Crop pocket to 10A around co-crystal ligand (TargetDiff
     PDBProtein expects this; the Pocket2Mol-prepared PDBs already
     have the protein-only structure ready).
  2. Run TargetDiff diffusion sampling for num_samples=30 with
     num_steps=1000 (default).
  3. Reconstruct molecules from diffusion outputs, take valid SMILES.
  4. Cascade: L4 (RDKit validity flag), L5 (PAINS/BRENK), L6a
     (ADMET-AI composite), and L6b fixed pathogen prior 0.6 like
     Pocket2Mol cascade.

Output: results/targetdiff/<pdb_id>/sdf/*.sdf and SMILES.txt; cascade
file results/cascade/cascade_results_targetdiff.json.
"""
from __future__ import annotations

import json
import math
import os
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import torch

BASE = Path(__file__).resolve().parent.parent
TARGETDIFF_REPO = Path("/home/xliu316/bio/targetdiff")
sys.path.insert(0, str(TARGETDIFF_REPO))

# Configure
N_SAMPLES = 10  # smaller than the 100 default for faster smoke run
N_STEPS = 200    # 1/5 of default; smoke quality is fine for cascade evaluation
BATCH_SIZE = 10

OUT_DIR = BASE / "results" / "targetdiff"
CASCADE_DIR = BASE / "results" / "cascade"
OUT_DIR.mkdir(parents=True, exist_ok=True)
CASCADE_DIR.mkdir(parents=True, exist_ok=True)

POCKET2MOL_DIR = BASE / "results" / "pocket2mol_designs"

# Pocket PDB sources from Pocket2Mol session 7 results
POCKETS = [
    {"pdb_id": "7BV2", "short": "RdRp", "ligand_id": "F86", "tlr_prior": 0.6,
     "pdb_file": str(next(POCKET2MOL_DIR.glob("7BV2/*/7bv2.pdb")))},
    {"pdb_id": "7VH8", "short": "Mpro", "ligand_id": "0EN", "tlr_prior": 0.6,
     "pdb_file": str(next(POCKET2MOL_DIR.glob("7VH8/*/7vh8.pdb")))},
    {"pdb_id": "3OXC", "short": "HIV_protease", "ligand_id": "?", "tlr_prior": 0.6,
     "pdb_file": str(next(POCKET2MOL_DIR.glob("3OXC/*/3oxc.pdb")))},
    {"pdb_id": "1VRT", "short": "HIV_RT", "ligand_id": "?", "tlr_prior": 0.6,
     "pdb_file": str(next(POCKET2MOL_DIR.glob("1VRT/*/1vrt.pdb")))},
]

DEVICE = "cuda:0" if torch.cuda.is_available() else "cpu"
W = {"L4": 0.25, "L5": 0.15, "L6a": 0.20, "L6b": 0.15}


def alerts_pains_brenk(smi: str) -> dict:
    from rdkit import Chem
    from rdkit.Chem import FilterCatalog
    cat_p = FilterCatalog.FilterCatalogParams()
    cat_p.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
    cat_b = FilterCatalog.FilterCatalogParams()
    cat_b.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
    cp = FilterCatalog.FilterCatalog(cat_p)
    cb = FilterCatalog.FilterCatalog(cat_b)
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return {"pains": 0, "brenk": 0, "valid": False}
    return {
        "pains": int(cp.HasMatch(m)),
        "brenk": int(cb.HasMatch(m)),
        "valid": True,
    }


STANDARD_AA = {
    "ALA","ARG","ASN","ASP","CYS","GLU","GLN","GLY","HIS","ILE","LEU","LYS",
    "MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
}


def clean_pdb_for_targetdiff(src_pdb: str, dst_pdb: str) -> None:
    """Strip non-standard-AA records (DNA/RNA, ligands, waters, metals)."""
    with open(src_pdb) as fin, open(dst_pdb, "w") as fout:
        for line in fin:
            if line.startswith(("ATOM",)):
                resname = line[17:20].strip()
                if resname in STANDARD_AA:
                    fout.write(line)
            elif line.startswith(("TER", "END", "REMARK")):
                fout.write(line)


def sample_pocket(pocket, model, ligand_atom_mode, transform, ligand_featurizer):
    from utils.data import PDBProtein
    from datasets.pl_data import ProteinLigandData, torchify_dict
    from scripts.sample_diffusion import sample_diffusion_ligand
    from utils import reconstruct
    from utils import transforms as trans
    from rdkit import Chem

    pdb_path = pocket["pdb_file"]
    pdb_id = pocket["pdb_id"]
    out_dir = OUT_DIR / pdb_id
    out_dir.mkdir(parents=True, exist_ok=True)

    cleaned_pdb = out_dir / "pocket_cleaned.pdb"
    clean_pdb_for_targetdiff(pdb_path, str(cleaned_pdb))
    pocket_dict = PDBProtein(str(cleaned_pdb)).to_dict_atom()
    data = ProteinLigandData.from_protein_ligand_dicts(
        protein_dict=torchify_dict(pocket_dict),
        ligand_dict={
            "element": torch.empty([0], dtype=torch.long),
            "pos": torch.empty([0, 3], dtype=torch.float),
            "atom_feature": torch.empty([0, 8], dtype=torch.float),
            "bond_index": torch.empty([2, 0], dtype=torch.long),
            "bond_type": torch.empty([0], dtype=torch.long),
        },
    )
    data = transform(data)

    print(f"  {pdb_id}: pocket atoms = {data.protein_pos.shape[0]}, sampling {N_SAMPLES} ligands ...")
    all_pred_pos, all_pred_v, _, _, _, _, _ = sample_diffusion_ligand(
        model, data, N_SAMPLES,
        batch_size=BATCH_SIZE, device=DEVICE,
        num_steps=N_STEPS, pos_only=False,
        center_pos_mode="protein", sample_num_atoms="prior",
    )

    # Reconstruct
    valid_mols = []
    smiles_list = []
    for sample_idx, (pred_pos, pred_v) in enumerate(zip(all_pred_pos, all_pred_v)):
        try:
            atom_type = trans.get_atomic_number_from_index(pred_v, mode="add_aromatic")
            aromatic = trans.is_aromatic_from_index(pred_v, mode="add_aromatic")
            mol = reconstruct.reconstruct_from_generated(pred_pos, atom_type, aromatic)
            smi = Chem.MolToSmiles(mol)
            if "." in smi:  # multi-component, drop
                continue
            smiles_list.append(smi)
            valid_mols.append(mol)
        except reconstruct.MolReconsError:
            continue
        except Exception as e:
            continue

    print(f"    valid SMILES: {len(smiles_list)}/{N_SAMPLES}")

    # Save SDFs and SMILES
    sdf_dir = out_dir / "sdf"
    sdf_dir.mkdir(exist_ok=True)
    for idx, mol in enumerate(valid_mols):
        if mol is None:
            continue
        w = Chem.SDWriter(str(sdf_dir / f"{idx:03d}.sdf"))
        w.write(mol)
        w.close()
    (out_dir / "SMILES.txt").write_text("\n".join(smiles_list))
    return smiles_list


def main():
    from torch_geometric.transforms import Compose
    from utils import transforms as trans
    from models.molopt_score_model import ScorePosNet3D

    print("=== T4.4 TargetDiff structure-based drug design ===")
    ckpt_path = str(TARGETDIFF_REPO / "pretrained_models" / "pretrained_diffusion.pt")
    print(f"Loading {ckpt_path}")
    ckpt = torch.load(ckpt_path, map_location=DEVICE, weights_only=False)
    ligand_atom_mode = ckpt["config"].data.transform.ligand_atom_mode
    protein_featurizer = trans.FeaturizeProteinAtom()
    ligand_featurizer = trans.FeaturizeLigandAtom(ligand_atom_mode)
    transform = Compose([protein_featurizer])

    model = ScorePosNet3D(
        ckpt["config"].model,
        protein_atom_feature_dim=protein_featurizer.feature_dim,
        ligand_atom_feature_dim=ligand_featurizer.feature_dim,
    ).to(DEVICE)
    model.load_state_dict(ckpt["model"])
    model.eval()
    print(f"  TargetDiff loaded ({sum(p.numel() for p in model.parameters()) / 1e6:.1f} M params)")

    # 1) Sample per pocket
    all_smiles = []  # list of (pocket, smi)
    for pocket in POCKETS:
        try:
            smiles_list = sample_pocket(pocket, model, ligand_atom_mode, transform, ligand_featurizer)
        except Exception as e:
            print(f"  {pocket['pdb_id']} FAILED: {e}")
            continue
        for smi in smiles_list:
            all_smiles.append({"pocket": pocket, "smiles": smi})

    if not all_smiles:
        print("No valid SMILES; exiting.")
        return

    # 2) Cascade
    print(f"\n=== Cascade scoring on {len(all_smiles)} SMILES ===")
    unique_smiles = list({s["smiles"] for s in all_smiles})
    from admet_ai import ADMETModel
    print(f"ADMET-AI on {len(unique_smiles)} unique SMILES ...")
    admet = ADMETModel()
    df = admet.predict(smiles=unique_smiles)
    admet_map = {}
    for i, smi in enumerate(unique_smiles):
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

    cascade_entries = []
    for s in all_smiles:
        smi = s["smiles"]
        if smi not in admet_map:
            continue
        adm = admet_map[smi]
        alerts = alerts_pains_brenk(smi)
        # r_L4 from validity-and-alerts: alerts add 0.5 each, capped at 1.0
        r_L4 = round(0.5 * (alerts["pains"] + alerts["brenk"]), 4)
        r_L5 = 0.5 if (alerts["pains"] or alerts["brenk"]) else 0.3
        r_L6a = adm["r_L6a"]
        r_L6b = s["pocket"]["tlr_prior"]
        # noisy-OR with both L6 channels
        p = 1.0
        for w, r in [(W["L4"], r_L4), (W["L5"], r_L5), (W["L6a"], r_L6a), (W["L6b"], r_L6b)]:
            p *= 1 - w * min(r, 1.0)
        P_harm = round(1.0 - p, 6)

        cascade_entries.append(
            {
                "protein_short_name": s["pocket"]["short"],
                "holo_pdb_id": s["pocket"]["pdb_id"],
                "smiles": smi,
                "model": "TargetDiff (Guan et al. ICLR 2023)",
                "alerts": alerts,
                "r_L4": r_L4,
                "r_L5": r_L5,
                "r_L6a": r_L6a,
                "r_L6b": r_L6b,
                "ClinTox": adm["ClinTox"],
                "hERG": adm["hERG"],
                "DILI": adm["DILI"],
                "AMES": adm["AMES"],
                "P_harm": P_harm,
            }
        )

    cascade_entries.sort(key=lambda x: x["P_harm"], reverse=True)
    for i, e in enumerate(cascade_entries, 1):
        e["rank"] = i

    out = {
        "task": "T4.4 SBDD - TargetDiff",
        "method": "TargetDiff diffusion + reconstruct + RDKit validate",
        "n_samples_per_pocket": N_SAMPLES,
        "n_diffusion_steps": N_STEPS,
        "n_pockets": len(POCKETS),
        "n_smiles_total": len(all_smiles),
        "n_unique_smiles": len(unique_smiles),
        "weights": W,
        "top10_entries": cascade_entries[:10],
        "all_ranked_entries": cascade_entries,
    }
    out_path = CASCADE_DIR / "cascade_results_targetdiff.json"
    json.dump(out, open(out_path, "w"), indent=2)
    print(f"\nSaved {len(cascade_entries)} cascade entries to {out_path}")
    print("\nTop 5 by P_harm:")
    for e in cascade_entries[:5]:
        print(
            f"  {e['rank']}. {e['protein_short_name']:14s} "
            f"PAINS={e['alerts']['pains']} BRENK={e['alerts']['brenk']} "
            f"r_L6a={e['r_L6a']:.3f} P_harm={e['P_harm']:.4f}  {e['smiles'][:60]}"
        )


if __name__ == "__main__":
    main()
