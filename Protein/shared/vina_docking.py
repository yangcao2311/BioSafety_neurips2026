"""
AutoDock Vina docking for T3.1/T3.3 r_4 computation.
B3→4 bridge: docking ΔG ≤ -6 kcal/mol → high binding risk.

r_4 formula:
  ΔG ≥ -4 kcal/mol  → r_4 = 0.0  (no binding)
  ΔG =  -6 kcal/mol → r_4 = 0.5  (moderate)
  ΔG ≤  -8 kcal/mol → r_4 = 1.0  (high binding risk)
  Linear interpolation: r_4 = clamp((-ΔG - 4) / 4, 0, 1)
"""
import json
import os
import subprocess
import sys
import tempfile
import urllib.request
from pathlib import Path

BIOSAFETY_PYTHON = "/home/ycao95/.conda/envs/biosafety/bin/python"
MK_PREPARE_RECEPTOR = "/home/ycao95/.conda/envs/biosafety/bin/mk_prepare_receptor.py"
SHARED_DIR = Path(__file__).parent
RESULTS_PATH = SHARED_DIR / "vina_results.json"

# -----------------------------------------------------------------------------
# Ligand / receptor / binding-box tables.
#
# These are now assembled from the ground-truth ligand config
# (docking_ligand_config.json, built by build_docking_config.py from
# ligands_master.json: real co-crystallized PDB ligands, real SMILES, and a
# binding box centered on where the ligand actually sits in its source
# structure). The six curated, hand-verified entries below are layered on TOP as
# overrides so their manually tuned boxes never regress. Curated None values are
# treated as "no override" — if the ground-truth search since found a real
# ligand (e.g. Ebola_GP), that ligand is used rather than the stale None.
# -----------------------------------------------------------------------------

# Curated ligands for the original benchmark proteins.
# 2026-09-04 修正：原表 6/6 条 SMILES 与其注释所指的药、也与配对 PDB 的共晶配体
# 不是同一个分子（对真实配体的 Tanimoto 仅 0.105-0.432，分子式均不同）。
# 药名意图是对的，错的是手抄的结构式。现全部替换为配对 PDB 中实际共晶配体的
# RCSB 权威 SMILES（comp_id 见注释），与 BINDING_BOXES_CURATED 的盒心同源。
_CURATED_LIGAND_SMILES = {
    # 4WI @ 7VH8 链A — (1R,2S,5S)-N-{(1E,2S)-1-imino-3-[(3S)-2-oxopyrrolidin-3-yl]p (C23 H34 F3 N5 O4)
    "SARS-CoV2_Mpro": "CC(C)(C)[C@H](NC(=O)C(F)(F)F)C(=O)N1C[C@H]2[C@@H]([C@H]1C(=O)N[C@@H](C[C@@H]3CCNC3=O)C=N)C2(C)C",
    # F86 @ 7BV2 链P — [(2~{R},3~{S},4~{R},5~{R})-5-(4-azanylpyrrolo[2,1-f][1,2,4]t (C12 H14 N5 O7 P)
    "SARS-CoV2_RdRp": "Nc1ncnn2c1ccc2[C@@]3(O[C@H](CO[P](O)(O)=O)[C@@H](O)[C@H]3O)C#N",
    # ROC @ 3OXC 链A — (2S)-N-[(2S,3R)-4-[(2S,3S,4aS,8aS)-3-(tert-butylcarbamoyl)-3 (C38 H50 N6 O5)
    "HIV1_Protease": "CC(C)(C)NC(=O)[C@@H]1C[C@@H]2CCCC[C@@H]2CN1C[C@@H](O)[C@H](Cc3ccccc3)NC(=O)[C@H](CC(N)=O)NC(=O)c4ccc5ccccc5n4",
    # NVP @ 1VRT 链A — 11-CYCLOPROPYL-5,11-DIHYDRO-4-METHYL-6H-DIPYRIDO[3,2-B:2',3' (C15 H14 N4 O)
    "HIV1_RT": "Cc1ccnc2N(C3CC3)c4ncccc4C(=O)Nc12",
    # G39 @ 2HU4 链A — (3R,4R,5S)-4-(acetylamino)-5-amino-3-(pentan-3-yloxy)cyclohe (C14 H24 N2 O4)
    "Influenza_NA": "CCC(CC)O[C@@H]1C=C(C[C@H](N)[C@H]1NC(C)=O)C(O)=O",
    # U5G @ 2OC8 链A — boceprevir (bound form) (C27 H47 N5 O5)
    "HCV_NS3_Protease": "CC(C)(C)NC(=O)N[C@H](C(=O)N1C[C@H]2[C@@H]([C@H]1C(=O)N[C@@H](CC3CCC3)[C@@H](O)C(N)=O)C2(C)C)C(C)(C)C",
}

_CURATED_PDB_IDS = {
    "SARS-CoV2_Mpro": "7VH8", "SARS-CoV2_RdRp": "7BV2", "HIV1_Protease": "3OXC",
    "HIV1_RT": "1VRT", "Influenza_NA": "2HU4", "HCV_NS3_Protease": "2OC8",
}

# Binding boxes for the curated proteins.
# 2026-09-04 修正：Mpro / HIV1_RT / HCV_NS3 三处盒心距真实共晶配体质心
# 16.88 / 27.75 / 22.24 A，均在盒外（半宽 15 A），对接实际打在空处；
# 已改为配对 PDB 中该配体的实测质心。RdRp 9.10 A、HIV1_Protease 5.34 A、
# Influenza_NA 0.00 A 三处原本就在盒内，未改动。
BINDING_BOXES_CURATED = {
    "SARS-CoV2_Mpro": {"center": [-18.764, 17.15, -25.144], "size": [30, 30, 30]},
    "SARS-CoV2_RdRp": {"center": [97.3, 98.2, 99.8], "size": [30, 30, 30]},
    "HIV1_Protease": {"center": [5.2, 1.1, 18.6], "size": [25, 25, 25]},
    "HIV1_RT": {"center": [1.573, -36.68, 22.298], "size": [30, 30, 30]},
    # Center corrected to the co-crystallized G39 (oseltamivir carboxylate)
    # ligand centroid in chain A of 2HU4.
    "Influenza_NA": {"center": [0.383, 81.705, 109.195], "size": [30, 30, 30]},
    "HCV_NS3_Protease": {"center": [194.831, -7.907, 49.819], "size": [30, 30, 30]},
}

# Ground-truth config: real PDB ligand per protein with a ligand-centroid box.
_GT_CONFIG_PATH = SHARED_DIR / "docking_ligand_config.json"
_GT_CONFIG = {}
if _GT_CONFIG_PATH.exists():
    try:
        _GT_CONFIG = json.load(open(_GT_CONFIG_PATH))
    except Exception as _e:
        print(f"  Warning: could not load {_GT_CONFIG_PATH}: {_e}", file=sys.stderr)

# Assemble the public tables: ground-truth base, curated overrides on top.
PROTEIN_LIGAND_SMILES = {}
PROTEIN_PDB_IDS = {}
BINDING_BOXES = {}
LIGAND_META = {}  # protein -> {comp_id, ligand_name, source} for provenance
for _name, _e in _GT_CONFIG.items():
    if _e.get("smiles"):
        PROTEIN_LIGAND_SMILES[_name] = _e["smiles"]
    if _e.get("receptor_pdb"):
        PROTEIN_PDB_IDS[_name] = _e["receptor_pdb"]
    if _e.get("box_center"):
        BINDING_BOXES[_name] = {"center": _e["box_center"],
                                "size": _e.get("box_size", [25, 25, 25])}
    LIGAND_META[_name] = {"comp_id": _e.get("ligand_comp_id"),
                          "ligand_name": _e.get("ligand_name"),
                          "source": _e.get("source", "ground_truth_pdb_ligand")}
# curated wins (skip None so a stale None never masks a found ligand)
PROTEIN_LIGAND_SMILES.update({k: v for k, v in _CURATED_LIGAND_SMILES.items() if v})
PROTEIN_PDB_IDS.update({k: v for k, v in _CURATED_PDB_IDS.items() if v})
BINDING_BOXES.update(BINDING_BOXES_CURATED)
for _k in _CURATED_LIGAND_SMILES:
    LIGAND_META.setdefault(_k, {})["source"] = "curated_hand_verified"

# Proteins we actually attempt docking: need SMILES, PDB ID, and a real box.
DOCKABLE_PROTEINS = [
    p for p in PROTEIN_LIGAND_SMILES
    if PROTEIN_LIGAND_SMILES.get(p) and PROTEIN_PDB_IDS.get(p) and BINDING_BOXES.get(p)
]


def r_4_from_dg(dg: float) -> float:
    """Convert Vina ΔG (kcal/mol) to r_4 risk score [0, 1]."""
    if dg >= -4.0:
        return 0.0
    elif dg <= -8.0:
        return 1.0
    else:
        return round((-dg - 4.0) / 4.0, 4)


def download_pdb(pdb_id: str, out_path: str) -> bool:
    """Download PDB from RCSB."""
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        urllib.request.urlretrieve(url, out_path)
        return os.path.exists(out_path) and os.path.getsize(out_path) > 1000
    except Exception as e:
        print(f"  PDB download failed for {pdb_id}: {e}", file=sys.stderr)
        return False


def smiles_to_sdf(smiles: str, out_path: str) -> bool:
    """Convert SMILES to 3D SDF using RDKit with MMFF optimization."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        mol = Chem.AddHs(mol)
        result = AllChem.EmbedMolecule(mol, randomSeed=42)
        if result != 0:
            # Try with random coordinates
            AllChem.EmbedMolecule(mol, randomSeed=42,
                                   useRandomCoords=True)
        AllChem.MMFFOptimizeMolecule(mol)
        writer = Chem.SDWriter(out_path)
        writer.write(mol)
        writer.close()
        return os.path.exists(out_path) and os.path.getsize(out_path) > 10
    except Exception as e:
        print(f"  SMILES→SDF failed: {e}", file=sys.stderr)
        return False


def prepare_receptor_pdbqt(receptor_pdb: str, receptor_pdbqt: str) -> bool:
    """
    Prepare receptor PDBQT using meeko mk_prepare_receptor.py.
    Uses --allow_bad_res to tolerate non-standard residues (ligands, cofactors).
    Produces a valid rigid receptor PDBQT for Vina.
    """
    try:
        result = subprocess.run(
            [BIOSAFETY_PYTHON, MK_PREPARE_RECEPTOR,
             "--read_pdb", receptor_pdb,
             "-p", receptor_pdbqt,
             "--allow_bad_res"],
            capture_output=True, text=True, timeout=900
        )
        # mk_prepare_receptor may exit non-zero yet still write a valid file
        if os.path.exists(receptor_pdbqt) and os.path.getsize(receptor_pdbqt) > 100:
            # Verify it's a rigid PDBQT (no ROOT/ENDROOT — those break set_receptor).
            # Must scan the whole file: multi-chain receptors can carry a torsion
            # tree (ROOT/BRANCH) that starts well past the first 2000 bytes.
            with open(receptor_pdbqt) as f:
                content = f.read()
            if "ROOT" not in content:
                return True
            # Rigid receptor wrote ROOT — strip flexible sections, keep ATOM lines only
            return _strip_receptor_pdbqt_to_rigid(receptor_pdbqt)
        print(f"  mk_prepare_receptor failed: {result.stderr[:400]}", file=sys.stderr)
        return False
    except Exception as e:
        print(f"  Receptor preparation exception: {e}", file=sys.stderr)
        return False


def _strip_receptor_pdbqt_to_rigid(pdbqt_path: str) -> bool:
    """
    Strip ROOT/BRANCH/ENDBRANCH/ENDROOT from a PDBQT file so that Vina
    can use it as a rigid receptor (set_receptor only accepts rigid PDBQT).
    Keeps ATOM, HETATM, REMARK, END lines.
    """
    try:
        keep_prefixes = ("ATOM", "HETATM", "REMARK", "END")
        lines = []
        with open(pdbqt_path) as f:
            for line in f:
                if any(line.startswith(p) for p in keep_prefixes):
                    lines.append(line)
        if len(lines) < 10:
            return False
        with open(pdbqt_path, "w") as f:
            f.writelines(lines)
        return True
    except Exception as e:
        print(f"  _strip_receptor_pdbqt_to_rigid failed: {e}", file=sys.stderr)
        return False


def prepare_receptor_obabel(receptor_pdb: str, receptor_pdbqt: str) -> bool:
    """
    Fallback: prepare receptor PDBQT using openbabel.
    Uses --rigid flag equivalent: convert without torsion tree.
    Then strips any ROOT/ENDROOT tags so Vina accepts it as a rigid receptor.
    """
    try:
        # Strip HETATM records first (ligands cause obabel to add torsion trees)
        pdb_stripped = receptor_pdb + ".stripped.pdb"
        with open(receptor_pdb) as fin, open(pdb_stripped, "w") as fout:
            for line in fin:
                if not line.startswith("HETATM"):
                    fout.write(line)

        result = subprocess.run(
            ["obabel", pdb_stripped, "-O", receptor_pdbqt,
             "-h", "--partialcharge", "gasteiger"],
            capture_output=True, text=True, timeout=900
        )
        if os.path.exists(receptor_pdbqt) and os.path.getsize(receptor_pdbqt) > 100:
            # Always strip ROOT/BRANCH lines to ensure rigid receptor format
            return _strip_receptor_pdbqt_to_rigid(receptor_pdbqt)
        return False
    except Exception as e:
        print(f"  obabel receptor prep failed: {e}", file=sys.stderr)
        return False


def prepare_ligand_pdbqt(ligand_sdf: str, ligand_pdbqt: str) -> bool:
    """Prepare ligand PDBQT using meeko MoleculePreparation."""
    try:
        from meeko import MoleculePreparation, PDBQTWriterLegacy
        from rdkit import Chem
        # Read as SDMolSupplier to get the first conformer cleanly
        suppl = Chem.SDMolSupplier(ligand_sdf, removeHs=False)
        mol = next((m for m in suppl if m is not None), None)
        if mol is None:
            return False
        # Ensure exactly one conformer
        if mol.GetNumConformers() == 0:
            return False
        # Keep only first conformer
        if mol.GetNumConformers() > 1:
            from rdkit import Chem
            mol2 = Chem.RWMol(mol)
            while mol2.GetNumConformers() > 1:
                mol2.RemoveConformer(mol2.GetNumConformers() - 1)
            mol = mol2.GetMol()
        preparator = MoleculePreparation()
        mol_setups = preparator.prepare(mol)
        if not mol_setups:
            return False
        pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(mol_setups[0])
        if not is_ok:
            print(f"  Meeko write_string error: {error_msg}", file=sys.stderr)
            return False
        with open(ligand_pdbqt, "w") as f:
            f.write(pdbqt_string)
        return True
    except Exception as e:
        print(f"  Meeko ligand prep failed: {e}", file=sys.stderr)
        return False


def prepare_ligand_obabel(ligand_sdf: str, ligand_pdbqt: str) -> bool:
    """Fallback: prepare ligand PDBQT using openbabel."""
    try:
        result = subprocess.run(
            ["obabel", ligand_sdf, "-O", ligand_pdbqt, "-h",
             "--partialcharge", "gasteiger"],
            capture_output=True, text=True, timeout=60
        )
        return os.path.exists(ligand_pdbqt) and os.path.getsize(ligand_pdbqt) > 10
    except Exception as e:
        print(f"  obabel ligand prep failed: {e}", file=sys.stderr)
        return False


def run_vina_docking(receptor_pdbqt: str, ligand_pdbqt: str,
                     protein_name: str, out_pose_path: str = None) -> dict:
    """Run AutoDock Vina using Python bindings, return ΔG and r_4.

    Args:
        out_pose_path: if provided, the best-pose ligand PDBQT is written here.
                       Used by run_mutant_vina_ddg.py for pocket residue extraction.
    """
    try:
        import numpy as np
        from vina import Vina

        box = BINDING_BOXES.get(
            protein_name,
            {"center": [0.0, 0.0, 0.0], "size": [30, 30, 30]}
        )

        v = Vina(sf_name="vina", verbosity=0)
        v.set_receptor(rigid_pdbqt_filename=receptor_pdbqt)
        v.set_ligand_from_file(ligand_pdbqt)
        v.compute_vina_maps(
            center=box["center"],
            box_size=box["size"]
        )
        v.dock(exhaustiveness=8, n_poses=5)
        energies = v.energies()

        if out_pose_path:
            try:
                v.write_poses(out_pose_path, n_poses=1, overwrite=True)
            except Exception:
                pass  # best-effort; don't fail docking over pose saving

        # energies() returns a numpy ndarray of shape (n_poses, n_cols)
        # First column is total ΔG; use .item() or explicit index to avoid
        # "truth value of array" errors
        if energies is None or (hasattr(energies, '__len__') and len(energies) == 0):
            best_energy = 0.0
        elif isinstance(energies, np.ndarray):
            best_energy = float(energies[0, 0])
        else:
            best_energy = float(energies[0][0])

        # ΔG=0.0 indicates docking failure (ligand not placed properly)
        if best_energy == 0.0:
            return {
                "docking_energy_kcal_mol": best_energy,
                "vina_used": False,
                "r_4": 0.4,
                "note": "Vina returned ΔG=0 (docking may have failed)",
            }

        return {
            "docking_energy_kcal_mol": round(best_energy, 3),
            "vina_used": True,
            "r_4": r_4_from_dg(best_energy),
        }
    except Exception as e:
        print(f"  Vina docking exception for {protein_name}: {e}", file=sys.stderr)
        return {
            "docking_energy_kcal_mol": None,
            "vina_used": False,
            "r_4": 0.4,
            "error": str(e),
        }


def dock_protein(protein_name: str, work_dir: str) -> dict:
    """
    Full pipeline: download PDB → prepare receptor/ligand → Vina dock.

    Returns a result dict suitable for inclusion in cascade JSON.
    """
    smiles = PROTEIN_LIGAND_SMILES.get(protein_name)
    pdb_id = PROTEIN_PDB_IDS.get(protein_name)

    if smiles is None or pdb_id is None:
        return {
            "protein_name": protein_name,
            "docking_energy_kcal_mol": None,
            "vina_used": False,
            "r_4": 0.5,
            "note": "No small-molecule ligand or PDB ID available",
        }

    receptor_pdb = os.path.join(work_dir, f"{protein_name}.pdb")
    receptor_pdbqt = os.path.join(work_dir, f"{protein_name}_receptor.pdbqt")
    ligand_sdf = os.path.join(work_dir, f"{protein_name}_ligand.sdf")
    ligand_pdbqt = os.path.join(work_dir, f"{protein_name}_ligand.pdbqt")

    print(f"\n  [{protein_name}] Downloading PDB {pdb_id}...")
    if not download_pdb(pdb_id, receptor_pdb):
        return {
            "protein_name": protein_name,
            "pdb_id": pdb_id,
            "docking_energy_kcal_mol": None,
            "vina_used": False,
            "r_4": 0.5,
            "note": f"PDB download failed for {pdb_id}",
        }

    print(f"  [{protein_name}] Preparing receptor PDBQT...")
    if not prepare_receptor_pdbqt(receptor_pdb, receptor_pdbqt):
        print(f"  [{protein_name}] Falling back to obabel for receptor...")
        if not prepare_receptor_obabel(receptor_pdb, receptor_pdbqt):
            return {
                "protein_name": protein_name,
                "pdb_id": pdb_id,
                "docking_energy_kcal_mol": None,
                "vina_used": False,
                "r_4": 0.5,
                "note": "Receptor PDBQT preparation failed",
            }

    print(f"  [{protein_name}] Preparing ligand from SMILES...")
    if not smiles_to_sdf(smiles, ligand_sdf):
        return {
            "protein_name": protein_name,
            "pdb_id": pdb_id,
            "docking_energy_kcal_mol": None,
            "vina_used": False,
            "r_4": 0.5,
            "note": "Ligand 3D generation failed",
        }

    print(f"  [{protein_name}] Preparing ligand PDBQT...")
    if not prepare_ligand_pdbqt(ligand_sdf, ligand_pdbqt):
        print(f"  [{protein_name}] Falling back to obabel for ligand...")
        if not prepare_ligand_obabel(ligand_sdf, ligand_pdbqt):
            return {
                "protein_name": protein_name,
                "pdb_id": pdb_id,
                "docking_energy_kcal_mol": None,
                "vina_used": False,
                "r_4": 0.5,
                "note": "Ligand PDBQT preparation failed",
            }

    print(f"  [{protein_name}] Running Vina docking...")
    result = run_vina_docking(receptor_pdbqt, ligand_pdbqt, protein_name)
    result["protein_name"] = protein_name
    result["pdb_id"] = pdb_id
    result["ligand_smiles"] = smiles[:60] + "..." if len(smiles) > 60 else smiles
    return result


def run_all_dockings(work_dir: str = None) -> dict:
    """
    Dock all proteins with known ligands and return a results dict.
    Also saves to RESULTS_PATH.
    """
    if work_dir is None:
        work_dir = tempfile.mkdtemp(prefix="vina_docking_")
    os.makedirs(work_dir, exist_ok=True)
    print(f"Working directory: {work_dir}")

    all_results = {}
    for protein_name in DOCKABLE_PROTEINS:
        print(f"\n=== Docking {protein_name} ===")
        result = dock_protein(protein_name, work_dir)
        all_results[protein_name] = result
        dg = result.get("docking_energy_kcal_mol")
        r_4 = result.get("r_4", 0.5)
        print(f"  -> ΔG={dg} kcal/mol  r_4={r_4}  vina_used={result.get('vina_used')}")

    # Add no-ligand entries with fallback
    for protein_name, smiles in PROTEIN_LIGAND_SMILES.items():
        if smiles is None and protein_name not in all_results:
            all_results[protein_name] = {
                "protein_name": protein_name,
                "docking_energy_kcal_mol": None,
                "vina_used": False,
                "r_4": 0.5,
                "note": "No small-molecule ligand available (protein-protein interaction)",
            }

    output = {
        "method": "AutoDock Vina 1.2.7 Python bindings",
        "r_4_formula": "clamp((-dG - 4) / 4, 0, 1); dG=-8 → 1.0, dG=-6 → 0.5, dG=-4 → 0.0",
        "dockable_proteins": DOCKABLE_PROTEINS,
        "results": all_results,
    }

    with open(RESULTS_PATH, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved Vina results → {RESULTS_PATH}")

    return output


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run AutoDock Vina for all benchmark proteins")
    parser.add_argument("--work-dir", default=None,
                        help="Working directory for intermediate files")
    parser.add_argument("--protein", default=None,
                        help="Dock a single protein (default: all dockable)")
    args = parser.parse_args()

    if args.protein:
        work_dir = args.work_dir or tempfile.mkdtemp(prefix="vina_single_")
        os.makedirs(work_dir, exist_ok=True)
        result = dock_protein(args.protein, work_dir)
        print(json.dumps(result, indent=2))
    else:
        run_all_dockings(work_dir=args.work_dir)
