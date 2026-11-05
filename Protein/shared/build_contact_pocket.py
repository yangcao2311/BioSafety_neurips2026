#!/usr/bin/env python3
"""
build_contact_pocket.py -- Rebuild binding-pocket residue sets from real
co-crystal structures, fixing two defects in the original
`vina_r_l4_upgrade.py::parse_pdb_pocket_residues()`:

  DEFECT 1 (silent cross-chain loss)
      The original loop never filtered by chain, yet stored the result under
      the name `pocket_residues_chain_A`.  Downstream `compute_pocket_mut_rate`
      indexes strictly by chain-A residue numbers, so every pocket residue
      living on another chain was silently dropped.  Measured impact:
      18/142 proteins affected, 83/935 residues (8.9%) lost.  Ebola_GP is the
      worst case -- 3 of its 4 pocket residues sit on GP2 (chain B), leaving
      only 1 residue actually driving r_4.

  DEFECT 2 (centroid proxy instead of real contacts)
      The original criterion was "CA within 8 A of the ligand *centroid*".
      For an elongated ligand this simultaneously misses genuine contacts at
      both ends and admits non-contacting residues near the middle.  The
      standard criterion -- and the one used here -- is "any protein heavy
      atom within CONTACT_CUTOFF_A of any ligand heavy atom".

Output is written to a NEW field/file; nothing already on disk is overwritten,
so previously reported r_4 numbers stay reproducible until we decide to
backfill.

No network access and no docking rerun: every structure is read from the
existing local cache at A_L3_Protein/shared/pdb_cache/.
"""

from __future__ import annotations

import json
import math
from collections import defaultdict
from pathlib import Path

SHARED = Path(__file__).resolve().parent
PDB_CACHE = SHARED / "pdb_cache"
LIGAND_CONFIG = SHARED / "docking_ligand_config.json"
VINA_BASELINES = SHARED / "vina_wt_baselines.json"
OUT_JSON = SHARED / "contact_pockets.json"
OUT_DIFF = SHARED / "contact_pocket_diff.md"

# Heavy-atom contact distance. 4.5 A is the common convention for
# protein-ligand contact residues (covers hydrogen bonds ~2.8-3.2 A and
# van der Waals contacts ~3.5-4.2 A with a small margin).
CONTACT_CUTOFF_A = 4.5

# Solvent / cryo-protectant / buffer components that are present in almost
# every structure and are never the biological ligand of interest.
SOLVENT_COMP_IDS = {
    "HOH", "DOD", "GOL", "EDO", "PEG", "PG4", "DMS", "MPD", "ACT", "SO4",
    "PO4", "CL", "NA", "MG", "CA", "ZN", "K", "IOD", "BR", "NO3", "TRS",
}


# ---------------------------------------------------------------------------
# PDB parsing
# ---------------------------------------------------------------------------

def parse_structure(pdb_path: Path):
    """Return (protein_atoms, hetatm_by_comp).

    protein_atoms : list of (x, y, z, chain, resnum, resname)  -- heavy atoms only
    hetatm_by_comp: dict comp_id -> list of (x, y, z, chain)   -- heavy atoms only
    """
    protein_atoms = []
    hetatm_by_comp = defaultdict(list)

    with open(pdb_path) as fh:
        for line in fh:
            rec = line[:6]
            if rec not in ("ATOM  ", "HETATM"):
                continue
            element = line[76:78].strip().upper()
            atom_name = line[12:16].strip()
            # Skip hydrogens (element column when present, else name heuristic).
            if element == "H" or (not element and atom_name.startswith("H")):
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            chain = line[21]
            resname = line[17:20].strip()

            if rec == "ATOM  ":
                try:
                    resnum = int(line[22:26])
                except ValueError:
                    continue
                protein_atoms.append((x, y, z, chain, resnum, resname))
            else:
                if resname in SOLVENT_COMP_IDS:
                    continue
                hetatm_by_comp[resname].append((x, y, z, chain))

    return protein_atoms, hetatm_by_comp


def contact_residues(protein_atoms, ligand_atoms, cutoff=CONTACT_CUTOFF_A):
    """Residues with >=1 heavy atom within `cutoff` of any ligand heavy atom.

    Returns dict chain -> sorted list of residue numbers, plus a flat
    (chain, resnum, resname) set for reporting.
    """
    if not ligand_atoms:
        return {}, set()

    cutoff_sq = cutoff * cutoff

    # Bucket ligand atoms into a coarse grid so the scan stays linear in
    # practice instead of |protein| x |ligand|.
    cell = cutoff
    grid = defaultdict(list)
    for lx, ly, lz, _ in ligand_atoms:
        key = (int(math.floor(lx / cell)), int(math.floor(ly / cell)), int(math.floor(lz / cell)))
        grid[key].append((lx, ly, lz))

    hits = set()
    for px, py, pz, chain, resnum, resname in protein_atoms:
        cx = int(math.floor(px / cell))
        cy = int(math.floor(py / cell))
        cz = int(math.floor(pz / cell))
        found = False
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for lx, ly, lz in grid.get((cx + dx, cy + dy, cz + dz), ()):
                        if (px - lx) ** 2 + (py - ly) ** 2 + (pz - lz) ** 2 <= cutoff_sq:
                            found = True
                            break
                    if found:
                        break
                if found:
                    break
            if found:
                break
        if found:
            hits.add((chain, resnum, resname))

    by_chain = defaultdict(list)
    for chain, resnum, _ in hits:
        by_chain[chain].append(resnum)
    return {c: sorted(set(v)) for c, v in sorted(by_chain.items())}, hits


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ligand_cfg = json.loads(LIGAND_CONFIG.read_text())
    baselines = json.loads(VINA_BASELINES.read_text())

    results = {}
    skipped = {}

    for name, cfg in sorted(ligand_cfg.items()):
        pdb_id = cfg.get("receptor_pdb")
        comp_id = cfg.get("ligand_comp_id")
        if not pdb_id or not comp_id:
            skipped[name] = "no receptor_pdb or ligand_comp_id"
            continue

        pdb_path = PDB_CACHE / f"{pdb_id}.pdb"
        if not pdb_path.exists():
            skipped[name] = f"pdb not cached: {pdb_id}"
            continue

        protein_atoms, hetatm_by_comp = parse_structure(pdb_path)
        ligand_atoms = hetatm_by_comp.get(comp_id, [])
        if not ligand_atoms:
            skipped[name] = f"ligand {comp_id} not found as HETATM in {pdb_id}"
            continue

        by_chain, hits = contact_residues(protein_atoms, ligand_atoms)
        n_res = sum(len(v) for v in by_chain.values())

        results[name] = {
            "receptor_pdb": pdb_id,
            "ligand_comp_id": comp_id,
            "ligand_heavy_atoms": len(ligand_atoms),
            "ligand_chains": sorted({c for *_, c in ligand_atoms}),
            "contact_pocket_by_chain": by_chain,
            "n_contact_residues": n_res,
            "method": f"cocrystal_heavy_atom_contact_{CONTACT_CUTOFF_A}A",
            "source": "real co-crystal structure (local pdb_cache), no network",
        }

    # ---- diff against the old centroid-based pockets -----------------------
    rows = []
    for name, res in sorted(results.items()):
        base = baselines.get(name)
        if not base or "pocket_residues_chain_A" not in base:
            continue
        old = set(base["pocket_residues_chain_A"])
        chain_a_nums = set(base.get("wt_seq_resnums_A", []))
        new_all = {(c, r) for c, nums in res["contact_pocket_by_chain"].items() for r in nums}
        new_nums = {r for _, r in new_all}

        # Which old residues would survive the downstream chain-A indexing?
        old_usable = old & chain_a_nums
        # NOTE: `old_covered_by_new` is a QA statistic for this migration only --
        # "how much of the old pocket list is contained in the new one".  It is
        # NOT a pipeline metric, is never scored, and has nothing to do with the
        # deprecated apo-pocket overlap metric.  It exists so we can see how far
        # the corrected pockets moved before deciding whether to backfill r_4.
        rows.append({
            "name": name,
            "pdb": res["receptor_pdb"],
            "ligand": res["ligand_comp_id"],
            "n_old": len(old),
            "n_old_usable_chainA": len(old_usable),
            "n_old_lost": len(old - chain_a_nums),
            "n_new": len(new_nums),
            "n_new_chains": len(res["contact_pocket_by_chain"]),
            "n_old_also_in_new": len(old & new_nums),
            "old_covered_by_new_pct": round(100 * len(old & new_nums) / len(old), 1) if old else None,
        })

    OUT_JSON.write_text(json.dumps({
        "contact_cutoff_A": CONTACT_CUTOFF_A,
        "n_proteins": len(results),
        "n_skipped": len(skipped),
        "note": (
            "Rebuilt from real co-crystal structures. Fixes two defects in "
            "vina_r_l4_upgrade.py::parse_pdb_pocket_residues(): (1) cross-chain "
            "residues silently dropped, (2) CA-to-ligand-centroid 8A proxy "
            "replaced by heavy-atom-to-heavy-atom contact. Written to a new "
            "file; nothing existing is overwritten."
        ),
        "pockets": results,
        "skipped": skipped,
        "diff_vs_old_centroid_pockets": rows,
    }, indent=2))

    # ---- human-readable diff ----------------------------------------------
    lost = [r for r in rows if r["n_old_lost"] > 0]
    lines = [
        "# Contact-pocket rebuild -- diff vs old centroid pockets",
        "",
        f"- Contact criterion: any protein heavy atom within **{CONTACT_CUTOFF_A} A** of any ligand heavy atom",
        f"- Proteins rebuilt: **{len(results)}**  (skipped {len(skipped)})",
        f"- Proteins whose OLD pocket silently lost cross-chain residues: **{len(lost)}**",
        "",
        "> The last column is a one-off QA statistic for this migration, not a",
        "> pipeline metric and unrelated to the deprecated apo-pocket overlap score.",
        "",
        "| protein | PDB | ligand | old n | old usable (chain A) | old lost | new n | new chains | old residues also in new | % of old kept |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in sorted(rows, key=lambda r: (-r["n_old_lost"], r["name"])):
        lines.append(
            f"| {r['name']} | {r['pdb']} | {r['ligand']} | {r['n_old']} | "
            f"{r['n_old_usable_chainA']} | {r['n_old_lost']} | {r['n_new']} | "
            f"{r['n_new_chains']} | {r['n_old_also_in_new']} | {r['old_covered_by_new_pct']} |"
        )
    if skipped:
        lines += ["", "## Skipped", ""]
        for k, v in sorted(skipped.items()):
            lines.append(f"- `{k}`: {v}")
    OUT_DIFF.write_text("\n".join(lines) + "\n")

    print(f"rebuilt : {len(results)} proteins")
    print(f"skipped : {len(skipped)}")
    print(f"old pockets losing cross-chain residues: {len(lost)}")
    print(f"-> {OUT_JSON}")
    print(f"-> {OUT_DIFF}")


if __name__ == "__main__":
    main()
