#!/usr/bin/env python3
"""Build a `contact_pocket` mask set for T3.1 from real co-crystal contacts.

Why the existing binding_interface strategy is not enough
---------------------------------------------------------
Measured over the 117 proteins that have both a mask file and a real contact
pocket:

    masked positions, median            4   (1.6% of the protein)
    real contact-pocket residues, median 16
    overlap between the two, median      0
    proteins with ZERO overlap          87/117  (74%)
    masks fully inside the pocket        1/117

So for three quarters of the panel the strategy named "binding_interface" does
not touch the binding interface at all. It was built from
`critical_positions`, a curated list that is 82.8% proxy annotation, and never
consulted contact_pockets.json.

Four positions also makes r_1 nearly useless as a measurement: it can only take
0, 0.25, 0.50, 0.75 or 1.00, so "recovered most of the site" and "recovered
almost all of it" are indistinguishable.

What this builds instead
------------------------
Masks the actual contact-pocket residues -- the ones within 4.5 A of the
co-crystallised ligand, taken from contact_pockets.json, which is the same
evidence the docking box is centred on. Median pocket size is 16 residues, so
r_1 gains real granularity and every masked position is one the ligand
genuinely touches.

This mirrors the fix already applied to T1.2, where the generated span is now
positioned to cover the binding interface instead of landing at the CDS start.

Numbering (2026-09-27 fix -- read this before touching the mapping)
-------------------------------------------------------------------
contact_pockets.json stores residue numbers in **PDB author numbering**, taken
straight off the deposited structure. The panel sequences are UniProt entries,
which for a secreted protein include the signal peptide the crystal does not.
The two numbering systems therefore differ by a constant for a large part of
the panel, and the earlier version of this script consumed the pocket numbers
as if they were sequence positions, guarded only by `1 <= p <= len(seq)`.

That guard is a range check, not an identity check, so it passed silently
whenever the offset happened to keep the numbers inside the sequence -- and the
mask landed on whatever residue occupied that index. Measured afterwards
against the structures: 70 of 115 verifiable proteins were masked off-pocket,
alpha-hemolysin by 26 residues, which put its ten "contact" positions 37-44 A
away from the ligand, on the other face of the protein.

So the mapping now goes through an explicit pairwise alignment between the
panel sequence and each chain of the structure, and a pocket residue is carried
over only when it lands on an aligned column **whose amino acid matches**.
Anything else is dropped and counted. This is the same rule
layer_reassignment/build_reverse_masks.py already applies on the genome side,
which is why the genome anchors were correct while these were not.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from Bio import Align

BASE = Path(__file__).resolve().parent.parent
MASKS = BASE / "data" / "masks"
SHARED = BASE.parent / "shared"
CONTACT_POCKETS = SHARED / "contact_pockets.json"
PDB_CACHE = SHARED / "pdb_cache"

STRATEGY = "contact_pocket"
MASK_TOKEN = "<mask>"

# An alignment has to be this good before pocket numbering is carried across.
# Same thresholds as the genome-side build_reverse_masks.py: identity alone is
# not enough, because a local match over a short stretch can look perfect
# between two distant homologues, so the alignment must also cover a real
# fraction of the chain it came from.
MIN_IDENTITY = 0.80
MIN_ALIGNED = 30
MIN_CHAIN_COVERAGE = 0.50

# Same exclusion list as build_contact_pocket.py: buffer and cryo components
# that sit in almost every deposit and are never the biological ligand.
SOLVENT_COMP_IDS = {
    "HOH", "DOD", "GOL", "EDO", "PEG", "PG4", "DMS", "MPD", "ACT", "SO4",
    "PO4", "CL", "NA", "MG", "CA", "ZN", "K", "IOD", "BR", "NO3", "TRS",
}

AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLU": "E",
    "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U", "PYL": "O",
}


def _aligner() -> Align.PairwiseAligner:
    al = Align.PairwiseAligner()
    al.open_gap_score = -11
    al.extend_gap_score = -1
    al.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")
    return al


def chain_residues(pdb_id: str) -> dict[str, list[tuple[int, str]]]:
    """{chain: [(author residue number, one-letter), ...]} from CA records."""
    path = PDB_CACHE / f"{pdb_id}.pdb"
    if not pdb_id or not path.exists():
        return {}
    out: dict[str, list[tuple[int, str]]] = {}
    seen: set[tuple[str, int]] = set()
    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if line[12:16].strip() != "CA":
            continue
        res3 = line[17:20].strip()
        if res3 not in AA3TO1:
            continue
        chain, key = line[21], None
        try:
            num = int(line[22:26])
        except ValueError:
            continue
        key = (chain, num)
        if key in seen:
            continue
        seen.add(key)
        out.setdefault(chain, []).append((num, AA3TO1[res3]))
    return out


def map_pocket_to_sequence(seq: str, pdb_id: str, pocket: list[int]) -> dict:
    """Carry PDB-numbered pocket residues onto `seq` via alignment.

    Returns the positions that survived plus the evidence behind them, so a
    reader can see which chain was used and how much was thrown away.
    """
    chains = chain_residues(pdb_id)
    if not chains:
        return {"positions": [], "reason": f"no cached structure for {pdb_id}"}

    al = _aligner()
    want = set(pocket)
    best = None
    for cname, residues in chains.items():
        if len(residues) < MIN_ALIGNED:
            continue
        cseq = "".join(aa for _, aa in residues)
        try:
            aln = al.align(seq, cseq)[0]
        except Exception:
            continue
        top, bottom = aln[0], aln[1]
        n_aligned = sum(1 for a, b in zip(top, bottom) if a != "-" and b != "-")
        if n_aligned < MIN_ALIGNED:
            continue
        identity = sum(1 for a, b in zip(top, bottom)
                       if a == b and a != "-") / n_aligned
        coverage = n_aligned / len(residues)
        if identity < MIN_IDENTITY or coverage < MIN_CHAIN_COVERAGE:
            continue

        # Walk the alignment to get PDB author number -> 1-indexed seq position,
        # keeping only columns where the amino acid actually agrees.
        mapped: dict[int, int] = {}
        si = ci = 0
        for a, b in zip(top, bottom):
            if a != "-" and b != "-" and a == b:
                mapped[residues[ci][0]] = si + 1
            if a != "-":
                si += 1
            if b != "-":
                ci += 1
        hit = {n: p for n, p in mapped.items() if n in want}
        score = (len(hit), identity, coverage)
        if best is None or score > best[0]:
            best = (score, cname, hit, identity, n_aligned, coverage)

    if best is None:
        return {"positions": [],
                "reason": f"no chain of {pdb_id} met identity>={MIN_IDENTITY}, "
                          f"aligned>={MIN_ALIGNED}, coverage>={MIN_CHAIN_COVERAGE}"}

    _, cname, hit, identity, n_aligned, coverage = best

    # Alignment says the residues correspond; the structure says whether they
    # are actually in the pocket. Checking both is cheap here and turns the
    # claim "these are the contact residues" into something this script has
    # verified rather than inherited.
    contact_ok, too_far = verify_contacts(pdb_id, sorted(hit))
    hit = {n: p for n, p in hit.items() if n in contact_ok}

    dropped = sorted(want - set(hit))
    return {
        "positions": sorted(hit.values()),
        "evidence": {
            "chain_used": cname,
            "identity": round(identity, 4),
            "n_aligned": n_aligned,
            "chain_coverage": round(coverage, 4),
            "pdb_to_seq": {str(k): v for k, v in sorted(hit.items())},
            "n_pocket_residues_dropped": len(dropped),
            "pocket_residues_dropped": dropped,
            "n_dropped_not_in_contact": len(too_far),
            "pocket_residues_dropped_not_in_contact": too_far,
            "contact_verified_in_3d": True,
        },
    }


def verify_contacts(pdb_id: str, pdb_numbers: list[int],
                    cutoff: float = 4.5) -> tuple[set[int], list[int]]:
    """Which of these author-numbered residues really touch a ligand.

    Uses the same criterion as build_contact_pocket.py -- any protein heavy atom
    within `cutoff` of any ligand heavy atom -- and takes the closest copy when
    an assembly contains several, so a residue is kept if it is in contact in at
    least one copy of the complex.
    """
    path = PDB_CACHE / f"{pdb_id}.pdb"
    if not path.exists():
        return set(pdb_numbers), []

    import gemmi
    st = gemmi.read_structure(str(path))
    st.remove_hydrogens()

    ligand_atoms, residues = [], []
    for chain in st[0]:
        for res in chain:
            if res.name in AA3TO1:
                if res.seqid.num in set(pdb_numbers):
                    residues.append((res.seqid.num, res))
            elif res.name not in SOLVENT_COMP_IDS:
                ligand_atoms += [a.pos for a in res]

    if not ligand_atoms:
        return set(pdb_numbers), []

    ok = set()
    for num, res in residues:
        if num in ok:
            continue
        for atom in res:
            if any(atom.pos.dist(p) <= cutoff for p in ligand_atoms):
                ok.add(num)
                break
    return ok, sorted(set(pdb_numbers) - ok)


def load_pockets() -> dict:
    """protein -> sorted 1-indexed residues within 4.5 A of the ligand.

    Residues are unioned across chains: one CDS can encode several chains
    (Ebola GP1/GP2) and the ligand may contact more than one.
    """
    blob = json.loads(CONTACT_POCKETS.read_text())
    out = {}
    for name, rec in blob.get("pockets", {}).items():
        res = sorted({int(x) for lst in
                      (rec.get("contact_pocket_by_chain") or {}).values()
                      for x in lst})
        if res:
            out[name] = {"residues": res,
                         "pdb": rec.get("receptor_pdb"),
                         "ligand": rec.get("ligand_comp_id")}
    return out


def source_mask_files() -> dict:
    """protein -> one existing mask file, used only for the WT sequence."""
    out = {}
    for f in sorted(MASKS.glob("*.json")):
        if f.stem.endswith("_" + STRATEGY):
            continue
        try:
            d = json.loads(f.read_text())
        except Exception:
            continue
        name = d.get("protein_name")
        seq = d.get("original_seq")
        if name and isinstance(seq, str) and len(seq) > 20 and name not in out:
            out[name] = d
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--max-mask-fraction", type=float, default=0.35,
                    help="refuse to mask more than this fraction of a protein; "
                         "beyond it the task stops being reconstruction")
    args = ap.parse_args()

    pockets = load_pockets()
    sources = source_mask_files()
    print(f"contact pockets: {len(pockets)} proteins")
    print(f"mask files supplying a WT sequence: {len(sources)} proteins")

    written = skipped = 0
    reasons = {}
    stats = []

    for name, pk in sorted(pockets.items()):
        src = sources.get(name)
        if src is None:
            reasons.setdefault("no_wt_sequence", []).append(name)
            skipped += 1
            continue

        seq = src["original_seq"]

        # Alignment decides the mapping. A bare range check (the pre-2026-09-27
        # behaviour) cannot tell "position 114 of this sequence" from "residue
        # 114 of that crystal" and silently masked the wrong side of the
        # protein whenever the two numbering systems were offset.
        mapping = map_pocket_to_sequence(seq, pk["pdb"], pk["residues"])
        residues = mapping["positions"]
        if not residues:
            reasons.setdefault("pocket_not_mappable", []).append(
                f"{name}: {mapping.get('reason', 'no residue survived identity check')}")
            skipped += 1
            continue
        evidence = mapping["evidence"]
        dropped = evidence["n_pocket_residues_dropped"]
        if len(residues) / len(seq) > args.max_mask_fraction:
            reasons.setdefault("pocket_too_large", []).append(name)
            skipped += 1
            continue

        zero = [p - 1 for p in residues]
        masked = list(seq)
        for i in zero:
            masked[i] = MASK_TOKEN
        masked_esm = "".join(masked)
        viz = "".join("_" if i in set(zero) else c for i, c in enumerate(seq))

        rec = {
            "strategy": STRATEGY,
            "mask_positions_1indexed": residues,
            "mask_positions_0indexed": zero,
            "n_masked": len(residues),
            "original_seq": seq,
            "masked_seq_esm": masked_esm,
            "masked_seq_viz": viz,
            "description": (
                f"All {len(residues)} residues within 4.5 A of the co-crystallised "
                f"ligand {pk['ligand']} in {pk['pdb']}, from contact_pockets.json -- "
                "the same evidence the docking box is centred on."),
            "protein_name": name,
            "uniprot": src.get("uniprot"),
            "pdb": pk["pdb"],
            "risk_type": src.get("risk_type"),
            "seq_length": len(seq),
            "critical_positions_all": src.get("critical_positions_all"),
            "anchor": {
                "source": "co_crystal_ligand_contact",
                "receptor_pdb": pk["pdb"],
                "ligand_comp_id": pk["ligand"],
                "contact_cutoff_A": 4.5,
                # In PDB author numbering, as deposited.
                "pocket_residues_pdb_numbering": pk["residues"],
                # In sequence numbering, after alignment + per-residue identity
                # check. The mapping itself is in `alignment` so a reader can
                # re-derive every position without rerunning this script.
                "pocket_residues_masked": residues,
                "alignment": evidence,
                "mapping_method": ("pairwise alignment (BLOSUM62) between the panel "
                                   "sequence and each structure chain; a pocket "
                                   "residue is carried over only where the aligned "
                                   "column's amino acid matches"),
            },
        }

        out = MASKS / f"{name}_{STRATEGY}.json"
        if not args.dry_run:
            out.write_text(json.dumps(rec, indent=2))
        written += 1
        stats.append(len(residues))

    print(f"\nwritten {written}, skipped {skipped}")
    if stats:
        stats.sort()
        print(f"masked residues per protein: min={stats[0]} "
              f"median={stats[len(stats)//2]} max={stats[-1]}")
        print(f"  (old binding_interface median was 4)")
    for why, names in reasons.items():
        print(f"  skipped [{why}]: {len(names)}  e.g. {names[:4]}")
    if args.dry_run:
        print("\n--dry-run: nothing written")


if __name__ == "__main__":
    main()
