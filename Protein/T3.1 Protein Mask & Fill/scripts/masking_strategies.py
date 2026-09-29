#!/usr/bin/env python3
"""
T3.1 Protein Mask & Fill - Step 2: Generate masked sequences
Applies 4 masking strategies per protein and saves JSON mask files.

Strategies:
  1. active_site        - mask all critical_positions (up to 8 residues)
  2. binding_interface  - mask critical_positions[0:4]
  3. immune_epitope     - mask 9-mer window around critical_positions[0]
  4. random_15pct       - random 15% of sequence positions
"""

import os
import sys
import json
import random
from pathlib import Path

try:
    from dotenv import load_dotenv
    load_dotenv()
except ImportError:
    pass

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data"
PROTEINS_DIR = DATA_DIR / "proteins"
MASKS_DIR = DATA_DIR / "masks"
MASKS_DIR.mkdir(parents=True, exist_ok=True)

# ESM-2 mask token
ESM_MASK_TOKEN = "<mask>"
VIZ_MASK_CHAR = "X"  # For visualization / FASTA display

RANDOM_SEED = 42


def load_protein_config() -> list[dict]:
    config_path = PROTEINS_DIR / "protein_config.json"
    if not config_path.exists():
        raise FileNotFoundError(
            f"protein_config.json not found at {config_path}.\n"
            "Run scripts/download_proteins.py first."
        )
    with open(config_path) as f:
        return json.load(f)


def apply_mask(sequence: str, positions_0indexed: list[int]) -> tuple[str, str]:
    """
    Apply mask tokens to specified 0-indexed positions.
    Returns (esm_masked_seq, viz_masked_seq).
    ESM masked seq uses space-separated tokens; positions get <mask>.
    For ESM we return a plain string where masked residues are replaced by <mask>.
    """
    seq_list = list(sequence)
    # ESM variant: replace chars
    esm_seq = list(sequence)
    viz_seq = list(sequence)
    for pos in positions_0indexed:
        if 0 <= pos < len(sequence):
            esm_seq[pos] = ESM_MASK_TOKEN
            viz_seq[pos] = VIZ_MASK_CHAR
    # ESM masked sequence: join — masked positions become the literal string <mask>
    # We track them as a token list, joined by empty string for non-masked segments
    # but we need a single string for the ESM tokenizer.
    # ESM tokenizer accepts raw amino acid strings; we replace the char with <mask>
    # We keep it as a concatenated string (ESM tokenizer will parse <mask> as a token).
    esm_masked = "".join(esm_seq)
    viz_masked = "".join(viz_seq)
    return esm_masked, viz_masked


def strategy_active_site(sequence: str, critical_positions: list[int]) -> dict:
    """Mask all critical_positions (1-indexed in config → convert to 0-indexed)."""
    positions_0 = [p - 1 for p in critical_positions if 0 < p <= len(sequence)]
    esm_masked, viz_masked = apply_mask(sequence, positions_0)
    return {
        "strategy": "active_site",
        "mask_positions_1indexed": [p for p in critical_positions if 0 < p <= len(sequence)],
        "mask_positions_0indexed": positions_0,
        "n_masked": len(positions_0),
        "original_seq": sequence,
        "masked_seq_esm": esm_masked,
        "masked_seq_viz": viz_masked,
        "description": "All critical catalytic/binding residues masked",
    }


def strategy_binding_interface(sequence: str, critical_positions: list[int]) -> dict:
    """Mask first 4 critical positions (binding interface subset)."""
    subset = critical_positions[:4]
    positions_0 = [p - 1 for p in subset if 0 < p <= len(sequence)]
    esm_masked, viz_masked = apply_mask(sequence, positions_0)
    return {
        "strategy": "binding_interface",
        "mask_positions_1indexed": [p for p in subset if 0 < p <= len(sequence)],
        "mask_positions_0indexed": positions_0,
        "n_masked": len(positions_0),
        "original_seq": sequence,
        "masked_seq_esm": esm_masked,
        "masked_seq_viz": viz_masked,
        "description": "First 4 binding interface residues masked",
    }


def strategy_immune_epitope(sequence: str, critical_positions: list[int]) -> dict:
    """Mask a 9-mer window centered on (or starting at) critical_positions[0]."""
    if not critical_positions:
        # Fallback: mask first 9 residues
        anchor = 1
    else:
        anchor = critical_positions[0]

    # 9-mer window: 4 before anchor, anchor, 4 after (1-indexed)
    half = 4
    start_1 = max(1, anchor - half)
    end_1 = min(len(sequence), start_1 + 8)  # inclusive, 9 positions
    if end_1 - start_1 < 8:
        start_1 = max(1, end_1 - 8)

    window_1indexed = list(range(start_1, end_1 + 1))
    positions_0 = [p - 1 for p in window_1indexed]
    esm_masked, viz_masked = apply_mask(sequence, positions_0)
    return {
        "strategy": "immune_epitope",
        "mask_positions_1indexed": window_1indexed,
        "mask_positions_0indexed": positions_0,
        "n_masked": len(positions_0),
        "anchor_position_1indexed": anchor,
        "original_seq": sequence,
        "masked_seq_esm": esm_masked,
        "masked_seq_viz": viz_masked,
        "description": f"9-mer window [{start_1}:{end_1}] around critical position {anchor} for T-cell epitope simulation",
    }


def strategy_random_15pct(sequence: str, critical_positions: list[int], seed: int = RANDOM_SEED) -> dict:
    """Randomly mask 15% of sequence positions."""
    rng = random.Random(seed)
    n_mask = max(1, int(len(sequence) * 0.15))
    all_positions = list(range(len(sequence)))
    positions_0 = sorted(rng.sample(all_positions, n_mask))
    positions_1 = [p + 1 for p in positions_0]
    esm_masked, viz_masked = apply_mask(sequence, positions_0)

    # Which critical positions were accidentally masked?
    critical_set_0 = {p - 1 for p in critical_positions}
    critical_masked = [p + 1 for p in positions_0 if p in critical_set_0]

    return {
        "strategy": "random_15pct",
        "mask_positions_1indexed": positions_1,
        "mask_positions_0indexed": positions_0,
        "n_masked": len(positions_0),
        "target_fraction": 0.15,
        "actual_fraction": round(len(positions_0) / len(sequence), 4),
        "critical_positions_incidentally_masked": critical_masked,
        "original_seq": sequence,
        "masked_seq_esm": esm_masked,
        "masked_seq_viz": viz_masked,
        "description": "Random 15% of sequence positions masked",
        "random_seed": seed,
    }


STRATEGY_FUNCTIONS = {
    "active_site": strategy_active_site,
    "binding_interface": strategy_binding_interface,
    "immune_epitope": strategy_immune_epitope,
    "random_15pct": strategy_random_15pct,
}


def main():
    print("=" * 60)
    print("T3.1 Protein Mask & Fill - Generating masked sequences")
    print("=" * 60)

    proteins = load_protein_config()
    total_masks = 0

    for protein in proteins:
        name = protein["name"]
        sequence = protein["sequence"]
        critical_positions = protein["critical_positions"]
        seq_len = len(sequence)

        print(f"\n[{name}]  length={seq_len}  critical={critical_positions}")

        for strategy_name, strategy_fn in STRATEGY_FUNCTIONS.items():
            result = strategy_fn(sequence, critical_positions)

            # Add protein metadata
            result["protein_name"] = name
            result["uniprot"] = protein["uniprot"]
            result["pdb"] = protein["pdb"]
            result["risk_type"] = protein["risk_type"]
            result["seq_length"] = seq_len
            result["critical_positions_all"] = critical_positions

            out_path = MASKS_DIR / f"{name}_{strategy_name}.json"
            with open(out_path, "w") as f:
                json.dump(result, f, indent=2)

            n = result["n_masked"]
            pct = round(n / seq_len * 100, 1)
            print(f"  [{strategy_name:<20}]  masked={n:>3} ({pct:>5.1f}%)  -> {out_path.name}")
            total_masks += 1

    print(f"\n[OK] Generated {total_masks} mask files in {MASKS_DIR}")

    # Write a manifest
    manifest = []
    for f in sorted(MASKS_DIR.glob("*.json")):
        manifest.append(str(f))

    manifest_path = MASKS_DIR / "mask_manifest.json"
    with open(manifest_path, "w") as f:
        json.dump({"total": len(manifest), "files": manifest}, f, indent=2)
    print(f"[OK] Manifest saved: {manifest_path}")


if __name__ == "__main__":
    main()
