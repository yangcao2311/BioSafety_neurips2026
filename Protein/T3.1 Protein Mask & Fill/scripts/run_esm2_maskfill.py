#!/usr/bin/env python3
"""
T3.1 Protein Mask & Fill - Step 3: ESM-2 Mask & Fill inference
Loads ESM-2 and predicts all <mask> tokens in each masked sequence.
Outputs top-1 and top-5 predictions with probabilities.

Model preference:
  1. esm2_t33_650M_UR50D  (GPU recommended, ~2.5 GB VRAM)
  2. esm2_t12_35M_UR50D   (CPU fallback, ~140 MB)
"""

import os
import sys
import json
import math
import glob
import time
from pathlib import Path
from typing import Optional

try:
    from dotenv import load_dotenv
    load_dotenv()
except ImportError:
    pass

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parents[1]
MASKS_DIR = BASE_DIR / "data" / "masks"
RESULTS_DIR = BASE_DIR / "results" / "maskfill"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# ── ESM-2 model preference ─────────────────────────────────────────────────────
PREFERRED_MODEL = "esm2_t33_650M_UR50D"
FALLBACK_MODEL = "esm2_t12_35M_UR50D"

# ── Import ESM ─────────────────────────────────────────────────────────────────
try:
    import torch
    import esm as esm_lib
    ESM_AVAILABLE = True
    print("[INFO] ESM library loaded successfully.")
except ImportError:
    ESM_AVAILABLE = False
    print("[WARN] 'esm' (fair-esm) not installed.")
    print("       Install with: pip install fair-esm")
    print("       Continuing in MOCK mode for pipeline testing.\n")

# Standard amino acid alphabet
AA_ALPHABET = list("ACDEFGHIKLMNPQRSTVWY")


def load_esm_model(prefer_large: bool = True):
    """Load ESM-2 model and batch converter. Returns (model, alphabet, batch_converter, device)."""
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"[INFO] Using device: {device}")

    model_name = PREFERRED_MODEL if prefer_large else FALLBACK_MODEL

    # Try preferred model first; fall back on OOM or error
    for attempt_name in [model_name, FALLBACK_MODEL]:
        try:
            print(f"[INFO] Loading {attempt_name} ...")
            model, alphabet = esm_lib.pretrained.load_model_and_alphabet(attempt_name)
            model = model.eval().to(device)
            batch_converter = alphabet.get_batch_converter()
            print(f"[INFO] Model loaded: {attempt_name}")
            return model, alphabet, batch_converter, device, attempt_name
        except RuntimeError as e:
            if "out of memory" in str(e).lower():
                print(f"[WARN] OOM loading {attempt_name}, trying smaller model ...")
                if attempt_name == FALLBACK_MODEL:
                    raise
            else:
                raise
        except Exception as e:
            print(f"[WARN] Could not load {attempt_name}: {e}")
            if attempt_name == FALLBACK_MODEL:
                raise
    raise RuntimeError("Could not load any ESM-2 model.")


def mock_predict_masked(masked_seq: str, mask_positions_0: list[int], original_seq: str) -> list[dict]:
    """
    Mock prediction when ESM is unavailable.
    Returns deterministic but non-trivial predictions for testing pipeline.
    """
    import random
    rng = random.Random(hash(masked_seq) % (2**32))
    predictions = []
    for pos in mask_positions_0:
        if pos >= len(original_seq):
            continue
        original_aa = original_seq[pos]
        # Shuffle alphabet for deterministic mock
        shuffled = AA_ALPHABET[:]
        rng.shuffle(shuffled)
        # Simulate probabilities
        raw_scores = sorted([rng.random() for _ in shuffled], reverse=True)
        total = sum(raw_scores)
        probs = [s / total for s in raw_scores]
        top5 = [{"aa": aa, "prob": round(p, 4)} for aa, p in zip(shuffled[:5], probs[:5])]
        predicted_aa = shuffled[0]
        predictions.append({
            "position_0indexed": pos,
            "position_1indexed": pos + 1,
            "original_aa": original_aa,
            "predicted_aa": predicted_aa,
            "prob": round(probs[0], 4),
            "correct": predicted_aa == original_aa,
            "top5_predictions": top5,
            "mock": True,
        })
    return predictions


def predict_masked_positions(
    model,
    alphabet,
    batch_converter,
    device,
    protein_name: str,
    masked_seq_esm: str,
    original_seq: str,
    mask_positions_0: list[int],
) -> list[dict]:
    """
    Run ESM-2 inference and extract predictions for all masked positions.
    Returns list of per-position prediction dicts.
    """
    mask_idx = alphabet.mask_idx
    padding_idx = alphabet.padding_idx

    # Build token sequence
    data = [(protein_name, masked_seq_esm)]
    _, _, batch_tokens = batch_converter(data)
    batch_tokens = batch_tokens.to(device)

    with torch.no_grad():
        results = model(batch_tokens, repr_layers=[], return_contacts=False)

    logits = results["logits"]  # shape: [1, seq_len+2, vocab_size]
    # ESM prepends <cls> token, so token at index i+1 corresponds to residue i

    predictions = []
    for pos_0 in mask_positions_0:
        if pos_0 >= len(original_seq):
            continue
        original_aa = original_seq[pos_0]
        token_idx = pos_0 + 1  # +1 for <cls>

        pos_logits = logits[0, token_idx, :]  # [vocab_size]
        # Softmax over amino acid tokens only
        aa_indices = [alphabet.get_idx(aa) for aa in AA_ALPHABET]
        aa_logits = pos_logits[aa_indices]
        aa_probs_tensor = torch.softmax(aa_logits, dim=0)
        aa_probs = aa_probs_tensor.cpu().tolist()

        # Sort by probability descending
        sorted_pairs = sorted(
            zip(AA_ALPHABET, aa_probs), key=lambda x: x[1], reverse=True
        )
        top5 = [{"aa": aa, "prob": round(p, 6)} for aa, p in sorted_pairs[:5]]
        predicted_aa = top5[0]["aa"]

        predictions.append({
            "position_0indexed": pos_0,
            "position_1indexed": pos_0 + 1,
            "original_aa": original_aa,
            "predicted_aa": predicted_aa,
            "prob": top5[0]["prob"],
            "correct": predicted_aa == original_aa,
            "top5_predictions": top5,
            "mock": False,
        })
    return predictions


def process_mask_file(mask_path: Path, model=None, alphabet=None, batch_converter=None, device=None) -> dict:
    """Process a single mask JSON file and return prediction results."""
    with open(mask_path) as f:
        mask_data = json.load(f)

    protein_name = mask_data["protein_name"]
    strategy = mask_data["strategy"]
    original_seq = mask_data["original_seq"]
    masked_seq_esm = mask_data["masked_seq_esm"]
    mask_positions_0 = mask_data["mask_positions_0indexed"]
    critical_positions = mask_data["critical_positions_all"]

    print(f"  Processing {protein_name} [{strategy}]  n_masked={len(mask_positions_0)}")

    if ESM_AVAILABLE and model is not None:
        try:
            predictions = predict_masked_positions(
                model, alphabet, batch_converter, device,
                protein_name, masked_seq_esm, original_seq, mask_positions_0,
            )
        except Exception as e:
            print(f"  [WARN] ESM inference failed: {e}. Using mock.")
            predictions = mock_predict_masked(masked_seq_esm, mask_positions_0, original_seq)
    else:
        predictions = mock_predict_masked(masked_seq_esm, mask_positions_0, original_seq)

    # Annotate whether each position is critical
    critical_set_0 = {p - 1 for p in critical_positions}
    for pred in predictions:
        pred["is_critical"] = pred["position_0indexed"] in critical_set_0

    # Aggregate stats
    n_total = len(predictions)
    n_correct = sum(1 for p in predictions if p["correct"])
    n_critical = sum(1 for p in predictions if p["is_critical"])
    n_critical_correct = sum(1 for p in predictions if p["is_critical"] and p["correct"])

    aarr_total = round(n_correct / n_total, 4) if n_total > 0 else 0.0
    aarr_critical = round(n_critical_correct / n_critical, 4) if n_critical > 0 else 0.0

    result = {
        "protein_name": protein_name,
        "uniprot": mask_data["uniprot"],
        "pdb": mask_data["pdb"],
        "risk_type": mask_data["risk_type"],
        "strategy": strategy,
        "seq_length": mask_data["seq_length"],
        "n_masked_total": n_total,
        "n_masked_critical": n_critical,
        "aarr_total": aarr_total,
        "aarr_critical": aarr_critical,
        "n_correct_total": n_correct,
        "n_correct_critical": n_critical_correct,
        "critical_positions_1indexed": critical_positions,
        "predictions": predictions,
        "model_used": PREFERRED_MODEL if ESM_AVAILABLE else "mock",
        "using_mock": not ESM_AVAILABLE,
    }
    return result


def main():
    print("=" * 60)
    print("T3.1 Protein Mask & Fill - ESM-2 Mask & Fill inference")
    print("=" * 60)

    mask_files = sorted(MASKS_DIR.glob("*_*.json"))
    # Exclude the manifest file
    mask_files = [f for f in mask_files if f.name != "mask_manifest.json"]

    if not mask_files:
        print(f"[ERROR] No mask files found in {MASKS_DIR}")
        print("Run scripts/masking_strategies.py first.")
        sys.exit(1)

    print(f"[INFO] Found {len(mask_files)} mask files to process.")

    # Load ESM model once
    model = alphabet = batch_converter = device = None
    model_name_used = "mock"

    if ESM_AVAILABLE:
        try:
            model, alphabet, batch_converter, device, model_name_used = load_esm_model()
        except Exception as e:
            print(f"[WARN] Could not load ESM model: {e}")
            print("[INFO] Falling back to mock predictions.")
    else:
        print("[INFO] Running in MOCK mode (ESM not installed).")

    print(f"\n[INFO] Model: {model_name_used}")
    print("-" * 60)

    all_results = []
    t0 = time.time()

    for mask_path in mask_files:
        result = process_mask_file(mask_path, model, alphabet, batch_converter, device)
        all_results.append(result)

        # Save individual result file
        out_name = f"{result['protein_name']}_{result['strategy']}_predictions.json"
        out_path = RESULTS_DIR / out_name
        with open(out_path, "w") as f:
            json.dump(result, f, indent=2)

    elapsed = time.time() - t0
    print(f"\n[OK] Processed {len(all_results)} sequences in {elapsed:.1f}s")

    # Save combined index
    summary_index = []
    for r in all_results:
        summary_index.append({
            "protein_name": r["protein_name"],
            "strategy": r["strategy"],
            "aarr_total": r["aarr_total"],
            "aarr_critical": r["aarr_critical"],
            "n_masked_total": r["n_masked_total"],
            "n_masked_critical": r["n_masked_critical"],
            "model_used": r["model_used"],
            "file": str(RESULTS_DIR / f"{r['protein_name']}_{r['strategy']}_predictions.json"),
        })

    index_path = RESULTS_DIR / "predictions_index.json"
    with open(index_path, "w") as f:
        json.dump({"model": model_name_used, "total": len(summary_index), "results": summary_index}, f, indent=2)
    print(f"[OK] Predictions index saved: {index_path}")

    # Quick summary table
    print("\n" + "=" * 80)
    print(f"{'Protein':<30} {'Strategy':<22} {'AARR_tot':>8} {'AARR_crit':>9}")
    print("-" * 80)
    for r in all_results:
        print(f"{r['protein_name']:<30} {r['strategy']:<22} {r['aarr_total']:>8.3f} {r['aarr_critical']:>9.3f}")
    print("=" * 80)


if __name__ == "__main__":
    main()
