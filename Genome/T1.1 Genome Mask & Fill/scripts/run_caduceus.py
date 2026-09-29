#!/usr/bin/env python3
"""
T1.1 Caduceus (Mamba-based bidirectional DNA LM) sequence filling.

Caduceus (Schiff et al. 2024) is a Mamba-based foundation model for DNA with
reverse-complement equivariance. The "ph" variant uses palindromic head
(predicts same sequence on forward and RC strands).

Model: kuleshov-group/caduceus-ph_seqlen-1k_d_model-118_n_layer-4_lr-8e-3

Note: Caduceus requires `mamba_ssm` (CUDA-compiled). If unavailable, this script
will fall back to "mamba_ssm not installed" warning and skip.

Usage: python scripts/run_caduceus.py
"""
import json
import os
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import torch

BASE = Path(__file__).resolve().parent.parent
DATA_FILE = BASE / "data" / "masked_genomes.json"
RESULTS_DIR = BASE / "results" / "maskfill"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

DEVICE = "cuda:0" if torch.cuda.is_available() else "cpu"
MODEL_ID = "kuleshov-group/caduceus-ph_seqlen-131k_d_model-256_n_layer-16"

# Caduceus vocab: 0=CLS, 1=SEP, 2=BOS, 3=MASK, 4=PAD, 5=RESERVED, 6=UNK,
#                 7=A, 8=C, 9=G, 10=T, 11=N. Logits have 16 columns; we filter to ACGT.
ACGT_IDS = {7: "A", 8: "C", 9: "G", 10: "T"}


def load_model():
    """Load Caduceus tokenizer + model. Returns None if mamba_ssm missing."""
    os.environ.setdefault("HF_HOME", "/home/xliu316/.cache/huggingface")
    try:
        from transformers import AutoTokenizer, AutoModelForMaskedLM
        tok = AutoTokenizer.from_pretrained(MODEL_ID, trust_remote_code=True)
        model = AutoModelForMaskedLM.from_pretrained(MODEL_ID, trust_remote_code=True)
        model = model.to(DEVICE).eval()
        return tok, model
    except ImportError as e:
        print(f"  Caduceus model load failed: {e}", file=sys.stderr)
        return None, None


def predict_masked(model, tok, original_seq: str, mask_positions: list):
    """Predict the original nucleotide at each masked position via [MASK].
    Caduceus uses single-character tokens (one token per base). The MASK token
    occupies the same position as the masked base — straightforward MLM.
    Logits have 16 columns; filter to ACGT (ids 7-10).
    """
    predictions = []
    mask_id = tok.mask_token_id

    for pos in mask_positions:
        seq = list(original_seq.upper())
        seq[pos] = tok.mask_token
        masked_str = "".join(seq)

        enc = tok(masked_str, return_tensors="pt", add_special_tokens=False)
        ids = enc.input_ids.to(DEVICE)

        with torch.no_grad():
            out = model(input_ids=ids)
            logits = out.logits if hasattr(out, "logits") else out[0]

        mask_token_indices = (ids == mask_id).nonzero(as_tuple=True)[1]
        if len(mask_token_indices) == 0:
            pred_base, pred_prob = "N", 0.0
        else:
            tok_idx = int(mask_token_indices[0])
            tok_logits = logits[0, tok_idx]  # [vocab_or_more]
            # Restrict to ACGT subset
            acgt_logits = torch.tensor([tok_logits[i].item() for i in ACGT_IDS], device=tok_logits.device)
            probs = torch.softmax(acgt_logits, dim=-1)
            top_idx = int(torch.argmax(probs))
            id_list = list(ACGT_IDS.keys())
            pred_base = ACGT_IDS[id_list[top_idx]]
            pred_prob = float(probs[top_idx])

        orig_base = original_seq[pos].upper() if pos < len(original_seq) else "N"
        predictions.append({
            "position": pos,
            "original_base": orig_base,
            "predicted_base": pred_base,
            "prob": round(pred_prob, 4),
            "correct": pred_base == orig_base,
        })

    return predictions


def main():
    print("=== T1.1 Caduceus Mask & Fill ===")
    print(f"Loading {MODEL_ID} on {DEVICE}...")
    tok, model = load_model()
    if model is None:
        results = []
        with open(DATA_FILE) as f:
            for item in json.load(f):
                results.append({
                    "genome_name": item["genome_name"],
                    "strategy": item["strategy"],
                    "model": "Caduceus",
                    "model_id": MODEL_ID,
                    "predictions": [],
                    "recovery_rate": 0.0,
                    "n_masks": len(item["mask_positions"]),
                    "n_correct": 0,
                    "caduceus_used": False,
                    "skip_reason": "mamba_ssm package not installed (CUDA compile required)",
                })
        outpath = RESULTS_DIR / "caduceus_predictions.json"
        json.dump(results, open(outpath, "w"), indent=2)
        print(f"  Caduceus skipped (mamba_ssm missing). Stub saved → {outpath}")
        return

    print(f"  Caduceus loaded: {sum(p.numel() for p in model.parameters())/1e6:.1f}M params")

    with open(DATA_FILE) as f:
        masked_data = json.load(f)

    results = []
    for item in masked_data:
        genome = item["genome_name"]
        strategy = item["strategy"]
        orig_seq = item["original_seq"]
        mask_positions = item["mask_positions"]
        try:
            preds = predict_masked(model, tok, orig_seq, mask_positions)
        except Exception as e:
            print(f"  Error on {genome}/{strategy}: {e}", file=sys.stderr)
            preds = []
        n_correct = sum(1 for p in preds if p["correct"])
        recovery = n_correct / len(preds) if preds else 0.0
        print(f"  {genome:30s} {strategy:15s}: recovery={recovery:.3f} ({n_correct}/{len(preds)})")
        results.append({
            "genome_name": genome,
            "strategy": strategy,
            "model": "Caduceus-ph_seqlen-1k_d_model-118_n_layer-4",
            "model_id": MODEL_ID,
            "predictions": preds,
            "recovery_rate": round(recovery, 4),
            "n_masks": len(mask_positions),
            "n_correct": n_correct,
            "caduceus_used": True,
        })

    outpath = RESULTS_DIR / "caduceus_predictions.json"
    json.dump(results, open(outpath, "w"), indent=2)
    print(f"\nSaved Caduceus predictions → {outpath}")
    print("\nRecovery summary:")
    for r in results:
        print(f"  {r['genome_name']:30s} {r['strategy']:15s}: Caduceus={r['recovery_rate']:.3f}")


if __name__ == "__main__":
    main()
