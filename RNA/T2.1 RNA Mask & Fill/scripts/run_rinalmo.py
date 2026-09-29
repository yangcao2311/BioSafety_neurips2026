#!/usr/bin/env python3
"""
T2.1 RiNALMo Mask & Fill — closes the model gap (RNA-FM ✓ + SpliceBERT ✓ + UTR-LM ✓ + RiNALMo ✓).

RiNALMo (Pénić et al. 2024, NeurIPS): RNA Language Model trained on 36M RNA
sequences with single-base tokenization. The "giga" variant is 650M params.

Source: multimolecule/rinalmo-giga (HF mirror of lbcb-sci/RiNALMo private weights)
"""
import json
import sys
import torch
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
DATA_FILE = BASE / "data" / "sequences" / "rna_sequences.json"
RESULTS_DIR = BASE / "results" / "maskfill"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

DEVICE = "cuda:0" if torch.cuda.is_available() else "cpu"
MODEL_ID = "multimolecule/rinalmo-giga"

BASES = set("AUGC")


def predict_masks_mlm(model, tokenizer, rna_seq: str, mask_positions: list) -> list:
    predictions = []
    mask_token = tokenizer.mask_token
    mask_token_id = tokenizer.mask_token_id

    for pos in mask_positions:
        seq_list = list(rna_seq)
        if pos < len(seq_list):
            seq_list[pos] = mask_token
        masked_str = "".join(seq_list)

        enc = tokenizer(masked_str, return_tensors="pt",
                        max_length=min(len(rna_seq) + 4, 4096),
                        truncation=True)
        input_ids = enc["input_ids"].to(DEVICE)

        with torch.no_grad():
            out = model(**{k: v.to(DEVICE) for k, v in enc.items()})
            logits = out.logits

        mask_idx = (input_ids[0] == mask_token_id).nonzero(as_tuple=True)[0]
        if len(mask_idx) == 0:
            predictions.append({"position": pos,
                                "original_base": rna_seq[pos].upper() if pos < len(rna_seq) else "N",
                                "predicted_base": "N", "prob": 0.0, "correct": False})
            continue

        tok_idx = int(mask_idx[0])
        tok_logits = logits[0, tok_idx]
        top_k = torch.topk(tok_logits, 10)
        top_tokens = tokenizer.convert_ids_to_tokens(top_k.indices.tolist())
        top_probs = torch.softmax(top_k.values, dim=-1).tolist()

        pred_base, pred_prob = "N", 0.0
        for tk, prob in zip(top_tokens, top_probs):
            tk_clean = tk.strip().upper().replace("T", "U")
            if tk_clean in BASES:
                pred_base, pred_prob = tk_clean, prob
                break
            elif len(tk_clean) >= 1 and tk_clean[0] in BASES:
                pred_base, pred_prob = tk_clean[0], prob
                break

        orig_base = rna_seq[pos].upper() if pos < len(rna_seq) else "N"
        predictions.append({
            "position": pos,
            "original_base": orig_base,
            "predicted_base": pred_base,
            "prob": round(pred_prob, 4),
            "correct": pred_base == orig_base,
        })
    return predictions


def main():
    print(f"=== T2.1 RiNALMo Mask & Fill ===")
    print(f"Loading {MODEL_ID} on {DEVICE}...")
    import os
    os.environ.setdefault("HF_HOME", "/home/xliu316/.cache/huggingface")
    from multimolecule import RnaTokenizer
    from transformers import AutoModelForMaskedLM
    tokenizer = RnaTokenizer.from_pretrained(MODEL_ID)
    model = AutoModelForMaskedLM.from_pretrained(MODEL_ID, trust_remote_code=True).to(DEVICE)
    model.eval()
    print(f"  RiNALMo loaded: {sum(p.numel() for p in model.parameters())/1e6:.1f}M params")

    with open(DATA_FILE) as f:
        masked_data = json.load(f)

    results = []
    for item in masked_data:
        rna = item.get("rna_name", item.get("genome_name", "unknown"))
        strategy = item["strategy"]
        rna_seq = item["original_seq"].upper().replace("T", "U")
        mask_positions = item.get("mask_positions", [])
        if mask_positions:
            preds = predict_masks_mlm(model, tokenizer, rna_seq, mask_positions)
        else:
            preds = []
        n_correct = sum(1 for p in preds if p["correct"])
        recovery = n_correct / len(preds) if preds else 0.0
        print(f"  {rna:30s} {strategy:18s}: recovery={recovery:.3f} ({n_correct}/{len(preds)})")
        results.append({
            "rna_name": rna,
            "strategy": strategy,
            "model": "RiNALMo-giga",
            "model_id": MODEL_ID,
            "predictions": preds,
            "recovery_rate": round(recovery, 4),
            "n_masks": len(mask_positions),
            "n_correct": n_correct,
            "rinalmo_used": True,
        })

    out = RESULTS_DIR / "rinalmo_predictions.json"
    json.dump(results, open(out, "w"), indent=2)
    print(f"\nSaved → {out}")
    print("\nSummary:")
    for r in results:
        print(f"  {r['rna_name']:30s} {r['strategy']:18s}: RiNALMo={r['recovery_rate']:.3f}")


if __name__ == "__main__":
    main()
