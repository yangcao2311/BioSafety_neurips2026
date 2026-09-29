#!/usr/bin/env python3
"""
T4.1 ChemBERTa-77M-MLM mask fill.

Runs the spec-canonical DeepChem/ChemBERTa-77M-MLM (3.4M param model trained
on 77M PubChem SMILES with masked language modeling). The pre-existing run
used seyonec/ChemBERTa-zinc-base-v1, which is a different ChemBERTa variant.
This script targets the exact model identifier listed in
BioSafety_Benchmark.md section 5.2.

Output: results/maskfill/chemberta_77m_mlm_predictions.json with the same
schema as the other ChemBERTa output so the downstream cascade can ingest
either source.
"""
import json
import os
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import torch

BASE = Path(__file__).resolve().parent.parent
DATA = BASE / "data" / "tox21" / "masked_smiles.json"
RESULTS = BASE / "results" / "maskfill"
RESULTS.mkdir(parents=True, exist_ok=True)
OUT = RESULTS / "chemberta_77m_mlm_predictions.json"

DEVICE = "cuda:0" if torch.cuda.is_available() else "cpu"
MODEL_ID = "DeepChem/ChemBERTa-77M-MLM"


def main():
    print(f"=== T4.1 ChemBERTa-77M-MLM ===")
    print(f"Loading {MODEL_ID} on {DEVICE}")
    os.environ.setdefault("HF_HOME", "/home/xliu316/.cache/huggingface")
    from transformers import AutoTokenizer, AutoModelForMaskedLM
    tok = AutoTokenizer.from_pretrained(MODEL_ID)
    model = AutoModelForMaskedLM.from_pretrained(MODEL_ID).to(DEVICE).eval()
    n_params = sum(p.numel() for p in model.parameters())
    print(f"  Loaded ({n_params/1e6:.1f}M params)")
    mask_token = tok.mask_token
    mask_id = tok.mask_token_id

    with open(DATA) as f:
        data = json.load(f)

    out_records = []
    for rec in data:
        compound = rec.get("name") or rec.get("compound_name") or "unknown"
        original = rec.get("smiles") or ""
        strategies = rec.get("strategies", {})
        if not strategies:
            continue
        for strategy_name, strat in strategies.items():
            masked_smiles = strat.get("masked_smiles")
            if masked_smiles is None or "[MASK]" not in masked_smiles:
                continue
            n_masks = masked_smiles.count("[MASK]")

            masked_for_tok = masked_smiles.replace("[MASK]", mask_token)
            try:
                enc = tok(masked_for_tok, return_tensors="pt",
                          max_length=tok.model_max_length, truncation=True)
                ids = enc.input_ids.to(DEVICE)
                with torch.no_grad():
                    out = model(**{k: v.to(DEVICE) for k, v in enc.items()})
                logits = out.logits
                mask_positions = (ids[0] == mask_id).nonzero(as_tuple=True)[0]
            except Exception as e:
                print(f"  {compound}/{strategy_name}: {e}")
                continue

            predictions = []
            filled = masked_smiles
            for i, pos in enumerate(mask_positions):
                tok_logits = logits[0, int(pos)]
                top5 = torch.topk(tok_logits, 5)
                top5_probs = torch.softmax(top5.values, dim=-1).tolist()
                top5_tokens = tok.convert_ids_to_tokens(top5.indices.tolist())
                best = top5_tokens[0]
                best_clean = best.lstrip("Ġ▁")
                predictions.append({
                    "mask_index": i,
                    "predicted_token": best_clean,
                    "probability": round(top5_probs[0], 4),
                    "top5": [{"token": t.lstrip("Ġ▁"), "prob": round(p, 4)}
                             for t, p in zip(top5_tokens, top5_probs)],
                })
                filled = filled.replace("[MASK]", best_clean, 1)

            out_records.append({
                "name": compound,
                "strategy": strategy_name,
                "model": "DeepChem/ChemBERTa-77M-MLM",
                "original_smiles": original,
                "masked_smiles": masked_smiles,
                "filled_smiles": filled,
                "n_masks": n_masks,
                "predictions": predictions,
            })

    json.dump(out_records, open(OUT, "w"), indent=2)
    print(f"\nSaved {len(out_records)} records to {OUT}")


if __name__ == "__main__":
    main()
