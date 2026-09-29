#!/usr/bin/env python3
"""
T2.1 RNA Mask & Fill — SpliceBERT + UTR-LM additional models.

SpliceBERT: RNA language model specialized for splice-site prediction.
UTR-LM: 3-mer tokenized UTR language model.

Both use BERT-style masked language modeling via multimolecule package.
"""
import json
import sys
import torch
from pathlib import Path

BASE = Path("/home/ycao95/BioSafety/Task/RNA/T2.1 RNA Mask & Fill")
DATA_FILE = BASE / "data/sequences/rna_sequences.json"
RESULTS_DIR = BASE / "results/maskfill"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

DEVICE = "cuda:1" if torch.cuda.is_available() and torch.cuda.device_count() > 1 else \
         "cuda:0" if torch.cuda.is_available() else "cpu"

SPLICEBERT_ID = "multimolecule/splicebert"
UTRLM_ID = "multimolecule/utrbert-3mer"

BASES = set("AUGC")


def predict_masks_mlm(model, tokenizer, rna_seq: str, mask_positions: list,
                      model_name: str) -> list:
    """Run masked LM inference on a list of positions."""
    predictions = []
    mask_token = tokenizer.mask_token or "<mask>"
    mask_token_id = tokenizer.mask_token_id

    for pos in mask_positions:
        seq_list = list(rna_seq)
        if pos < len(seq_list):
            seq_list[pos] = mask_token
        masked_str = "".join(seq_list)

        enc = tokenizer(
            masked_str,
            return_tensors='pt',
            max_length=min(len(rna_seq) + 4, 512),
            truncation=True
        )
        input_ids = enc['input_ids'].to(DEVICE)

        with torch.no_grad():
            output = model(**{k: v.to(DEVICE) for k, v in enc.items()})
            logits = output.logits

        # Find mask token position
        mask_positions_in_tokens = (input_ids[0] == mask_token_id).nonzero(as_tuple=True)[0]
        if len(mask_positions_in_tokens) == 0:
            predictions.append({
                "position": pos,
                "original_base": rna_seq[pos] if pos < len(rna_seq) else 'N',
                "predicted_base": 'N',
                "prob": 0.0,
                "correct": False,
            })
            continue

        tok_idx = mask_positions_in_tokens[0].item()
        tok_logits = logits[0, tok_idx]

        top_k = torch.topk(tok_logits, 10)
        top_tokens = tokenizer.convert_ids_to_tokens(top_k.indices.tolist())
        top_probs = torch.softmax(top_k.values, dim=-1).tolist()

        # Find valid RNA base prediction
        pred_base = None
        pred_prob = 0.0
        for tok, prob in zip(top_tokens, top_probs):
            tok_clean = tok.strip().upper().replace('T', 'U')
            if tok_clean in BASES:
                pred_base = tok_clean
                pred_prob = prob
                break
            elif len(tok_clean) >= 1 and tok_clean[0] in BASES:
                pred_base = tok_clean[0]
                pred_prob = prob
                break

        if pred_base is None:
            pred_base = 'N'
            pred_prob = 0.0

        orig_base = rna_seq[pos].upper() if pos < len(rna_seq) else 'N'
        predictions.append({
            "position": pos,
            "original_base": orig_base,
            "predicted_base": pred_base,
            "prob": round(pred_prob, 4),
            "correct": pred_base == orig_base,
        })

    return predictions


def run_model(model_id, model_class, tokenizer_class, model_name, masked_data):
    print(f"\n{'='*60}")
    print(f"Loading {model_name} ({model_id})")
    try:
        tokenizer = tokenizer_class.from_pretrained(model_id)
        model = model_class.from_pretrained(model_id).to(DEVICE)
        model.eval()
        print(f"  {model_name} loaded on {DEVICE}")
        model_ok = True
    except Exception as e:
        print(f"  Failed to load {model_name}: {e}")
        model_ok = False
        model = tokenizer = None

    results = []
    for item in masked_data:
        genome = item.get('rna_name', item.get('genome_name', 'unknown'))
        strategy = item['strategy']
        rna_seq = item['original_seq'].upper().replace('T', 'U')
        mask_positions = item.get('mask_positions', [])

        if model_ok and mask_positions:
            try:
                preds = predict_masks_mlm(model, tokenizer, rna_seq, mask_positions, model_name)
            except Exception as e:
                print(f"  Error on {genome}/{strategy}: {e}", file=sys.stderr)
                preds = []
        else:
            preds = []

        n_correct = sum(1 for p in preds if p.get('correct', False))
        recovery = n_correct / len(preds) if preds else 0.0

        print(f"  {genome:35s} {strategy:15s}: recovery={recovery:.3f} "
              f"({n_correct}/{len(preds)})")

        results.append({
            'genome_name': genome,
            'strategy': strategy,
            'model': model_name,
            'model_id': model_id,
            'predictions': preds,
            'recovery_rate': round(recovery, 4),
            'n_masks': len(mask_positions),
            'n_correct': n_correct,
            'model_used': model_ok and len(preds) > 0,
        })

    if model_ok:
        del model
        torch.cuda.empty_cache()

    return results


def main():
    print("=== T2.1 Additional RNA Models: SpliceBERT + UTR-LM ===")
    with open(DATA_FILE) as f:
        masked_data = json.load(f)
    print(f"Loaded {len(masked_data)} RNA sequences")

    from multimolecule import SpliceBertForMaskedLM, UtrLmForMaskedLM, AutoTokenizer

    # Run SpliceBERT
    splicebert_results = run_model(
        SPLICEBERT_ID, SpliceBertForMaskedLM, AutoTokenizer,
        "SpliceBERT", masked_data
    )
    out_splice = RESULTS_DIR / "splicebert_predictions.json"
    with open(out_splice, 'w') as f:
        json.dump(splicebert_results, f, indent=2)
    print(f"Saved SpliceBERT → {out_splice}")

    # Run UTR-LM
    utrlm_results = run_model(
        UTRLM_ID, UtrLmForMaskedLM, AutoTokenizer,
        "UTR-LM", masked_data
    )
    out_utr = RESULTS_DIR / "utrlm_predictions.json"
    with open(out_utr, 'w') as f:
        json.dump(utrlm_results, f, indent=2)
    print(f"Saved UTR-LM → {out_utr}")

    # Summary comparison
    print("\n=== Multi-Model RNA Recovery Comparison ===")

    # Load existing RNA-FM results
    rnafm_path = RESULTS_DIR / "rnafm_predictions.json"
    if rnafm_path.exists():
        with open(rnafm_path) as f:
            rnafm_results = json.load(f)
        rnafm_avg = sum(r.get('recovery_rate', 0) for r in rnafm_results) / max(1, len(rnafm_results))
        print(f"  RNA-FM average recovery:    {rnafm_avg:.3f}")

    splice_avg = sum(r.get('recovery_rate', 0) for r in splicebert_results) / max(1, len(splicebert_results))
    utr_avg = sum(r.get('recovery_rate', 0) for r in utrlm_results) / max(1, len(utrlm_results))
    print(f"  SpliceBERT average recovery: {splice_avg:.3f}")
    print(f"  UTR-LM average recovery:     {utr_avg:.3f}")

    # Save combined summary
    summary = {
        "models": {
            "RNA-FM": {
                "model_id": "multimolecule/rnafm",
                "avg_recovery": round(rnafm_avg, 4) if rnafm_path.exists() else None
            },
            "SpliceBERT": {
                "model_id": SPLICEBERT_ID,
                "avg_recovery": round(splice_avg, 4)
            },
            "UTR-LM": {
                "model_id": UTRLM_ID,
                "avg_recovery": round(utr_avg, 4)
            }
        }
    }
    out_summary = RESULTS_DIR / "multimodel_comparison.json"
    with open(out_summary, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"Saved multi-model comparison → {out_summary}")


if __name__ == "__main__":
    main()
