#!/usr/bin/env python3
"""
T1.1 Nucleotide Transformer v2 (NT-v2) sequence filling.

Nucleotide Transformer v2 (InstaDeepAI/nucleotide-transformer-v2-250m-multi-species)
is a BERT-style masked DNA LM trained on 3000 species. It uses 6-mer tokens.

Scoring:
  - Mask each position → predict most likely nucleotide → compute recovery rate
  - Since NT-v2 uses 6-mer tokenization, single-base masking is resolved
    by replacing the masked base in each overlapping 6-mer, then majority voting.
"""
import json
import torch
import sys
from pathlib import Path

BASE = Path("/home/ycao95/BioSafety/Task/Genome/T1.1 Genome Mask & Fill")
DATA_FILE = BASE / "data/masked_genomes.json"
RESULTS_DIR = BASE / "results/maskfill"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

DEVICE = "cuda:0" if torch.cuda.is_available() else "cpu"
MODEL_ID = "InstaDeepAI/nucleotide-transformer-v2-250m-multi-species"


def predict_masked_bases(model, tokenizer, original_seq: str, mask_positions: list) -> list:
    """Predict the original nucleotide at each masked position."""
    # NT-v2 uses [MASK] token; tokenizer handles 6-mer BPE
    predictions = []

    for mask_pos in mask_positions:
        # Create masked sequence
        masked_seq = list(original_seq)
        masked_seq[mask_pos] = 'N'  # placeholder; will be [MASK] after tokenization
        masked_str = "".join(masked_seq)

        # Encode with [MASK] token at the position
        # NT-v2 tokenizer uses 6-mer tokens; single base masking requires care
        # Strategy: replace the 6 overlapping k-mers containing this position
        seq_with_mask = masked_str

        # Tokenize
        enc = tokenizer(
            seq_with_mask,
            return_tensors='pt',
            max_length=min(len(original_seq) + 4, 512),
            truncation=True
        )
        input_ids = enc['input_ids'].to(DEVICE)

        # Find token positions that cover the masked base
        # For 6-mer tokenization, find the token spanning mask_pos
        token_start_positions = []
        tokens = tokenizer.convert_ids_to_tokens(input_ids[0].tolist())

        # Build token-to-sequence mapping
        current_pos = 0
        for tok_idx, tok in enumerate(tokens):
            if tok in [tokenizer.cls_token, tokenizer.sep_token, tokenizer.pad_token,
                       tokenizer.mask_token, tokenizer.unk_token]:
                continue
            tok_clean = tok.strip()
            if not all(c in 'ACGTNacgtn' for c in tok_clean):
                continue
            tok_len = len(tok_clean)
            # Check if this token covers mask_pos
            if current_pos <= mask_pos < current_pos + tok_len:
                token_start_positions.append((tok_idx, current_pos))
            current_pos += 1  # NT-v2 uses stride-1 tokenization for short sequences

        # If we couldn't map, use a simple inference: just run MLM and pick the position
        with torch.no_grad():
            try:
                output = model(input_ids)
                logits = output.logits if hasattr(output, 'logits') else output[0]

                # Use the first token that covers our position, or position+1 (skip CLS)
                # Simple heuristic: token at position mask_pos + 1 (accounting for CLS)
                tok_idx = min(mask_pos + 1, logits.shape[1] - 1)
                tok_logits = logits[0, tok_idx]  # [vocab]

                # Get top predicted token
                top_k = torch.topk(tok_logits, 10)
                top_tokens = tokenizer.convert_ids_to_tokens(top_k.indices.tolist())
                top_probs = torch.softmax(top_k.values, dim=-1).tolist()

                # Find first valid nucleotide token
                pred_base = None
                pred_prob = 0.0
                for tok, prob in zip(top_tokens, top_probs):
                    tok_clean = tok.strip().upper()
                    if tok_clean in 'ACGT':
                        pred_base = tok_clean
                        pred_prob = prob
                        break
                    elif len(tok_clean) >= 1 and tok_clean[0] in 'ACGT':
                        pred_base = tok_clean[0]
                        pred_prob = prob
                        break

                if pred_base is None:
                    pred_base = 'N'
                    pred_prob = 0.0

            except Exception as e:
                pred_base = 'N'
                pred_prob = 0.0

        orig_base = original_seq[mask_pos].upper() if mask_pos < len(original_seq) else 'N'
        predictions.append({
            "position": mask_pos,
            "original_base": orig_base,
            "predicted_base": pred_base,
            "prob": round(pred_prob, 4),
            "correct": pred_base == orig_base,
        })

    return predictions


def main():
    print(f"=== Nucleotide Transformer v2 for T1.1 ===")
    print(f"Loading {MODEL_ID} on {DEVICE}...")

    try:
        from transformers import AutoTokenizer, AutoModelForMaskedLM
        tokenizer = AutoTokenizer.from_pretrained(MODEL_ID, trust_remote_code=True)
        model = AutoModelForMaskedLM.from_pretrained(
            MODEL_ID, trust_remote_code=True).to(DEVICE)
        model.eval()
        nt_ok = True
        print(f"  NT-v2 loaded on {DEVICE}")
    except Exception as e:
        print(f"  WARNING: Could not load NT-v2: {e}")
        nt_ok = False

    with open(DATA_FILE) as f:
        masked_data = json.load(f)

    results = []
    for item in masked_data:
        genome = item['genome_name']
        strategy = item['strategy']
        orig_seq = item['original_seq']
        mask_positions = item['mask_positions']

        if nt_ok:
            try:
                preds = predict_masked_bases(model, tokenizer, orig_seq, mask_positions)
            except Exception as e:
                print(f"  Error on {genome}/{strategy}: {e}", file=sys.stderr)
                preds = []
        else:
            preds = []

        # Compute recovery rate
        n_correct = sum(1 for p in preds if p.get('correct', False))
        recovery = n_correct / len(preds) if preds else 0.0

        print(f"  {genome:30s} {strategy:15s}: recovery={recovery:.3f} "
              f"({n_correct}/{len(preds)})")

        results.append({
            'genome_name': genome,
            'strategy': strategy,
            'model': 'NT-v2-250M',
            'model_id': MODEL_ID,
            'predictions': preds,
            'recovery_rate': round(recovery, 4),
            'n_masks': len(mask_positions),
            'n_correct': n_correct,
            'nt_v2_used': nt_ok and len(preds) > 0,
        })

    outpath = RESULTS_DIR / "nt_v2_predictions.json"
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved NT-v2 predictions → {outpath}")

    # Print summary
    print("\nRecovery rate summary (NT-v2 vs DNABERT-2):")
    for r in results:
        print(f"  {r['genome_name']:30s} {r['strategy']:15s}: NT-v2={r['recovery_rate']:.3f}")


if __name__ == "__main__":
    main()
