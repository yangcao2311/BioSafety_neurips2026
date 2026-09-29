#!/usr/bin/env python3
"""
T3.1 Additional protein models: ProtBERT (MLM), ProtT5 (encoder), Ankh (encoder)
- ProtBERT: BERT-style → AutoModelForMaskedLM
- ProtT5: T5 encoder → T5EncoderModel + linear head for residue prediction
- Ankh: T5-like encoder → same approach as ProtT5
"""
import json, torch
from pathlib import Path

DEVICE = "cuda:0" if torch.cuda.is_available() else "cpu"
BASE = Path("/home/ycao95/BioSafety/Task/Protein/T3.1 Protein Mask & Fill")
MASKS_DIR  = BASE / "data/masks"
RESULTS_DIR = BASE / "results/maskfill"

def load_masks():
    return [json.load(open(f)) for f in sorted(MASKS_DIR.glob("*.json"))
            if f.name != "mask_manifest.json"]

# ──────────────────────────────────────────────────────────────────────────────
# ProtBERT  (BERT-style MLM)
# ──────────────────────────────────────────────────────────────────────────────
def run_protbert(masks):
    from transformers import BertForMaskedLM, BertTokenizer
    print("\n" + "="*60 + "\nLoading ProtBERT (Rostlab/prot_bert)")
    tokenizer = BertTokenizer.from_pretrained("Rostlab/prot_bert", do_lower_case=False)
    model = BertForMaskedLM.from_pretrained("Rostlab/prot_bert").to(DEVICE)
    model.eval()
    mask_token = tokenizer.mask_token   # "[MASK]"
    out_dir = RESULTS_DIR / "ProtBERT"; out_dir.mkdir(exist_ok=True)

    for entry in masks:
        protein  = entry.get("protein_name", entry.get("name", "unknown"))
        strategy = entry.get("strategy", "unknown")
        orig_seq = entry.get("original_seq", "")
        mask_pos = entry.get("mask_positions_0indexed", [])

        # ProtBERT needs space-separated single-char AA
        spaced = list(orig_seq)
        for p in mask_pos:
            if p < len(spaced):
                spaced[p] = mask_token
        seq_str = " ".join(spaced)

        inputs = tokenizer(seq_str, return_tensors="pt",
                           truncation=True, max_length=512).to(DEVICE)
        with torch.no_grad():
            logits = model(**inputs).logits[0]   # [L, vocab]

        input_ids = inputs["input_ids"][0]
        mask_id   = tokenizer.mask_token_id
        mask_tok_positions = (input_ids == mask_id).nonzero(as_tuple=True)[0].tolist()

        predictions = []
        for tok_idx, orig_pos in zip(mask_tok_positions[:len(mask_pos)], mask_pos):
            top5  = torch.topk(logits[tok_idx], 5)
            toks  = [tokenizer.convert_ids_to_tokens([i.item()])[0].strip()
                     for i in top5.indices]
            probs = top5.values.softmax(-1).tolist()
            orig_aa = orig_seq[orig_pos] if orig_pos < len(orig_seq) else "?"
            predictions.append({
                "position": orig_pos,
                "original_aa": orig_aa,
                "predicted_aa": toks[0],
                "prob": round(probs[0], 4),
                "correct": toks[0].upper() == orig_aa.upper(),
                "top5": [{"aa": t, "prob": round(p, 4)} for t, p in zip(toks, probs)]
            })

        result = {"protein_name": protein, "strategy": strategy,
                  "model": "ProtBERT", "predictions": predictions}
        json.dump(result, open(out_dir / f"{protein}_{strategy}_predictions.json","w"), indent=2)
        print(f"  [{protein}/{strategy}] {len(predictions)} preds", flush=True)

    del model; torch.cuda.empty_cache()
    print("ProtBERT done")

# ──────────────────────────────────────────────────────────────────────────────
# ProtT5  (T5 encoder → per-residue logits via shared embedding)
# ──────────────────────────────────────────────────────────────────────────────
def run_prot_t5(masks):
    from transformers import T5EncoderModel, T5Tokenizer
    print("\n" + "="*60 + "\nLoading ProtT5 (Rostlab/prot_t5_xl_uniref50)")
    tokenizer = T5Tokenizer.from_pretrained("Rostlab/prot_t5_xl_uniref50",
                                             do_lower_case=False, legacy=True)
    model = T5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_uniref50").to(DEVICE)
    model.eval()
    embed_dim = model.config.d_model   # 1024 for XL
    # Use shared embedding matrix as linear classifier (dimension: vocab × d_model)
    shared_weight = model.shared.weight.detach()   # [vocab, d_model]
    out_dir = RESULTS_DIR / "ProtT5"; out_dir.mkdir(exist_ok=True)

    AA_VOCAB = list("ACDEFGHIKLMNPQRSTVWY")

    for entry in masks:
        protein  = entry.get("protein_name", entry.get("name","unknown"))
        strategy = entry.get("strategy","unknown")
        orig_seq = entry.get("original_seq","")
        mask_pos = entry.get("mask_positions_0indexed",[])

        # ProtT5 format: space-separated AAs; replace masked positions with X
        spaced = list(orig_seq)
        for p in mask_pos:
            if p < len(spaced):
                spaced[p] = "X"   # unknown token used as mask proxy
        seq_str = " ".join(spaced)

        inputs = tokenizer(seq_str, return_tensors="pt",
                           truncation=True, max_length=512,
                           add_special_tokens=True).to(DEVICE)
        with torch.no_grad():
            hidden = model(**inputs).last_hidden_state[0]  # [L, 1024]

        # Score each masked position against shared embedding
        predictions = []
        # offset 1 for <s> special token
        for orig_pos in mask_pos:
            tok_idx = orig_pos + 1   # account for BOS
            if tok_idx >= hidden.shape[0]:
                continue
            h = hidden[tok_idx]  # [d_model]
            scores = shared_weight @ h   # [vocab]
            top5_idx = torch.topk(scores, 5).indices.tolist()
            top5_probs = torch.softmax(torch.topk(scores, 5).values, -1).tolist()
            top5_toks = [tokenizer.convert_ids_to_tokens([i])[0].strip()
                         for i in top5_idx]
            # Filter to single AA tokens only
            aa_pred = next((t for t in top5_toks if t.upper() in AA_VOCAB), top5_toks[0])
            orig_aa = orig_seq[orig_pos] if orig_pos < len(orig_seq) else "?"
            predictions.append({
                "position": orig_pos,
                "original_aa": orig_aa,
                "predicted_aa": aa_pred,
                "prob": round(top5_probs[0], 4),
                "correct": aa_pred.upper() == orig_aa.upper(),
                "top5": [{"aa": t, "prob": round(p, 4)}
                         for t, p in zip(top5_toks, top5_probs)]
            })

        result = {"protein_name": protein, "strategy": strategy,
                  "model": "ProtT5", "predictions": predictions}
        json.dump(result, open(out_dir / f"{protein}_{strategy}_predictions.json","w"), indent=2)
        print(f"  [{protein}/{strategy}] {len(predictions)} preds", flush=True)

    del model; torch.cuda.empty_cache()
    print("ProtT5 done")

# ──────────────────────────────────────────────────────────────────────────────
# Ankh  (T5-like, smaller, better AA tokenisation)
# ──────────────────────────────────────────────────────────────────────────────
def run_ankh(masks):
    from transformers import T5EncoderModel, AutoTokenizer
    print("\n" + "="*60 + "\nLoading Ankh (ElnaggarLab/ankh-base)")
    tokenizer = AutoTokenizer.from_pretrained("ElnaggarLab/ankh-base", use_fast=True)
    model = T5EncoderModel.from_pretrained("ElnaggarLab/ankh-base").to(DEVICE)
    model.eval()
    shared_weight = model.shared.weight.detach()
    AA_VOCAB = list("ACDEFGHIKLMNPQRSTVWY")
    out_dir = RESULTS_DIR / "Ankh"; out_dir.mkdir(exist_ok=True)

    for entry in masks:
        protein  = entry.get("protein_name", entry.get("name","unknown"))
        strategy = entry.get("strategy","unknown")
        orig_seq = entry.get("original_seq","")
        mask_pos = entry.get("mask_positions_0indexed",[])

        aa_seq = list(orig_seq)
        for p in mask_pos:
            if p < len(aa_seq):
                aa_seq[p] = "X"
        seq_str = " ".join(aa_seq)

        inputs = tokenizer(seq_str, return_tensors="pt",
                           truncation=True, max_length=512).to(DEVICE)
        with torch.no_grad():
            hidden = model(**inputs).last_hidden_state[0]

        predictions = []
        for orig_pos in mask_pos:
            tok_idx = orig_pos + 1
            if tok_idx >= hidden.shape[0]:
                continue
            h = hidden[tok_idx]
            scores = shared_weight @ h
            top5_idx = torch.topk(scores, 5).indices.tolist()
            top5_probs = torch.softmax(torch.topk(scores, 5).values, -1).tolist()
            top5_toks = [tokenizer.convert_ids_to_tokens([i])[0].strip() for i in top5_idx]
            aa_pred = next((t for t in top5_toks if t.upper() in AA_VOCAB), top5_toks[0])
            orig_aa = orig_seq[orig_pos] if orig_pos < len(orig_seq) else "?"
            predictions.append({
                "position": orig_pos,
                "original_aa": orig_aa,
                "predicted_aa": aa_pred,
                "prob": round(top5_probs[0], 4),
                "correct": aa_pred.upper() == orig_aa.upper(),
                "top5": [{"aa": t, "prob": round(p, 4)}
                         for t, p in zip(top5_toks, top5_probs)]
            })

        result = {"protein_name": protein, "strategy": strategy,
                  "model": "Ankh", "predictions": predictions}
        json.dump(result, open(out_dir / f"{protein}_{strategy}_predictions.json","w"), indent=2)
        print(f"  [{protein}/{strategy}] {len(predictions)} preds", flush=True)

    del model; torch.cuda.empty_cache()
    print("Ankh done")

# ──────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    masks = load_masks()
    print(f"Loaded {len(masks)} mask files | Device: {DEVICE}")
    run_protbert(masks)
    run_prot_t5(masks)
    run_ankh(masks)
    print("\nAll extra protein models DONE")
