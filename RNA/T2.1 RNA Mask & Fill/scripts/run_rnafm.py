#!/usr/bin/env python3
"""
T2.1 RNA-FM Mask & Fill using multimolecule package.
"""
import json
import os
import subprocess
import tempfile
import sys

BASE = "/home/ycao95/BioSafety/Task/RNA/T2.1 RNA Mask & Fill"
DATA_FILE = os.path.join(BASE, "data/sequences/rna_sequences.json")
RESULTS_DIR = os.path.join(BASE, "results/maskfill")
os.makedirs(RESULTS_DIR, exist_ok=True)

BIOSAFETY_PYTHON = "/home/ycao95/.conda/envs/biosafety/bin/python"

def run_rnafm_subprocess(data_file, results_dir):
    """Run RNA-FM inference in biosafety conda env."""
    script = f"""
import json
import os
import torch

DATA_FILE = {json.dumps(data_file)}
RESULTS_DIR = {json.dumps(results_dir)}

from multimolecule import RnaTokenizer, RnaFmForMaskedLM

DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
print(f"RNA-FM on {{DEVICE}}")

tokenizer = RnaTokenizer.from_pretrained("multimolecule/rnafm")
model = RnaFmForMaskedLM.from_pretrained("multimolecule/rnafm").to(DEVICE)
model.eval()

with open(DATA_FILE) as f:
    data = json.load(f)

RNA_BASES = set("AUGCN")

results = []
for item in data:
    rna_name = item["rna_name"]
    strategy = item["strategy"]
    orig_seq = item["original_seq"].upper().replace("T", "U")
    mask_positions = item["mask_positions"]
    
    # Build masked sequence
    seq_list = list(orig_seq)
    for pos in mask_positions:
        if pos < len(seq_list):
            seq_list[pos] = tokenizer.mask_token
    masked_seq = "".join(seq_list)
    
    # Tokenize (truncate if too long)
    max_len = min(len(masked_seq) + 10, 512)
    enc = tokenizer(masked_seq, return_tensors="pt", max_length=max_len, truncation=True).to(DEVICE)
    
    mask_id = tokenizer.mask_token_id
    input_ids = enc["input_ids"][0]
    mask_token_indices = (input_ids == mask_id).nonzero(as_tuple=True)[0]
    
    with torch.no_grad():
        logits = model(**enc).logits[0]
    
    probs_all = logits.softmax(-1)
    vocab = tokenizer.get_vocab()
    
    # Get RNA base token IDs
    base_to_id = {{}}
    for b in "AUGCN":
        if b in vocab:
            base_to_id[b] = vocab[b]
    
    predictions = []
    for i, tok_idx in enumerate(mask_token_indices):
        if i >= len(mask_positions):
            break
        pos = mask_positions[i]
        if pos >= len(orig_seq):
            continue
        orig_base = orig_seq[pos]
        
        probs = probs_all[tok_idx]
        base_probs = [(b, float(probs[bid])) for b, bid in base_to_id.items() if bid < len(probs)]
        base_probs.sort(key=lambda x: x[1], reverse=True)
        
        if not base_probs:
            pred_base = orig_base
            pred_prob = 0.25
        else:
            pred_base, pred_prob = base_probs[0]
        
        predictions.append({{
            "position": pos,
            "original_base": orig_base,
            "predicted_base": pred_base,
            "prob": round(pred_prob, 4),
            "correct": pred_base.upper() == orig_base.upper(),
            "top5": [{{"base": b, "prob": round(p, 4)}} for b, p in base_probs[:5]],
        }})
    
    n_correct = sum(1 for p in predictions if p["correct"])
    total = len(predictions)
    recovery = n_correct / total if total > 0 else 0.0
    
    print(f"  {{rna_name}}/{{strategy}}: {{total}} masks, recovery={{recovery:.1%}}")
    
    results.append({{
        "rna_name": rna_name,
        "description": item["description"],
        "risk_type": item["risk_type"],
        "organism": item["organism"],
        "strategy": strategy,
        "n_masks": len(mask_positions),
        "rnafm_used": True,
        "model": "multimolecule/rnafm",
        "recovery_rate": round(recovery, 4),
        "predictions": predictions,
        "original_seq": orig_seq,
        "functional_sites": item.get("functional_sites", []),
    }})

outpath = os.path.join(RESULTS_DIR, "rnafm_predictions.json")
with open(outpath, "w") as f:
    json.dump(results, f, indent=2)
print(f"Saved {{len(results)}} results to {{outpath}}")
"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.py', delete=False) as tf:
        tf.write(script)
        script_path = tf.name
    try:
        result = subprocess.run(
            [BIOSAFETY_PYTHON, script_path],
            capture_output=False, text=True, timeout=600
        )
        return result.returncode == 0
    finally:
        os.unlink(script_path)

def main():
    print("=== T2.1 RNA-FM Mask & Fill ===")
    success = run_rnafm_subprocess(DATA_FILE, RESULTS_DIR)
    if not success:
        print("RNA-FM subprocess failed")
        sys.exit(1)
    print("RNA-FM inference complete")

if __name__ == "__main__":
    main()
