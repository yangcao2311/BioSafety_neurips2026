"""
run_dnabert2.py
Run DNABERT-2 on T1.1 Masked Genome dataset.
"""
import json
import torch
from pathlib import Path
from transformers import BertForMaskedLM, AutoTokenizer

DEVICE = "cuda:0" if torch.cuda.is_available() else "cpu"
BASE = Path("/home/ycao95/BioSafety/Task/Genome/T1.1 Genome Mask & Fill")
DATA_FILE = BASE / "data/masked_genomes.json"
RESULTS_DIR = BASE / "results/maskfill"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

MODEL_ID = "zhihan1996/DNABERT-2-117M"

def main():
    print(f"Loading {MODEL_ID} on {DEVICE}...")
    tokenizer = AutoTokenizer.from_pretrained(MODEL_ID, trust_remote_code=True)
    model = BertForMaskedLM.from_pretrained(MODEL_ID).to(DEVICE)
    model.eval()
    
    with open(DATA_FILE) as f:
        data = json.load(f)
        
    print(f"Loaded {len(data)} masked sequences.")
    
    results = []
    
    for item in data:
        genome_name = item["genome_name"]
        strategy = item["strategy"]
        original_seq = item["original_seq"]
        mask_positions = item["mask_positions"]
        
        # Build masked sequence
        # DNABERT-2 uses [MASK] token. Since it uses BPE, masking arbitrary single bases can disrupt tokenization.
        # A common approach for single-base masking in BPE DNA models is to just replace the base with N, 
        # or force space-separated chars if the model supports it. DNABERT-2 tokenizer handles BPE.
        # If we insert [MASK], the tokenizer will treat it as a single token.
        
        # However, to be precise, if we replace contiguous blocks with [MASK], it works well.
        # For random 15%, we have many single base masks.
        
        masked_seq = list(original_seq)
        for pos in mask_positions:
            masked_seq[pos] = "[MASK]"
            
        # Join into string. Note: DNABERT-2 doesn't need spaces between characters
        masked_seq_str = "".join(masked_seq)
        # Because we put "[MASK]" directly in the string, we might need to be careful with BPE.
        # Actually, let's tokenize it and see if [MASK] is preserved.
        
        inputs = tokenizer(masked_seq_str, return_tensors="pt", max_length=512, truncation=True).to(DEVICE)
        
        with torch.no_grad():
            outputs = model(**inputs)
            
        logits = outputs.logits[0]
        input_ids = inputs["input_ids"][0]
        
        mask_token_id = tokenizer.mask_token_id or tokenizer.convert_tokens_to_ids("[MASK]")
        
        mask_indices = (input_ids == mask_token_id).nonzero(as_tuple=True)[0]
        
        predictions = []
        for i, tok_idx in enumerate(mask_indices):
            if i >= len(mask_positions):
                break # safety
            
            orig_pos = mask_positions[i]
            orig_base = original_seq[orig_pos]
            
            top5 = torch.topk(logits[tok_idx], 5)
            top5_probs = top5.values.softmax(-1).tolist()
            top5_tokens = tokenizer.convert_ids_to_tokens(top5.indices.tolist())
            
            predictions.append({
                "position": orig_pos,
                "original_base": orig_base,
                "predicted_base": top5_tokens[0],
                "prob": round(top5_probs[0], 4),
                "top5": [{"base": t, "prob": round(p, 4)} for t, p in zip(top5_tokens, top5_probs)]
            })
            
        results.append({
            "genome_name": genome_name,
            "strategy": strategy,
            "model": "DNABERT-2",
            "n_masks": len(mask_positions),
            "predictions": predictions
        })
        
        print(f"Processed {genome_name} ({strategy}) - found {len(mask_indices)} mask tokens.")
        
    out_file = RESULTS_DIR / "dnabert2_predictions.json"
    with open(out_file, "w") as f:
        json.dump(results, f, indent=2)
    print(f"Saved results to {out_file}")

if __name__ == "__main__":
    main()
