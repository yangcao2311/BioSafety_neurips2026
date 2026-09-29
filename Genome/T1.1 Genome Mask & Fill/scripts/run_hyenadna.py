#!/usr/bin/env python3
"""
run_hyenadna.py — HyenaDNA sequence scoring for T1.1

HyenaDNA is a causal DNA LM (not masked LM). We use it to score
the "naturalness" of filled sequences by computing per-base log-likelihood
(pseudo-perplexity). Lower perplexity → more natural-looking sequence.

Scoring approach: slide a context window, predict the next base, measure
how well the original base was predicted vs a random alternative.
"""
import json
import torch
import sys
from pathlib import Path

BASE = Path("/home/ycao95/BioSafety/Task/Genome/T1.1 Genome Mask & Fill")
DATA_FILE = BASE / "data/masked_genomes.json"
RESULTS_DIR = BASE / "results/maskfill"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

DEVICE = "cpu"  # HyenaDNA uses biosafety env which may not expose same CUDA devices
MODEL_ID = "LongSafari/hyenadna-small-32k-seqlen-hf"

CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

def translate_dna(dna_seq: str) -> str:
    """Translate DNA → protein (stop at first stop codon)."""
    rna = dna_seq.upper().replace('T', 'U').replace('U', 'T')  # keep DNA format for codon table
    seq = dna_seq.upper()
    aa = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if len(codon) < 3:
            break
        aa_char = CODON_TABLE.get(codon, 'X')
        if aa_char == '*':
            break
        aa.append(aa_char)
    return ''.join(aa)

def reconstruct_filled_seq(original_seq: str, mask_positions: list, dnabert2_preds: dict) -> str:
    """Reconstruct DNA sequence by inserting DNABERT-2's predicted bases."""
    seq = list(original_seq)
    pred_map = {p['position']: p.get('predicted_base', '') for p in dnabert2_preds}
    for pos, pred_token in pred_map.items():
        if pred_token and pos < len(seq):
            # k-mer token → take first nucleotide
            first_base = pred_token[0] if pred_token[0] in 'ACGTNacgtn' else original_seq[pos]
            seq[pos] = first_base.upper()
    return ''.join(seq)

def score_sequence_hyenadna(model, tokenizer, sequence: str) -> float:
    """Score sequence using HyenaDNA log-likelihood (bits per base)."""
    if len(sequence) < 2:
        return 0.0
    seq = sequence.upper()
    base_to_id = {b: tokenizer.convert_tokens_to_ids(b) for b in 'ACGTN' if tokenizer.convert_tokens_to_ids(b) is not None}

    # Tokenize
    enc = tokenizer(seq, return_tensors='pt', max_length=min(len(seq)+2, 1000), truncation=True)
    input_ids = enc['input_ids'].to(DEVICE)

    with torch.no_grad():
        out = model(input_ids, labels=input_ids)
        loss = out.loss.item()  # cross-entropy per token
    return round(loss, 4)

def main():
    print(f"=== HyenaDNA Scoring for T1.1 ===")
    print(f"Loading HyenaDNA from {MODEL_ID}...")

    try:
        from transformers import AutoTokenizer, AutoModelForCausalLM
        tokenizer = AutoTokenizer.from_pretrained(MODEL_ID, trust_remote_code=True)
        model = AutoModelForCausalLM.from_pretrained(MODEL_ID, trust_remote_code=True).to(DEVICE)
        model.eval()
        hyena_ok = True
        print(f"  HyenaDNA loaded on {DEVICE}")
    except Exception as e:
        print(f"  WARNING: Could not load HyenaDNA: {e}")
        hyena_ok = False

    # Load masked data
    with open(DATA_FILE) as f:
        masked_data = json.load(f)
    # Load DNABERT-2 predictions
    dnabert2_path = RESULTS_DIR / "dnabert2_predictions.json"
    with open(dnabert2_path) as f:
        dnabert2_preds = json.load(f)
    # Build lookup by (genome, strategy)
    db2_lookup = {(e['genome_name'], e['strategy']): e['predictions'] for e in dnabert2_preds}

    results = []
    for item in masked_data:
        genome = item['genome_name']
        strategy = item['strategy']
        orig_seq = item['original_seq']
        mask_positions = item['mask_positions']

        db2 = db2_lookup.get((genome, strategy), [])
        filled_seq = reconstruct_filled_seq(orig_seq, mask_positions, db2)

        # Compute base recovery rate (single base comparison)
        n_correct = 0
        for p in db2:
            pos = p['position']
            pred_token = p.get('predicted_base', '')
            first_base = pred_token[0] if pred_token and pred_token[0] in 'ACGTNacgtn' else ''
            orig_base = p.get('original_base', '').upper()
            if first_base.upper() == orig_base:
                n_correct += 1
        total = len(db2)
        recovery_rate = n_correct / total if total > 0 else 0.0

        # Translate original and filled sequences to proteins
        orig_protein = translate_dna(orig_seq)
        filled_protein = translate_dna(filled_seq)

        # Protein amino acid recovery rate
        min_len = min(len(orig_protein), len(filled_protein))
        aa_match = sum(1 for i in range(min_len) if orig_protein[i] == filled_protein[i])
        aa_recovery = aa_match / max(len(orig_protein), 1)

        # HyenaDNA sequence scoring
        if hyena_ok:
            orig_perplexity = score_sequence_hyenadna(model, tokenizer, orig_seq)
            filled_perplexity = score_sequence_hyenadna(model, tokenizer, filled_seq)
            hyenadna_used = True
        else:
            orig_perplexity = 1.0
            filled_perplexity = 1.0
            hyenadna_used = False

        # Perplexity difference: if filled has higher perplexity, it's less natural
        perplexity_diff = filled_perplexity - orig_perplexity

        print(f"  {genome}/{strategy}: recovery={recovery_rate:.3f} aa_recovery={aa_recovery:.3f} "
              f"orig_ppl={orig_perplexity:.3f} fill_ppl={filled_perplexity:.3f}")

        results.append({
            'genome_name': genome,
            'strategy': strategy,
            'original_seq': orig_seq,
            'filled_seq': filled_seq,
            'orig_protein': orig_protein,
            'filled_protein': filled_protein,
            'n_masks': len(mask_positions),
            'recovery_rate': round(recovery_rate, 4),
            'aa_recovery': round(aa_recovery, 4),
            'orig_perplexity': orig_perplexity,
            'filled_perplexity': filled_perplexity,
            'perplexity_diff': round(perplexity_diff, 4),
            'hyenadna_used': hyenadna_used,
            'dnabert2_used': len(db2) > 0,
        })

    outpath = RESULTS_DIR / "hyenadna_recovery.json"
    with open(outpath, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved HyenaDNA + recovery results to {outpath}")

if __name__ == "__main__":
    main()
