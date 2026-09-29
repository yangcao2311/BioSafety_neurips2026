"""
run_chemberta_maskfill.py
T4.1 Small Molecule Mask & Fill – Run ChemBERTa masked-language model inference.

Loads seyonec/ChemBERTa-zinc-base-v1 (primary) or seyonec/ChemBERTa-77d-MTR
(secondary) from HuggingFace and fills [MASK] tokens in SMILES strings.

Falls back to a statistical baseline (most-frequent-atom substitution) when the
HuggingFace model is unavailable or the environment lacks transformers.
"""

import os
import sys
import json
import random
import re
import time

random.seed(0)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_TOX21    = os.path.join(BASE_DIR, "data", "tox21")
RESULTS_DIR   = os.path.join(BASE_DIR, "results", "maskfill")
INPUT_JSON    = os.path.join(DATA_TOX21, "masked_smiles.json")
OUTPUT_JSON   = os.path.join(RESULTS_DIR, "chemberta_predictions.json")

os.makedirs(RESULTS_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# Model identifiers (tried in order)
# ---------------------------------------------------------------------------
MODEL_CANDIDATES = [
    "seyonec/ChemBERTa-zinc-base-v1",
    "seyonec/ChemBERTa-77d-MTR",
    "seyonec/ChemBERTa-zinc-base-v1-finetuned-BBBP",
]

MASK_TOKEN = "[MASK]"

# ---------------------------------------------------------------------------
# Atom vocabulary for statistical fallback
# ---------------------------------------------------------------------------
COMMON_ATOMS = ["C", "N", "O", "S", "F", "Cl", "Br", "c", "n", "o", "s", "P", "I"]
ATOM_WEIGHTS  = [0.40, 0.12, 0.15, 0.05, 0.08, 0.07, 0.04, 0.03, 0.02, 0.01, 0.01, 0.01, 0.01]


# ---------------------------------------------------------------------------
# Statistical fallback predictor
# ---------------------------------------------------------------------------
class StatisticalFallback:
    """
    Replaces each [MASK] with the most probable atom drawn from COMMON_ATOMS
    weighted distribution. Records a fixed pseudo-probability.
    """
    def __init__(self):
        self.name = "statistical_fallback"

    def fill_masks(self, masked_smiles: str):
        """
        Return (filled_smiles, list_of_predictions).
        Each prediction = {'position_index': i, 'predicted_token': t, 'probability': p}.
        """
        tokens = masked_smiles.split(MASK_TOKEN)
        n_masks = len(tokens) - 1
        predictions = []
        filled = tokens[0]
        for i in range(n_masks):
            chosen = random.choices(COMMON_ATOMS, weights=ATOM_WEIGHTS, k=1)[0]
            prob = ATOM_WEIGHTS[COMMON_ATOMS.index(chosen)]
            predictions.append({
                "mask_index": i,
                "predicted_token": chosen,
                "probability": round(prob, 4),
                "top5": [
                    {"token": COMMON_ATOMS[j], "prob": round(ATOM_WEIGHTS[j], 4)}
                    for j in sorted(range(len(COMMON_ATOMS)),
                                    key=lambda x: ATOM_WEIGHTS[x], reverse=True)[:5]
                ],
            })
            filled += chosen + tokens[i + 1]
        return filled, predictions


# ---------------------------------------------------------------------------
# HuggingFace MLM predictor
# ---------------------------------------------------------------------------
class ChemBERTaPredictor:
    """Wraps a HuggingFace AutoModelForMaskedLM for SMILES mask-filling."""

    def __init__(self, model_name: str):
        from transformers import AutoTokenizer, AutoModelForMaskedLM
        import torch

        self.model_name = model_name
        print(f"  Loading tokenizer: {model_name}")
        self.tokenizer = AutoTokenizer.from_pretrained(model_name)
        print(f"  Loading model: {model_name}")
        self.model = AutoModelForMaskedLM.from_pretrained(model_name)
        self.model.eval()
        self.torch = torch
        print(f"  Model loaded successfully ({model_name})")

    def fill_masks(self, masked_smiles: str):
        """
        Fill all [MASK] tokens in masked_smiles.

        Returns (filled_smiles, list_of_predictions).
        """
        import torch

        # The ChemBERTa tokenizer may use different mask token internally
        mask_token = self.tokenizer.mask_token  # typically [MASK]
        # Replace our mask token with the model's mask token
        model_input = masked_smiles.replace(MASK_TOKEN, mask_token)

        encoding = self.tokenizer(
            model_input,
            return_tensors="pt",
            padding=True,
            truncation=True,
            max_length=512,
        )
        input_ids     = encoding["input_ids"]
        attention_mask = encoding["attention_mask"]

        # Find mask positions in the token-id sequence
        mask_token_id = self.tokenizer.mask_token_id
        mask_positions = (input_ids[0] == mask_token_id).nonzero(as_tuple=True)[0].tolist()

        with torch.no_grad():
            outputs = self.model(**encoding)
        logits = outputs.logits  # (1, seq_len, vocab_size)

        import torch.nn.functional as F
        probs = F.softmax(logits[0], dim=-1)  # (seq_len, vocab_size)

        predictions = []
        filled_ids  = input_ids[0].clone()

        for rank, pos in enumerate(mask_positions):
            pos_probs = probs[pos]
            top_ids   = pos_probs.topk(5).indices.tolist()
            top_probs = pos_probs.topk(5).values.tolist()
            best_id   = top_ids[0]
            best_token = self.tokenizer.convert_ids_to_tokens(best_id)
            filled_ids[pos] = best_id
            predictions.append({
                "mask_index":      rank,
                "token_position":  pos,
                "predicted_token": best_token,
                "probability":     round(top_probs[0], 6),
                "top5": [
                    {
                        "token": self.tokenizer.convert_ids_to_tokens(top_ids[j]),
                        "prob":  round(top_probs[j], 6),
                    }
                    for j in range(len(top_ids))
                ],
            })

        filled_smiles = self.tokenizer.decode(
            filled_ids,
            skip_special_tokens=True,
            clean_up_tokenization_spaces=False,
        )
        # Remove extra spaces that may be introduced by some tokenizers
        filled_smiles = filled_smiles.replace(" ", "")
        return filled_smiles, predictions


# ---------------------------------------------------------------------------
# Load model (with fallback)
# ---------------------------------------------------------------------------
def load_model():
    """Try ChemBERTa variants; fall back to statistical model."""
    for candidate in MODEL_CANDIDATES:
        try:
            predictor = ChemBERTaPredictor(candidate)
            return predictor, candidate
        except ImportError:
            print("  transformers not installed – using statistical fallback")
            break
        except Exception as exc:
            print(f"  Failed to load {candidate}: {exc}")
            continue
    print("  Falling back to StatisticalFallback predictor")
    return StatisticalFallback(), "statistical_fallback"


# ---------------------------------------------------------------------------
# Run inference
# ---------------------------------------------------------------------------
def run_inference(records: list, predictor, model_name: str) -> list:
    results = []
    n_total = sum(len(r["strategies"]) for r in records)
    processed = 0

    for record in records:
        compound_name = record["name"]
        source        = record.get("source", "unknown")
        original_smiles = record["smiles"]

        compound_result = {
            "name":            compound_name,
            "source":          source,
            "original_smiles": original_smiles,
            "model":           model_name,
            "strategies":      {},
        }

        for strategy, strategy_data in record["strategies"].items():
            masked_smiles = strategy_data["masked_smiles"]
            n_masks       = strategy_data["n_masks"]

            t0 = time.time()
            try:
                filled_smiles, predictions = predictor.fill_masks(masked_smiles)
                elapsed = round(time.time() - t0, 3)
                success = True
                error_msg = None
            except Exception as exc:
                elapsed = round(time.time() - t0, 3)
                filled_smiles = masked_smiles.replace(MASK_TOKEN, "C")  # naive fill
                predictions   = []
                success       = False
                error_msg     = str(exc)
                print(f"    ERROR ({compound_name}/{strategy}): {exc}")

            compound_result["strategies"][strategy] = {
                "masked_smiles":   masked_smiles,
                "filled_smiles":   filled_smiles,
                "n_masks":         n_masks,
                "predictions":     predictions,
                "inference_time_s": elapsed,
                "success":         success,
                "error":           error_msg,
            }
            processed += 1
            if processed % 20 == 0:
                print(f"  Progress: {processed}/{n_total}")

        results.append(compound_result)

    return results


# ---------------------------------------------------------------------------
# Fallback data (used if masked_smiles.json is missing)
# ---------------------------------------------------------------------------
MINIMAL_FALLBACK = [
    {
        "name": "Chlorambucil",
        "source": "fallback",
        "smiles": "OC(=O)CCCc1ccc(N(CCCl)CCCl)cc1",
        "tox21_labels": {"NR-AR": 1, "SR-MMP": 1, "SR-p53": 1},
        "strategies": {
            "toxic_fg_mask": {
                "masked_smiles": "OC(=O)CCCc1ccc(N(CC[MASK])CC[MASK])cc1",
                "n_masks": 2,
                "masked_positions": [14, 17],
            },
            "scaffold_mask": {
                "masked_smiles": "OC(=O)CCCc1[MASK][MASK][MASK](N(CCCl)CCCl)[MASK][MASK]1",
                "n_masks": 5,
                "masked_positions": [5, 6, 7, 8, 9],
            },
            "linker_mask": {
                "masked_smiles": "OC(=O)[MASK][MASK][MASK]c1ccc(N(CCCl)CCCl)cc1",
                "n_masks": 3,
                "masked_positions": [2, 3, 4],
            },
            "random_15pct": {
                "masked_smiles": "OC(=O)CCC[MASK]1ccc(N(CCCl)CCCl)cc1",
                "n_masks": 1,
                "masked_positions": [4],
            },
        },
    },
]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print("=== Phase 3a: ChemBERTa Mask & Fill ===")

    # Load input
    if os.path.exists(INPUT_JSON):
        with open(INPUT_JSON) as f:
            records = json.load(f)
        print(f"  Loaded {len(records)} compounds from {INPUT_JSON}")
    else:
        print(f"  WARNING: {INPUT_JSON} not found – using minimal fallback data")
        records = MINIMAL_FALLBACK

    # Load model
    predictor, model_name = load_model()
    print(f"  Using model: {model_name}")

    # Run
    results = run_inference(records, predictor, model_name)

    # Save
    output = {
        "model":        model_name,
        "n_compounds":  len(results),
        "predictions":  results,
    }
    with open(OUTPUT_JSON, "w") as f:
        json.dump(output, f, indent=2)

    print(f"\n  Saved {len(results)} compound predictions to {OUTPUT_JSON}")

    # Summary
    n_success = sum(
        1
        for r in results
        for s_data in r["strategies"].values()
        if s_data.get("success", True)
    )
    n_total = sum(len(r["strategies"]) for r in results)
    print(f"  Success rate: {n_success}/{n_total} ({100*n_success/max(n_total,1):.1f}%)")
