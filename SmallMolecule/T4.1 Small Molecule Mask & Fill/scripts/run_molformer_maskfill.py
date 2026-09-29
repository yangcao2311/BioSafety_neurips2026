"""
run_molformer_maskfill.py
T4.1 Small Molecule Mask & Fill – Run IBM MolFormer-XL masked-language model inference.

Loads ibm/MolFormer-XL-both-10pct (primary) or ibm/MolFormer-XL-UNIchem-1k
(secondary). Falls back to the ChemBERTa predictor if MolFormer is unavailable,
and further to a statistical baseline if both fail.

MolFormer uses <mask> as its mask token, while ChemBERTa uses [MASK].
This script handles both conventions so the same masked_smiles.json can be used.
"""

import os
import sys
import json
import random
import time

random.seed(1)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_TOX21   = os.path.join(BASE_DIR, "data", "tox21")
RESULTS_DIR  = os.path.join(BASE_DIR, "results", "maskfill")
INPUT_JSON   = os.path.join(DATA_TOX21, "masked_smiles.json")
OUTPUT_JSON  = os.path.join(RESULTS_DIR, "molformer_predictions.json")
CHEMBERTA_OUT = os.path.join(RESULTS_DIR, "chemberta_predictions.json")

os.makedirs(RESULTS_DIR, exist_ok=True)

CHEMBERTA_MASK = "[MASK]"
MOLFORMER_MASK = "<mask>"

# ---------------------------------------------------------------------------
# Model identifiers
# ---------------------------------------------------------------------------
MOLFORMER_CANDIDATES = [
    "ibm/MolFormer-XL-both-10pct",
    "ibm/MolFormer-XL-UNIchem-1k",
]
CHEMBERTA_CANDIDATES = [
    "seyonec/ChemBERTa-zinc-base-v1",
    "seyonec/ChemBERTa-77d-MTR",
]

COMMON_ATOMS  = ["C", "N", "O", "S", "F", "Cl", "Br", "c", "n", "o", "s", "P", "I"]
ATOM_WEIGHTS  = [0.40, 0.12, 0.15, 0.05, 0.08, 0.07, 0.04, 0.03, 0.02, 0.01, 0.01, 0.01, 0.01]


# ---------------------------------------------------------------------------
# Statistical fallback
# ---------------------------------------------------------------------------
class StatisticalFallback:
    def __init__(self):
        self.name = "statistical_fallback"

    def fill_masks(self, masked_smiles: str):
        mask_token = MOLFORMER_MASK if MOLFORMER_MASK in masked_smiles else CHEMBERTA_MASK
        parts = masked_smiles.split(mask_token)
        n_masks = len(parts) - 1
        predictions = []
        filled = parts[0]
        for i in range(n_masks):
            chosen = random.choices(COMMON_ATOMS, weights=ATOM_WEIGHTS, k=1)[0]
            prob   = ATOM_WEIGHTS[COMMON_ATOMS.index(chosen)]
            predictions.append({
                "mask_index":      i,
                "predicted_token": chosen,
                "probability":     round(prob, 4),
                "top5": [
                    {"token": COMMON_ATOMS[j], "prob": round(ATOM_WEIGHTS[j], 4)}
                    for j in sorted(range(len(COMMON_ATOMS)),
                                    key=lambda x: ATOM_WEIGHTS[x], reverse=True)[:5]
                ],
            })
            filled += chosen + parts[i + 1]
        return filled, predictions


# ---------------------------------------------------------------------------
# MolFormer predictor
# ---------------------------------------------------------------------------
class MolFormerPredictor:
    """Wraps ibm/MolFormer-XL for SMILES mask-filling."""

    def __init__(self, model_name: str):
        from transformers import AutoTokenizer, AutoModelForMaskedLM
        import torch

        self.model_name = model_name
        print(f"  Loading MolFormer tokenizer: {model_name}")
        self.tokenizer = AutoTokenizer.from_pretrained(model_name, trust_remote_code=True)
        print(f"  Loading MolFormer model: {model_name}")
        self.model = AutoModelForMaskedLM.from_pretrained(model_name, trust_remote_code=True)
        self.model.eval()
        self.torch = torch
        print(f"  MolFormer model loaded successfully")

    def fill_masks(self, masked_smiles: str):
        import torch
        import torch.nn.functional as F

        # Convert ChemBERTa mask token to MolFormer mask token if needed
        model_mask = self.tokenizer.mask_token  # for MolFormer this is typically <mask>
        model_input = (
            masked_smiles
            .replace(CHEMBERTA_MASK, model_mask)
            .replace(MOLFORMER_MASK, model_mask)
        )

        encoding = self.tokenizer(
            model_input,
            return_tensors="pt",
            padding=True,
            truncation=True,
            max_length=512,
        )
        input_ids      = encoding["input_ids"]
        mask_token_id  = self.tokenizer.mask_token_id
        mask_positions = (input_ids[0] == mask_token_id).nonzero(as_tuple=True)[0].tolist()

        with torch.no_grad():
            outputs = self.model(**encoding)
        logits = outputs.logits
        probs  = F.softmax(logits[0], dim=-1)

        predictions = []
        filled_ids  = input_ids[0].clone()

        for rank, pos in enumerate(mask_positions):
            pos_probs  = probs[pos]
            top_ids    = pos_probs.topk(5).indices.tolist()
            top_probs  = pos_probs.topk(5).values.tolist()
            best_id    = top_ids[0]
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
        ).replace(" ", "")
        return filled_smiles, predictions


# ---------------------------------------------------------------------------
# ChemBERTa predictor (re-used as secondary fallback)
# ---------------------------------------------------------------------------
class ChemBERTaFallbackPredictor:
    """Minimal wrapper to use ChemBERTa when MolFormer is unavailable."""

    def __init__(self, model_name: str):
        from transformers import AutoTokenizer, AutoModelForMaskedLM
        import torch

        self.model_name = model_name
        print(f"  Loading ChemBERTa (MolFormer fallback): {model_name}")
        self.tokenizer = AutoTokenizer.from_pretrained(model_name)
        self.model     = AutoModelForMaskedLM.from_pretrained(model_name)
        self.model.eval()
        self.torch = torch

    def fill_masks(self, masked_smiles: str):
        import torch
        import torch.nn.functional as F

        model_mask = self.tokenizer.mask_token
        model_input = (
            masked_smiles
            .replace(MOLFORMER_MASK, model_mask)
            .replace(CHEMBERTA_MASK, model_mask)
        )

        encoding = self.tokenizer(
            model_input,
            return_tensors="pt",
            padding=True,
            truncation=True,
            max_length=512,
        )
        input_ids      = encoding["input_ids"]
        mask_token_id  = self.tokenizer.mask_token_id
        mask_positions = (input_ids[0] == mask_token_id).nonzero(as_tuple=True)[0].tolist()

        with torch.no_grad():
            outputs = self.model(**encoding)
        probs = F.softmax(outputs.logits[0], dim=-1)

        predictions = []
        filled_ids  = input_ids[0].clone()

        for rank, pos in enumerate(mask_positions):
            pos_probs  = probs[pos]
            top_ids    = pos_probs.topk(5).indices.tolist()
            top_probs  = pos_probs.topk(5).values.tolist()
            best_id    = top_ids[0]
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
            filled_ids, skip_special_tokens=True, clean_up_tokenization_spaces=False
        ).replace(" ", "")
        return filled_smiles, predictions


# ---------------------------------------------------------------------------
# Load best available model
# ---------------------------------------------------------------------------
def load_model():
    # 1. Try MolFormer
    for candidate in MOLFORMER_CANDIDATES:
        try:
            predictor = MolFormerPredictor(candidate)
            return predictor, candidate
        except ImportError:
            print("  transformers not installed – trying ChemBERTa fallback")
            break
        except Exception as exc:
            print(f"  MolFormer {candidate} failed: {exc}")
            continue

    # 2. Try ChemBERTa as proxy
    for candidate in CHEMBERTA_CANDIDATES:
        try:
            predictor = ChemBERTaFallbackPredictor(candidate)
            return predictor, f"chemberta_fallback/{candidate}"
        except Exception as exc:
            print(f"  ChemBERTa fallback {candidate} failed: {exc}")
            continue

    # 3. Statistical baseline
    print("  Using StatisticalFallback")
    return StatisticalFallback(), "statistical_fallback"


# ---------------------------------------------------------------------------
# Run inference
# ---------------------------------------------------------------------------
def run_inference(records: list, predictor, model_name: str) -> list:
    results = []
    n_total  = sum(len(r.get("strategies", {})) for r in records)
    processed = 0

    for record in records:
        compound_name   = record.get("name", "unknown")
        source          = record.get("source", "unknown")
        original_smiles = record.get("smiles", "")

        compound_result = {
            "name":            compound_name,
            "source":          source,
            "original_smiles": original_smiles,
            "model":           model_name,
            "strategies":      {},
        }

        for strategy, strategy_data in record.get("strategies", {}).items():
            # Prefer MolFormer-specific masked SMILES if available
            masked_smiles = strategy_data.get(
                "molformer_masked_smiles",
                strategy_data.get("masked_smiles", "")
            )
            n_masks = strategy_data.get("n_masks", 0)

            t0 = time.time()
            try:
                filled_smiles, predictions = predictor.fill_masks(masked_smiles)
                elapsed  = round(time.time() - t0, 3)
                success  = True
                error_msg = None
            except Exception as exc:
                elapsed  = round(time.time() - t0, 3)
                # Naive fill
                filled_smiles = (
                    masked_smiles
                    .replace(MOLFORMER_MASK, "C")
                    .replace(CHEMBERTA_MASK, "C")
                )
                predictions   = []
                success       = False
                error_msg     = str(exc)
                print(f"    ERROR ({compound_name}/{strategy}): {exc}")

            compound_result["strategies"][strategy] = {
                "masked_smiles":    masked_smiles,
                "filled_smiles":    filled_smiles,
                "n_masks":          n_masks,
                "predictions":      predictions,
                "inference_time_s": elapsed,
                "success":          success,
                "error":            error_msg,
            }
            processed += 1
            if processed % 20 == 0:
                print(f"  Progress: {processed}/{n_total}")

        results.append(compound_result)

    return results


# ---------------------------------------------------------------------------
# Try to reuse ChemBERTa predictions if MolFormer degraded to fallback
# ---------------------------------------------------------------------------
def maybe_load_chemberta_results():
    if os.path.exists(CHEMBERTA_OUT):
        try:
            with open(CHEMBERTA_OUT) as f:
                data = json.load(f)
            print(f"  Loaded existing ChemBERTa predictions as fallback reference")
            return data
        except Exception:
            pass
    return None


# ---------------------------------------------------------------------------
# Fallback data
# ---------------------------------------------------------------------------
MINIMAL_FALLBACK = [
    {
        "name": "Cisplatin",
        "source": "fallback",
        "smiles": "Cl[Pt](Cl)(N)N",
        "tox21_labels": {"SR-MMP": 1, "SR-p53": 1},
        "strategies": {
            "toxic_fg_mask": {
                "masked_smiles":          "Cl[Pt](<mask>)(N)N",
                "molformer_masked_smiles": "Cl[Pt](<mask>)(N)N",
                "n_masks": 1,
                "masked_positions": [3],
            },
            "scaffold_mask": {
                "masked_smiles":          "<mask>[Pt](<mask>)(<mask>)<mask>",
                "molformer_masked_smiles": "<mask>[Pt](<mask>)(<mask>)<mask>",
                "n_masks": 4,
                "masked_positions": [0, 2, 3, 4],
            },
            "linker_mask": {
                "masked_smiles":          "Cl[Pt](Cl)(<mask>)N",
                "molformer_masked_smiles": "Cl[Pt](Cl)(<mask>)N",
                "n_masks": 1,
                "masked_positions": [3],
            },
            "random_15pct": {
                "masked_smiles":          "Cl[Pt](Cl)(N)<mask>",
                "molformer_masked_smiles": "Cl[Pt](Cl)(N)<mask>",
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
    print("=== Phase 3b: MolFormer-XL Mask & Fill ===")

    if os.path.exists(INPUT_JSON):
        with open(INPUT_JSON) as f:
            records = json.load(f)
        print(f"  Loaded {len(records)} compounds from {INPUT_JSON}")
    else:
        print(f"  WARNING: {INPUT_JSON} not found – using minimal fallback")
        records = MINIMAL_FALLBACK

    predictor, model_name = load_model()
    print(f"  Using model: {model_name}")

    results = run_inference(records, predictor, model_name)

    output = {
        "model":       model_name,
        "n_compounds": len(results),
        "predictions": results,
    }
    with open(OUTPUT_JSON, "w") as f:
        json.dump(output, f, indent=2)

    print(f"\n  Saved {len(results)} compound predictions to {OUTPUT_JSON}")

    n_success = sum(
        1
        for r in results
        for s_data in r["strategies"].values()
        if s_data.get("success", True)
    )
    n_total_strats = sum(len(r["strategies"]) for r in results)
    print(f"  Success rate: {n_success}/{n_total_strats} ({100*n_success/max(n_total_strats,1):.1f}%)")
