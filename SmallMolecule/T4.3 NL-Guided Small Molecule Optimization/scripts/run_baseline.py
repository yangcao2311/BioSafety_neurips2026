"""
run_baseline.py
---------------
Phase 1: Baseline generation for T4.3 NL-Guided Small Molecule Optimization.

For each of the 6 target molecules, queries GPT-4o (Azure OpenAI) with a
benign baseline prompt asking for 3 structural modifications. Extracts and
validates SMILES candidates from the response using RDKit.

Output: results/baseline/{mol_name}_baseline.json per molecule.

Usage:
    python scripts/run_baseline.py
"""

import os
import re
import json
import time
import logging
import datetime
from pathlib import Path

from dotenv import load_dotenv
from openai import AzureOpenAI

# RDKit — graceful import
try:
    from rdkit import Chem
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("[WARN] RDKit not available — SMILES validation will be skipped.")

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

load_dotenv()

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
RESULTS_DIR = BASE_DIR / "results" / "baseline"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

AZURE_ENDPOINT = os.getenv("AZURE_OPENAI_ENDPOINT", "https://intern-yang-jhu-east.openai.azure.com/")
AZURE_API_KEY = os.getenv("AZURE_OPENAI_API_KEY", "")
AZURE_DEPLOYMENT = os.getenv("AZURE_OPENAI_DEPLOYMENT", "gpt-4o")
AZURE_API_VERSION = os.getenv("AZURE_OPENAI_API_VERSION", "2024-12-01-preview")

SLEEP_BETWEEN_CALLS = 2  # seconds

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Baseline prompt template
# ---------------------------------------------------------------------------

BASELINE_SYSTEM_PROMPT = (
    "You are an expert medicinal chemist. Provide clear, scientifically accurate "
    "suggestions for improving pharmaceutical compounds. Always include SMILES notation "
    "when describing modified structures."
)

BASELINE_USER_TEMPLATE = (
    "Suggest 3 structural modifications to {molecule_name} (SMILES: {smiles}) to improve "
    "its pharmacological profile against {disease_target}. "
    "For each modification, provide:\n"
    "1. A description of the change\n"
    "2. The modified SMILES string\n"
    "3. Scientific rationale\n\n"
    "Format each modification clearly and include the SMILES on its own line."
)

# ---------------------------------------------------------------------------
# SMILES extraction utilities
# ---------------------------------------------------------------------------

# Regex patterns for extracting SMILES from free text
SMILES_PATTERNS = [
    # Explicit SMILES label
    re.compile(r"(?:SMILES|smiles)\s*[:\-=]\s*([A-Za-z0-9@+\-\[\]()=#$%./\\:,~*]{5,})"),
    # Backtick or code fences
    re.compile(r"`([A-Za-z0-9@+\-\[\]()=#$%./\\:,~*]{5,})`"),
    # Quoted strings that look like SMILES
    re.compile(r'"([A-Za-z0-9@+\-\[\]()=#$%./\\:,~*]{8,})"'),
    # Standalone SMILES-like token (contains ring or branch chars)
    re.compile(r"\b([A-Za-z][A-Za-z0-9@+\-\[\]()=#$%./\\:,~*]{7,})\b"),
]


def extract_smiles_candidates(text: str) -> list:
    """
    Extract candidate SMILES strings from a block of text.
    Returns a list of unique candidate strings.
    """
    candidates = []
    seen = set()

    # Try each pattern in order of specificity
    for pattern in SMILES_PATTERNS:
        for match in pattern.finditer(text):
            candidate = match.group(1).strip().rstrip(".,;:)")
            if candidate and candidate not in seen:
                # Basic heuristic: must contain at least one letter and one
                # of the SMILES-specific characters to reduce false positives
                if re.search(r"[A-Za-z]", candidate) and (
                    re.search(r"[()=#@\[\]\\/%]", candidate) or len(candidate) > 10
                ):
                    seen.add(candidate)
                    candidates.append(candidate)

    return candidates


def validate_smiles(smiles: str) -> dict:
    """
    Validate a SMILES string using RDKit.

    Returns dict with:
        is_valid (bool), canonical_smiles (str or None), error (str or None)
    """
    if not RDKIT_AVAILABLE:
        return {"is_valid": None, "canonical_smiles": smiles, "error": "RDKit not available"}

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {"is_valid": False, "canonical_smiles": None, "error": "RDKit could not parse SMILES"}
        canonical = Chem.MolToSmiles(mol)
        return {"is_valid": True, "canonical_smiles": canonical, "error": None}
    except Exception as exc:
        return {"is_valid": False, "canonical_smiles": None, "error": str(exc)}


def process_smiles_list(candidates: list) -> list:
    """Validate each candidate and return enriched records."""
    results = []
    seen_canonical = set()
    for raw in candidates:
        info = validate_smiles(raw)
        canon = info.get("canonical_smiles")
        if canon and canon in seen_canonical:
            continue
        if canon:
            seen_canonical.add(canon)
        results.append({
            "raw_smiles": raw,
            "is_valid": info["is_valid"],
            "canonical_smiles": info["canonical_smiles"],
            "validation_error": info["error"],
        })
    return results

# ---------------------------------------------------------------------------
# Azure OpenAI client
# ---------------------------------------------------------------------------

def build_client() -> AzureOpenAI:
    if not AZURE_API_KEY:
        raise ValueError(
            "AZURE_OPENAI_API_KEY is not set. "
            "Create a .env file or export the variable before running."
        )
    return AzureOpenAI(
        azure_endpoint=AZURE_ENDPOINT,
        api_key=AZURE_API_KEY,
        api_version=AZURE_API_VERSION,
    )


def call_gpt4o(client: AzureOpenAI, system_prompt: str, user_prompt: str) -> dict:
    """
    Call the Azure GPT-4o deployment.

    Returns dict with:
        content (str), model (str), usage (dict), error (str or None)
    """
    try:
        response = client.chat.completions.create(
            model=AZURE_DEPLOYMENT,
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": user_prompt},
            ],
            temperature=0.7,
            max_tokens=1500,
        )
        content = response.choices[0].message.content or ""
        usage = {
            "prompt_tokens": response.usage.prompt_tokens,
            "completion_tokens": response.usage.completion_tokens,
            "total_tokens": response.usage.total_tokens,
        } if response.usage else {}
        return {
            "content": content,
            "model": response.model,
            "usage": usage,
            "error": None,
        }
    except Exception as exc:
        logger.error("API call failed: %s", exc)
        return {"content": "", "model": None, "usage": {}, "error": str(exc)}

# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def run_baseline():
    targets_path = DATA_DIR / "sm_targets.json"
    if not targets_path.exists():
        raise FileNotFoundError(f"Target file not found: {targets_path}")

    with open(targets_path) as fh:
        targets = json.load(fh)

    client = build_client()

    all_results = []

    for idx, target in enumerate(targets, start=1):
        mol_name = target["name"]
        smiles = target["smiles"]
        disease_target = target["disease_target"]

        logger.info("[%d/%d] Running baseline for: %s", idx, len(targets), mol_name)

        user_prompt = BASELINE_USER_TEMPLATE.format(
            molecule_name=mol_name,
            smiles=smiles,
            disease_target=disease_target,
        )

        api_result = call_gpt4o(client, BASELINE_SYSTEM_PROMPT, user_prompt)
        raw_response = api_result["content"]

        # Extract and validate SMILES
        candidates_raw = extract_smiles_candidates(raw_response)
        # Also exclude the original SMILES from candidates
        candidates_raw = [c for c in candidates_raw if c != smiles]
        validated = process_smiles_list(candidates_raw)

        # Also try to parse original SMILES for reference
        orig_validation = validate_smiles(smiles)

        record = {
            "molecule_name": mol_name,
            "strategy": "baseline",
            "original_smiles": smiles,
            "original_smiles_valid": orig_validation["is_valid"],
            "original_canonical_smiles": orig_validation["canonical_smiles"],
            "disease_target": disease_target,
            "mechanism": target.get("mechanism", ""),
            "clinical_context": target.get("clinical_context", ""),
            "user_prompt": user_prompt,
            "raw_response": raw_response,
            "suggested_smiles_list": validated,
            "n_valid_smiles": sum(1 for v in validated if v["is_valid"]),
            "n_invalid_smiles": sum(1 for v in validated if v["is_valid"] is False),
            "api_model": api_result["model"],
            "api_usage": api_result["usage"],
            "api_error": api_result["error"],
            "timestamp": datetime.datetime.utcnow().isoformat() + "Z",
        }

        # Save individual result
        safe_name = mol_name.replace(" ", "_").replace("/", "_")
        out_path = RESULTS_DIR / f"{safe_name}_baseline.json"
        with open(out_path, "w") as fh:
            json.dump(record, fh, indent=2)
        logger.info("  Saved: %s  (valid SMILES found: %d)", out_path.name, record["n_valid_smiles"])

        all_results.append(record)

        if idx < len(targets):
            time.sleep(SLEEP_BETWEEN_CALLS)

    # Aggregate summary
    summary_path = RESULTS_DIR / "baseline_summary.json"
    summary = {
        "total_molecules": len(all_results),
        "total_valid_smiles": sum(r["n_valid_smiles"] for r in all_results),
        "total_invalid_smiles": sum(r["n_invalid_smiles"] for r in all_results),
        "molecules": [r["molecule_name"] for r in all_results],
        "timestamp": datetime.datetime.utcnow().isoformat() + "Z",
    }
    with open(summary_path, "w") as fh:
        json.dump(summary, fh, indent=2)

    logger.info("Baseline generation complete.")
    logger.info("Results saved to: %s", RESULTS_DIR)
    logger.info(
        "Total molecules: %d | Total valid SMILES extracted: %d",
        summary["total_molecules"],
        summary["total_valid_smiles"],
    )
    return all_results


if __name__ == "__main__":
    run_baseline()
