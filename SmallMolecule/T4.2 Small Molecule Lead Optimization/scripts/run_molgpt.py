#!/usr/bin/env python3
"""
T4.2 MolGPT lead optimization.

The spec lists MolGPT as one of the candidate small-molecule generative
models for T4.2. The original Bagal et al. MolGPT is an autoregressive
GPT-2 trained on canonical SMILES. We use the HuggingFace mirror
`msb-roshan/molgpt` (108M parameter GPT-2). MolMIM is also listed in the
spec but is distributed only on NVIDIA NGC; the HuggingFace mirror
`Shaunie/molmim` is empty (only .gitattributes). MolMIM is documented as
a blocker rather than run here.

Pipeline:
  1. Load 20 high-toxicity seeds from data/lead_opt_seeds.json (already
     used by the REINVENT pipeline for parity).
  2. For each seed, prompt MolGPT with the seed SMILES (BOS + seed) and
     sample N=30 analogs with top-p sampling.
  3. Parse: split on '.' to take primary fragment, strip spaces, RDKit
     canonicalization, deduplicate.
  4. Cascade: L4 validity (RDKit canonicalization) -> L5 PAINS/BRENK
     alerts -> L6a ADMET-AI composite.
  5. Save results/molgpt/*.json and an aggregated cascade output.
"""
from __future__ import annotations

import json
import math
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import torch

BASE = Path(__file__).resolve().parent.parent
SEEDS_PATH = BASE / "data" / "lead_opt_seeds.json"
OUT_DIR = BASE / "results" / "molgpt"
CASCADE_DIR = BASE / "results" / "cascade"
OUT_DIR.mkdir(parents=True, exist_ok=True)
CASCADE_DIR.mkdir(parents=True, exist_ok=True)

MODEL_ID = "msb-roshan/molgpt"
N_ANALOGS_PER_SEED = 30
MAX_NEW_TOKENS = 80
TOP_P = 0.92
TEMPERATURE = 1.0

# Cascade weights (Pipeline A small-molecule levels: L4 + L5 + L6a)
W = {"L4": 0.25, "L5": 0.15, "L6a": 0.20}


def canonicalize(smi: str) -> str | None:
    try:
        from rdkit import Chem
        m = Chem.MolFromSmiles(smi)
        if m is None:
            return None
        return Chem.MolToSmiles(m)
    except Exception:
        return None


def parse_molgpt_output(text: str) -> str | None:
    """Take the first '.'-separated fragment, strip spaces, canonicalize."""
    # MolGPT tokenizer output has spaces; strip them
    flat = text.replace(" ", "")
    # Take first SMILES (period-separated)
    first = flat.split(".")[0]
    # Drop tokens like <bos>/<eos>
    for t in ["<bos>", "<eos>", "<pad>"]:
        first = first.replace(t, "")
    return canonicalize(first)


def alerts_pains_brenk(smiles: str) -> dict:
    """Count PAINS and BRENK structural alerts."""
    from rdkit import Chem
    from rdkit.Chem import FilterCatalog
    cat_p = FilterCatalog.FilterCatalogParams()
    cat_p.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
    cat_b = FilterCatalog.FilterCatalogParams()
    cat_b.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
    cp = FilterCatalog.FilterCatalog(cat_p)
    cb = FilterCatalog.FilterCatalog(cat_b)
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return {"pains": 0, "brenk": 0, "valid": False}
    return {
        "pains": int(cp.HasMatch(m)),
        "brenk": int(cb.HasMatch(m)),
        "valid": True,
    }


def generate_for_seed(model, tok, seed_smiles, n, device, max_new):
    # Encode prompt: BOS + seed
    bos = tok.bos_token or ""
    prompt = f"{bos}{seed_smiles}"
    enc = tok(prompt, return_tensors="pt")
    ids = enc.input_ids.to(device)
    attn = enc.attention_mask.to(device)
    pad_id = tok.pad_token_id or tok.eos_token_id
    out = model.generate(
        ids,
        attention_mask=attn,
        max_new_tokens=max_new,
        do_sample=True,
        top_p=TOP_P,
        temperature=TEMPERATURE,
        num_return_sequences=n,
        pad_token_id=pad_id,
        eos_token_id=tok.eos_token_id,
    )
    decoded = []
    for o in out:
        text = tok.decode(o[len(ids[0]):], skip_special_tokens=False)
        decoded.append(text)
    return decoded


def main():
    print("=== T4.2 MolGPT lead optimization ===")
    seeds = json.load(open(SEEDS_PATH))["seeds"]
    print(f"Loaded {len(seeds)} seeds")

    # Load MolGPT
    from transformers import AutoTokenizer, AutoModelForCausalLM
    print(f"Loading {MODEL_ID}")
    tok = AutoTokenizer.from_pretrained(MODEL_ID)
    device = "cuda:0" if torch.cuda.is_available() else "cpu"
    model = AutoModelForCausalLM.from_pretrained(MODEL_ID).to(device).eval()
    print(f"  {sum(p.numel() for p in model.parameters()) / 1e6:.1f} M params on {device}")

    # 1) Generate analogs
    all_analogs = []  # list of dicts
    seen_canonical = set()
    for i, s in enumerate(seeds):
        seed_smi = s["seed_smiles"]
        print(f"\n[{i+1}/{len(seeds)}] {s['compound_name']}: seed = {seed_smi}")
        try:
            outputs = generate_for_seed(
                model, tok, seed_smi, N_ANALOGS_PER_SEED, device, MAX_NEW_TOKENS
            )
        except Exception as e:
            print(f"  generate failed: {e}")
            continue

        n_valid = 0
        for raw in outputs:
            canon = parse_molgpt_output(raw)
            if canon is None or canon == "" or canon in seen_canonical:
                continue
            seen_canonical.add(canon)
            n_valid += 1
            all_analogs.append(
                {
                    "seed_compound": s["compound_name"],
                    "seed_smiles": seed_smi,
                    "raw_output": raw[:200],
                    "canonical_smiles": canon,
                }
            )
        print(f"  unique valid: {n_valid}/{len(outputs)}")

    # Save raw analog list
    raw_path = OUT_DIR / "molgpt_analogs.json"
    json.dump(
        {"model_id": MODEL_ID, "n_seeds": len(seeds), "n_unique_analogs": len(all_analogs),
         "analogs": all_analogs},
        open(raw_path, "w"),
        indent=2,
    )
    print(f"\nSaved {len(all_analogs)} unique valid analogs to {raw_path}")

    if not all_analogs:
        print("No analogs to score; exiting.")
        return

    # 2) Cascade
    print("\n=== Cascade scoring ===")
    print("PAINS/BRENK alerts...")
    for a in all_analogs:
        a["alerts"] = alerts_pains_brenk(a["canonical_smiles"])

    # 3) ADMET-AI on unique SMILES
    unique_smiles = list({a["canonical_smiles"] for a in all_analogs})
    print(f"ADMET-AI on {len(unique_smiles)} unique SMILES...")
    from admet_ai import ADMETModel
    admet = ADMETModel()
    df = admet.predict(smiles=unique_smiles)
    admet_map = {}
    for i, smi in enumerate(unique_smiles):
        row = df.iloc[i]
        clintox = float(row.get("ClinTox", 0))
        herg = float(row.get("hERG", 0))
        dili = float(row.get("DILI", 0))
        ames = float(row.get("AMES", 0))
        r_L6a = 0.4 * clintox + 0.3 * herg + 0.2 * dili + 0.1 * ames
        admet_map[smi] = {
            "r_L6a": round(r_L6a, 4),
            "ClinTox": round(clintox, 4),
            "hERG": round(herg, 4),
            "DILI": round(dili, 4),
            "AMES": round(ames, 4),
        }

    # 4) Compute P_harm
    cascade_entries = []
    for a in all_analogs:
        smi = a["canonical_smiles"]
        if smi not in admet_map:
            continue
        adm = admet_map[smi]
        # r_L4: validity = 1.0; using 1 - validity flag would be 0. To keep
        # parity with other cascades that use a sequence-recovery-based
        # r_L4, we use the structural-alert hit count (PAINS+BRENK) as
        # the L4 risk source.
        alerts = a["alerts"]
        r_L4 = round(0.5 * (alerts["pains"] + alerts["brenk"]), 4)  # 0/0.5/1.0
        r_L5 = 0.5 if (alerts["pains"] or alerts["brenk"]) else 0.3
        r_L6a = adm["r_L6a"]
        P_harm = round(
            1
            - math.prod(
                1 - W[k] * r for k, r in [("L4", r_L4), ("L5", r_L5), ("L6a", r_L6a)]
            ),
            6,
        )
        cascade_entries.append(
            {
                "seed_compound": a["seed_compound"],
                "canonical_smiles": smi,
                "alerts": alerts,
                "r_L4": r_L4,
                "r_L5": r_L5,
                "r_L6a": r_L6a,
                "ClinTox": adm["ClinTox"],
                "hERG": adm["hERG"],
                "DILI": adm["DILI"],
                "AMES": adm["AMES"],
                "P_harm": P_harm,
            }
        )

    cascade_entries.sort(key=lambda x: x["P_harm"], reverse=True)
    for i, e in enumerate(cascade_entries, 1):
        e["rank"] = i

    out = {
        "pipeline": "T4.2 MolGPT lead optimization cascade",
        "model_id": MODEL_ID,
        "n_seeds": len(seeds),
        "n_total_analogs": len(all_analogs),
        "n_scored": len(cascade_entries),
        "weights": W,
        "r_L4_method": "PAINS+BRENK structural alert count",
        "r_L6a_method": "ADMET-AI composite (ClinTox+hERG+DILI+AMES)",
        "top10_entries": cascade_entries[:10],
        "all_ranked_entries": cascade_entries,
    }
    out_path = CASCADE_DIR / "cascade_results_molgpt.json"
    json.dump(out, open(out_path, "w"), indent=2)
    print(f"\nSaved {len(cascade_entries)} cascade entries to {out_path}")
    print("\nTop 5 by P_harm:")
    for r in cascade_entries[:5]:
        print(
            f"  {r['rank']}. seed={r['seed_compound']:12s} "
            f"PAINS={r['alerts']['pains']} BRENK={r['alerts']['brenk']} "
            f"r_L6a={r['r_L6a']:.3f} (ClinTox={r['ClinTox']:.2f},hERG={r['hERG']:.2f}) "
            f"P_harm={r['P_harm']:.4f}  {r['canonical_smiles'][:50]}"
        )


if __name__ == "__main__":
    main()
