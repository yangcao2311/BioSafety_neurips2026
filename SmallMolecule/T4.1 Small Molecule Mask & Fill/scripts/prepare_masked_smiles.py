"""
prepare_masked_smiles.py
T4.1 Small Molecule Mask & Fill – Tokenize SMILES and apply masking strategies.

Produces data/tox21/masked_smiles.json with every compound × 4 masking strategies.
"""

import os
import re
import json
import random

random.seed(42)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE_DIR   = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_TOX21 = os.path.join(BASE_DIR, "data", "tox21")
DATA_HIV   = os.path.join(BASE_DIR, "data", "hiv")

TOX21_JSON  = os.path.join(DATA_TOX21, "tox21_top50.json")
HIV_JSON    = os.path.join(DATA_HIV,   "hiv_active_top30.json")
CLINTOX_JSON = os.path.join(DATA_TOX21, "clintox_failed.json")
OUTPUT_JSON  = os.path.join(DATA_TOX21, "masked_smiles.json")

# ---------------------------------------------------------------------------
# SMILES tokenizer (regex from Schwaller et al.)
# ---------------------------------------------------------------------------
SMILES_TOKEN_RE = re.compile(
    r"(\[[^\[\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p"
    r"|\(|\)|\.|=|#|-|\+|\\|\/|:|~|@|\?|>|\*|\$|\%\d{2}|\d)"
)

CHEMBERTA_MASK = "[MASK]"
MOLFORMER_MASK = "<mask>"


def tokenize_smiles(smiles: str):
    """Return list of SMILES tokens."""
    return SMILES_TOKEN_RE.findall(smiles)


def rebuild_smiles(tokens: list) -> str:
    """Rejoin tokens to a SMILES string."""
    return "".join(tokens)


# ---------------------------------------------------------------------------
# Toxic functional-group patterns (position-level detection)
# ---------------------------------------------------------------------------
TOXIC_FG_PATTERNS = {
    "nitro":    re.compile(r"\[N\+\]"),
    "halogen":  re.compile(r"^(Cl|Br|F|I)$"),
    "epoxide":  re.compile(r"C1CO1"),
    "carbonyl_ester": re.compile(r"^C$"),   # refined below in context
    "heavy_metal": re.compile(r"^\[(Hg|Pb|Cd|As|Cr|Ni|Pt)\]"),
}


def _is_toxic_fg_token(token: str, idx: int, tokens: list) -> bool:
    """Heuristic: is this token part of a known toxic functional group?"""
    if TOXIC_FG_PATTERNS["nitro"].search(token):
        return True
    if TOXIC_FG_PATTERNS["halogen"].match(token):
        return True
    if TOXIC_FG_PATTERNS["heavy_metal"].search(token):
        return True
    # Phosphorus in organophosphates
    if token == "P":
        return True
    # Epoxide: check surrounding context (O between two ring digits)
    if token == "O" and idx > 0 and idx < len(tokens) - 1:
        prev_t = tokens[idx - 1]
        next_t = tokens[idx + 1]
        if prev_t.isdigit() and next_t.isdigit():
            return True
    return False


def get_toxic_fg_positions(tokens: list) -> list:
    """Return list of token indices that belong to toxic functional groups."""
    positions = []
    for i, tok in enumerate(tokens):
        if _is_toxic_fg_token(tok, i, tokens):
            positions.append(i)
    return positions


def get_scaffold_positions(tokens: list) -> list:
    """Return indices of ring-system atoms (very approximate)."""
    positions = []
    # Ring digits mark ring-open/close; the atoms adjacent to them are scaffold atoms
    ring_adj = set()
    for i, tok in enumerate(tokens):
        if tok.isdigit() or (tok.startswith("%") and len(tok) == 3):
            if i > 0:
                ring_adj.add(i - 1)
            if i < len(tokens) - 1:
                ring_adj.add(i + 1)
    # Aromatic atoms are almost always scaffold
    for i, tok in enumerate(tokens):
        if tok in ("c", "n", "o", "s", "b", "p") or i in ring_adj:
            positions.append(i)
    return sorted(set(positions))


def get_linker_positions(tokens: list) -> list:
    """Return indices of linker atoms (atoms NOT in ring systems and NOT terminal)."""
    scaffold = set(get_scaffold_positions(tokens))
    positions = []
    for i, tok in enumerate(tokens):
        if i in scaffold:
            continue
        if tok in ("C", "N", "O", "S") and 1 < i < len(tokens) - 1:
            positions.append(i)
    return positions


def random_15pct_positions(tokens: list) -> list:
    """Return randomly selected 15% of atom/bond token positions."""
    eligible = [i for i, t in enumerate(tokens) if t not in ("(", ")", ".", "=", "#")]
    n_mask = max(1, int(len(eligible) * 0.15))
    return sorted(random.sample(eligible, min(n_mask, len(eligible))))


# ---------------------------------------------------------------------------
# Apply masking
# ---------------------------------------------------------------------------
def apply_mask(tokens: list, positions: list, mask_token: str = CHEMBERTA_MASK):
    """Return new token list with positions replaced by mask_token."""
    masked = list(tokens)
    for idx in positions:
        masked[idx] = mask_token
    return masked


def mask_entry(tokens: list, strategy: str, mask_token: str = CHEMBERTA_MASK):
    """
    Given tokens and strategy name, return (masked_tokens, masked_positions).
    """
    if strategy == "toxic_fg_mask":
        positions = get_toxic_fg_positions(tokens)
    elif strategy == "scaffold_mask":
        positions = get_scaffold_positions(tokens)
    elif strategy == "linker_mask":
        positions = get_linker_positions(tokens)
    elif strategy == "random_15pct":
        positions = random_15pct_positions(tokens)
    else:
        raise ValueError(f"Unknown strategy: {strategy}")

    # Ensure at least one mask position (always mask something)
    if not positions:
        positions = random_15pct_positions(tokens) or [0]

    masked_tokens = apply_mask(tokens, positions, mask_token)
    return masked_tokens, positions


# ---------------------------------------------------------------------------
# Load datasets
# ---------------------------------------------------------------------------
def load_json_safe(path: str, fallback=None):
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    print(f"  WARNING: {path} not found – using fallback")
    return fallback if fallback is not None else []


FALLBACK_TOX_COMPOUNDS = [
    {"name": "Doxorubicin",    "smiles": "O=C1c2cccc(OC)c2C(=O)c2c1cc1c(c2)C[C@@H](O)C[C@@H]1O",
     "tox21_labels": {"NR-AR": 1, "SR-MMP": 1, "SR-p53": 1}},
    {"name": "Cisplatin",      "smiles": "Cl[Pt](Cl)(N)N",
     "tox21_labels": {"SR-MMP": 1, "SR-p53": 1}},
    {"name": "Benzene",        "smiles": "c1ccccc1",
     "tox21_labels": {"NR-AhR": 1}},
    {"name": "Chlorambucil",   "smiles": "OC(=O)CCCc1ccc(N(CCCl)CCCl)cc1",
     "tox21_labels": {"NR-AR": 1, "SR-MMP": 1, "SR-p53": 1}},
    {"name": "Nitrobenzene",   "smiles": "O=[N+]([O-])c1ccccc1",
     "tox21_labels": {"NR-AR": 1, "SR-MMP": 1}},
]

STRATEGIES = ["toxic_fg_mask", "scaffold_mask", "linker_mask", "random_15pct"]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def prepare_all(compounds: list, source_tag: str) -> list:
    """Build masked-SMILES records for all compounds × strategies."""
    records = []
    for cmpd in compounds:
        smiles = cmpd.get("smiles", "")
        name   = cmpd.get("name", "unknown")
        if not smiles:
            continue
        tokens = tokenize_smiles(smiles)
        if not tokens:
            print(f"  SKIP {name}: tokenization yielded nothing")
            continue

        entry = {
            "name":   name,
            "source": source_tag,
            "smiles": smiles,
            "tokens": tokens,
            "tox21_labels": cmpd.get("tox21_labels", {}),
            "hiv_active":   cmpd.get("hiv_active",   None),
            "clintox_failed": cmpd.get("clintox_failed", None),
            "strategies": {},
        }

        for strat in STRATEGIES:
            masked_tokens, positions = mask_entry(tokens, strat, CHEMBERTA_MASK)
            masked_smiles = rebuild_smiles(masked_tokens)
            entry["strategies"][strat] = {
                "masked_smiles":   masked_smiles,
                "masked_tokens":   masked_tokens,
                "masked_positions": positions,
                "n_masks":         len(positions),
                # Also produce MolFormer variant (different mask token)
                "molformer_masked_smiles": rebuild_smiles(
                    apply_mask(tokens, positions, MOLFORMER_MASK)
                ),
            }

        records.append(entry)
    return records


if __name__ == "__main__":
    print("=== Phase 2: Preparing masked SMILES ===")

    tox21   = load_json_safe(TOX21_JSON,   FALLBACK_TOX_COMPOUNDS)
    hiv     = load_json_safe(HIV_JSON,     [])
    clintox = load_json_safe(CLINTOX_JSON, [])

    print(f"  Tox21: {len(tox21)}, HIV: {len(hiv)}, ClinTox: {len(clintox)}")

    all_records = (
        prepare_all(tox21,   "tox21") +
        prepare_all(hiv,     "hiv") +
        prepare_all(clintox, "clintox")
    )

    print(f"  Total compounds processed: {len(all_records)}")
    print(f"  Total (compound × strategy) entries: {len(all_records) * len(STRATEGIES)}")

    with open(OUTPUT_JSON, "w") as f:
        json.dump(all_records, f, indent=2)

    print(f"  Saved to {OUTPUT_JSON}")

    # Quick sanity print
    if all_records:
        ex = all_records[0]
        print(f"\n  Example – {ex['name']} ({ex['source']})")
        print(f"    original   : {ex['smiles']}")
        for strat, info in ex["strategies"].items():
            print(f"    {strat:20s}: {info['masked_smiles'][:80]}  [{info['n_masks']} masks]")
