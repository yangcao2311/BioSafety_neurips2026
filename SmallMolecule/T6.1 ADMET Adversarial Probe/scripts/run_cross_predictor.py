#!/usr/bin/env python3
"""
T6.1 Cross-Predictor Robustness Probe.

Implements MolE-style and DeepTox-style predictors as ADMET-AI alternatives,
trained on Tox21, evaluated on the same 221 ClinTox+T4.x molecules used by ADMET-AI.

- MolE-style:  molfeat ECFP4 (1024-bit) + sklearn LogisticRegression — captures
               canonical molecular embedding + linear classifier, the design
               family of MolE / Mendez-Lucio et al. predictors.
- DeepTox-style: Morgan radius-2 FP (2048-bit) + sklearn MLP — mirrors the
                 Mayr et al. 2016 DeepTox architecture (FP descriptors → DNN).
- Trained on Tox21 SR-p53 endpoint (most aligned with general toxicity).

Outputs:
- results/predictions/cross_predictor_results.json
- Per-predictor false-safe / false-toxic rates
- Pairwise agreement matrix (ADMET-AI vs MolE-style vs DeepTox-style)
"""
import json
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score

BASE = Path(__file__).resolve().parent.parent
PRED_DIR = BASE / "results" / "predictions"
PRED_DIR.mkdir(parents=True, exist_ok=True)

# Tox21 CSV (full) lives in Yang's old workspace (read-only)
TOX21_CSV = Path("/home/ycao95/BioSafety/code/Task_2.1/data/tox21.csv")

ADMET_PRED = PRED_DIR / "admet_predictions.json"

ENDPOINT = "SR-p53"   # general toxicity endpoint, most aligned with ClinTox CT_TOX


def smiles_to_morgan(smi, radius=2, nbits=2048):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits)
    arr = np.zeros((nbits,), dtype=np.int8)
    AllChem.DataStructs.ConvertToNumpyArray(fp, arr)
    return arr


def smiles_to_ecfp4(smi, radius=2, nbits=1024):
    return smiles_to_morgan(smi, radius=radius, nbits=nbits)


def load_tox21_endpoint(endpoint=ENDPOINT):
    """Load Tox21 SMILES and binary labels for one endpoint, drop NA."""
    df = pd.read_csv(TOX21_CSV)
    sub = df[["smiles", endpoint]].dropna()
    sub = sub[sub[endpoint].isin([0, 1])]
    print(f"  Tox21 [{endpoint}]: {len(sub)} compounds, "
          f"{int(sub[endpoint].sum())} positive ({100*sub[endpoint].mean():.1f}%)")
    return sub


def featurize(df, fp_func, label_col):
    feats, labels, idx = [], [], []
    for i, r in df.iterrows():
        f = fp_func(r["smiles"])
        if f is None:
            continue
        feats.append(f)
        labels.append(int(r[label_col]))
        idx.append(i)
    return np.array(feats), np.array(labels)


def train_mole_style(X, y):
    """MolE-style: ECFP4 + LogisticRegression."""
    sc = StandardScaler(with_mean=False)
    Xs = sc.fit_transform(X)
    clf = LogisticRegression(C=1.0, max_iter=500, class_weight="balanced", n_jobs=-1)
    clf.fit(Xs, y)
    return ("MolE-style", sc, clf)


def train_deeptox_style(X, y):
    """DeepTox-style: Morgan-2048 + MLP (small, 3-layer)."""
    sc = StandardScaler(with_mean=False)
    Xs = sc.fit_transform(X)
    clf = MLPClassifier(hidden_layer_sizes=(512, 128), activation="relu",
                        learning_rate_init=1e-3, max_iter=80, batch_size=128,
                        early_stopping=True, validation_fraction=0.1,
                        random_state=42)
    clf.fit(Xs, y)
    return ("DeepTox-style", sc, clf)


def predict_proba(model_tuple, X):
    name, sc, clf = model_tuple
    Xs = sc.transform(X)
    p = clf.predict_proba(Xs)[:, 1]
    return p


def main():
    print("=== T6.1 Cross-Predictor Robustness Probe ===")
    print(f"Endpoint for training: Tox21 {ENDPOINT}")

    # 1) Load existing T6.1 ADMET-AI predictions
    if not ADMET_PRED.exists():
        raise SystemExit(f"Missing {ADMET_PRED}; run run_admet_probe.py first")
    admet = json.load(open(ADMET_PRED))
    rows = admet["all_results"]
    print(f"  Loaded {len(rows)} ADMET-AI predictions")

    # 2) Load Tox21 + featurize
    tox = load_tox21_endpoint(ENDPOINT)
    X_morgan, y = featurize(tox, smiles_to_morgan, ENDPOINT)
    X_ecfp4, _ = featurize(tox, smiles_to_ecfp4, ENDPOINT)
    print(f"  Tox21 features: Morgan-2048 {X_morgan.shape}, ECFP4-1024 {X_ecfp4.shape}")

    # 3) 80/20 split for held-out auroc reporting (sanity-check predictor quality)
    rng = np.random.default_rng(42)
    n = len(y)
    perm = rng.permutation(n)
    split = int(0.8 * n)
    tr, te = perm[:split], perm[split:]

    # MolE-style on ECFP4
    mole_train = train_mole_style(X_ecfp4[tr], y[tr])
    p_te = predict_proba(mole_train, X_ecfp4[te])
    auc_mole = roc_auc_score(y[te], p_te)
    # DeepTox-style on Morgan
    dt_train = train_deeptox_style(X_morgan[tr], y[tr])
    p_te_dt = predict_proba(dt_train, X_morgan[te])
    auc_dt = roc_auc_score(y[te], p_te_dt)
    print(f"  Held-out AUC: MolE-style={auc_mole:.3f}, DeepTox-style={auc_dt:.3f}")

    # Re-train on all data for final inference
    mole_full = train_mole_style(X_ecfp4, y)
    dt_full = train_deeptox_style(X_morgan, y)

    # 4) Apply to T6.1 test set
    print(f"\n  Predicting on {len(rows)} T6.1 molecules...")
    cross_results = []
    n_skip = 0
    for r in rows:
        smi = r["smiles"]
        ec = smiles_to_ecfp4(smi)
        mg = smiles_to_morgan(smi)
        if ec is None or mg is None:
            n_skip += 1
            continue
        p_mole = float(predict_proba(mole_full, ec.reshape(1, -1))[0])
        p_dt = float(predict_proba(dt_full, mg.reshape(1, -1))[0])
        admet_pred_toxic = bool(r.get("predicted_toxic", False))
        admet_clintox = float(r.get("ClinTox", 0))
        known_toxic = bool(r.get("known_toxic", False))

        # Threshold at 0.5 for binary
        mole_pred = p_mole > 0.5
        dt_pred = p_dt > 0.5

        # Adversarial: predicted-safe but known-toxic
        cross_results.append({
            "smiles": smi,
            "source": r.get("source", "unknown"),
            "known_toxic": known_toxic,
            "label": r.get("label", "TOXIC" if known_toxic else "SAFE"),
            # ADMET-AI
            "admet_clintox": round(admet_clintox, 4),
            "admet_pred_toxic": admet_pred_toxic,
            "admet_false_safe": (known_toxic and not admet_pred_toxic),
            # MolE-style
            "mole_score": round(p_mole, 4),
            "mole_pred_toxic": mole_pred,
            "mole_false_safe": (known_toxic and not mole_pred),
            # DeepTox-style
            "deeptox_score": round(p_dt, 4),
            "deeptox_pred_toxic": dt_pred,
            "deeptox_false_safe": (known_toxic and not dt_pred),
        })
    print(f"  Skipped {n_skip} unparsable SMILES")

    # 5) Aggregate stats
    n_total = len(cross_results)
    n_toxic = sum(1 for r in cross_results if r["known_toxic"])
    n_safe = n_total - n_toxic

    def rate(field):
        return sum(1 for r in cross_results if r[field])

    def predictor_rates(prefix):
        n_pred_toxic = sum(1 for r in cross_results if r[f"{prefix}_pred_toxic"])
        # Recall on toxic: TP / (TP + FN)
        tp = sum(1 for r in cross_results if r["known_toxic"] and r[f"{prefix}_pred_toxic"])
        fn = sum(1 for r in cross_results if r["known_toxic"] and not r[f"{prefix}_pred_toxic"])
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0
        # Specificity on safe: TN / (TN + FP)
        tn = sum(1 for r in cross_results if not r["known_toxic"] and not r[f"{prefix}_pred_toxic"])
        fp = sum(1 for r in cross_results if not r["known_toxic"] and r[f"{prefix}_pred_toxic"])
        spec = tn / (tn + fp) if (tn + fp) > 0 else 0
        return {
            "n_pred_toxic": n_pred_toxic,
            "tp": tp, "fp": fp, "tn": tn, "fn": fn,
            "recall_on_toxic": round(recall, 4),
            "specificity_on_safe": round(spec, 4),
            "false_safe_count": fn,
            "false_safe_rate": round(fn / n_toxic, 4) if n_toxic > 0 else 0,
        }

    summary = {
        "n_total": n_total,
        "n_known_toxic": n_toxic,
        "n_known_safe": n_safe,
        "training_endpoint": ENDPOINT,
        "training_set_size": int(len(y)),
        "training_auc_holdout": {
            "MolE-style (ECFP4 + LogReg)": round(auc_mole, 4),
            "DeepTox-style (Morgan-2048 + MLP)": round(auc_dt, 4),
        },
        "predictor_stats": {
            "ADMET-AI": predictor_rates("admet"),
            "MolE-style": predictor_rates("mole"),
            "DeepTox-style": predictor_rates("deeptox"),
        },
    }

    # Pairwise agreement (Cohen's kappa not needed; use raw agreement)
    pairs = [("admet", "mole"), ("admet", "deeptox"), ("mole", "deeptox")]
    pairwise = {}
    for a, b in pairs:
        agree = sum(1 for r in cross_results
                    if r[f"{a}_pred_toxic"] == r[f"{b}_pred_toxic"])
        pairwise[f"{a}_vs_{b}"] = round(agree / n_total, 4)
    summary["pairwise_agreement"] = pairwise

    # All-three-agree on false-safe (indicates systematic blind spot)
    all_three_false_safe = sum(
        1 for r in cross_results
        if r["known_toxic"] and not r["admet_pred_toxic"]
        and not r["mole_pred_toxic"] and not r["deeptox_pred_toxic"]
    )
    summary["all_three_false_safe"] = all_three_false_safe
    summary["all_three_false_safe_rate"] = round(
        all_three_false_safe / n_toxic, 4) if n_toxic > 0 else 0

    print("\n=== Predictor stats ===")
    for k, v in summary["predictor_stats"].items():
        print(f"  {k}: recall_toxic={v['recall_on_toxic']:.3f}, "
              f"spec_safe={v['specificity_on_safe']:.3f}, "
              f"false_safe={v['false_safe_count']}/{n_toxic} ({v['false_safe_rate']*100:.1f}%)")
    print("\n=== Pairwise agreement ===")
    for k, v in pairwise.items():
        print(f"  {k}: {v*100:.1f}%")
    print(f"\nAll-three false_safe: {all_three_false_safe}/{n_toxic} "
          f"({summary['all_three_false_safe_rate']*100:.1f}%) — systematic blind spot")

    out = {
        "task": "T6.1 Cross-Predictor Robustness Extension",
        "summary": summary,
        "results": cross_results,
    }
    out_path = PRED_DIR / "cross_predictor_results.json"
    json.dump(out, open(out_path, "w"), indent=2)
    print(f"\nSaved → {out_path}")


if __name__ == "__main__":
    main()
