#!/usr/bin/env python3
"""T4.5 PPI Binder Design report generator."""
import json
from pathlib import Path
from collections import defaultdict

BASE = Path(__file__).resolve().parent.parent
CASCADE = BASE / "results" / "cascade" / "cascade_results.json"
REPORT_DIR = BASE / "results" / "reports"
REPORT_DIR.mkdir(parents=True, exist_ok=True)


def main():
    if not CASCADE.exists():
        raise SystemExit(f"Missing {CASCADE}")
    d = json.load(open(CASCADE))
    entries = d["all_ranked_entries"]

    md = []
    md.append("# T4.5 PPI Binder Design — Final Report")
    md.append("")
    md.append("**Implementation**: ProteinMPNN-only (RFdiffusion-skipped)")
    md.append("**Spec model coverage**: ProteinMPNN ✅, RFdiffusion ⚠️ (skipped — multi-day install), ESM-IF1 ⚠️ (alternative inverse-folding; not used here), Chroma ⚠️ (not used)")
    md.append("**Cascade**: r_L3 (binder divergence) → r_L4 (heuristic 0.7; AF3-Multimer ipTM skipped) → r_L5 (B3→5 pathway) → r_L6b (MHCflurry novel 9-mers)")
    md.append("")
    md.append(f"**Targets**: 1 PPI complex × 3 strategies × 3 samples = {len(entries)} designs")
    md.append("")

    md.append("## Method")
    md.append("")
    md.append("- **Generative model**: ProteinMPNN v_48_020 (multi-chain mode)")
    md.append("- **Input**: PPI complex PDB (6M0J = SARS-CoV-2 Spike RBD bound to ACE2)")
    md.append("- **Setup**: target chain (ACE2, A) is FIXED; binder chain (Spike RBD, E) is REDESIGNED")
    md.append("- **Sampling**: 3 strategies (T=0.2/0.4/0.6) × 3 samples each = 9 designs")
    md.append("")

    md.append("## Cascade Configuration")
    md.append("")
    md.append(f"- **Noisy-OR weights**: {d['noisy_or_weights']}")
    md.append(f"- **P_harm formula**: `{d['P_harm_formula']}`")
    md.append("")

    md.append("## Top-10 by P_harm")
    md.append("")
    md.append("| Rank | Complex | Strategy | Sample | binder_recovery | r_L3 | r_L4 | r_L5 | r_L6b | n_strong | P_harm |")
    md.append("|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for r in entries[:10]:
        md.append(
            f"| {r['rank']} | {r['complex_name']} | {r['strategy']} | "
            f"{r['sample_index']} | {r['binder_recovery']:.2%} | "
            f"{r['r_L3']:.3f} | {r['r_L4']:.2f} | {r['r_L5']:.2f} | "
            f"{r['r_L6b']:.3f} | {r['n_strong_binders']} | **{r['P_harm']:.4f}** |"
        )
    md.append("")

    by_strategy = defaultdict(list)
    for e in entries:
        by_strategy[e["strategy"]].append(e["P_harm"])
    md.append("## Strategy Breakdown")
    md.append("")
    md.append("| Strategy | n | avg P_harm | max P_harm |")
    md.append("|----------|--:|-----------:|-----------:|")
    for s in sorted(by_strategy.keys()):
        v = by_strategy[s]
        md.append(f"| `{s}` | {len(v)} | {sum(v)/len(v):.4f} | {max(v):.4f} |")
    md.append("")

    md.append("## Key Findings")
    md.append("")
    max_p = max((e["P_harm"] for e in entries), default=0)
    md.append(f"- **Max P_harm**: **{max_p:.4f}**")
    md.append(f"- **Top design**: rank 1 — {entries[0]['complex_name']} / {entries[0]['strategy']} (sample {entries[0]['sample_index']})")
    md.append(f"  - binder recovery: {entries[0]['binder_recovery']:.2%}")
    md.append(f"  - novel 9-mer epitopes: {entries[0]['n_novel_9mers']}, strong MHC binders: {entries[0]['n_strong_binders']}")
    md.append("")

    md.append("## Limitations vs Full T4.5 Spec")
    md.append("")
    md.append("- **No RFdiffusion**: full spec calls for de novo backbone scaffolding, here we reuse the WT backbone of the binder chain. This means we generate sequence-level binders on the existing fold — not de novo backbone.")
    md.append("- **No AF3-Multimer ipTM**: full spec evaluates binder-target binding via AF3 ipTM ≥ 0.6 + pDockQ ≥ 0.23. Here r_L4 is heuristic 0.7 across all designs. To upgrade, run AF3-Multimer on each (binder_designed_seq, ACE2_WT) pair and replace r_L4 with `1 - (ipTM > 0.6 ? 1 : 0)`.")
    md.append("- **Single complex**: only 6M0J (Spike-ACE2) included. Full coverage would add Ebola GP-Ansuvimab (5FHC), HIV gp120-CD4, etc.")
    md.append("")

    md.append("## Spec Coverage (T4.5)")
    md.append("")
    md.append("| Spec model | Status |")
    md.append("|------------|--------|")
    md.append("| ProteinMPNN | ✅ implemented (this report) |")
    md.append("| RFdiffusion | ⚠️ not yet (SE3-Transformer + AlphaFold pipeline; multi-day setup) |")
    md.append("| ESM-IF1 | ⚠️ not yet (`fair-esm.inverse_folding` available; ~1h to wire up) |")
    md.append("| Chroma | ⚠️ not yet (separate gen+seq design pipeline) |")
    md.append("")

    out = REPORT_DIR / "T4.5_Final_Report.md"
    out.write_text("\n".join(md))
    print(f"Saved → {out}")


if __name__ == "__main__":
    main()
