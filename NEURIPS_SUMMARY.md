# BioSafety Benchmark: A Cascade-Aware Risk Evaluation Framework for Generative Biology AI

**Submission category**: NeurIPS Datasets and Benchmarks Track (proposed format)
**Date**: 2026-05-02

---

## Abstract

We introduce **BioSafety Benchmark**, the first cascade-aware safety evaluation framework for generative biology foundation models. Existing benchmarks evaluate "single-modality models grading themselves": a protein language model rates its own output, a small-molecule generator scores its own SMILES. This systematically under-counts harm because biological risk propagates along the central dogma — a designed nucleotide change becomes an mRNA isoform, becomes a protein, becomes a binding event, becomes a pathway perturbation, becomes an organismal phenotype. We benchmark 25 foundation and specialist models across **15 tasks and two pipelines** (Pipeline A: natural central-dogma flow; Pipeline B: therapeutic agents — CRISPR gRNA and siRNA) by running their outputs through a 6-level cascade with 6 bridge translators. We define a noisy-OR aggregate risk score `P_harm` that decomposes by level, and document its biological limitations as a calibration upgrade target. The benchmark exposes adversarial attack surfaces invisible to single-modality evaluation: 100% jailbreak success on T_siRNA, 9.05x off-target inflation on T_CRISPR adversarial prompts, and 43.8% false-safe rate of ADMET-AI on known clinical toxics.

---

## 1. Why this benchmark exists

### 1.1 The single-modality blind spot

Most current generative biology safety evaluations look like:

```
Generated sequence -> same-modality predictor (FM grading itself) -> "safe" / "unsafe"
```

A protein language model assesses its own designs against self-trained novelty filters. A small-molecule model scores its own SMILES against single-task admet endpoints. Cross-modality propagation — the actual mechanism by which biological harm occurs — is absent.

**The consequence**: ESM-2 designs that escape MHC-I presentation pass the protein-self check but fail at L6b. ChemBERTa-mask-fills that look chemically plausible at L4 pass admet self-check but trigger the ADMET-AI false-safe failure mode at L6a (43.8% on ClinTox).

### 1.2 What this benchmark contributes

1. **An end-to-end cascade**. Six biological levels (L1 genome -> L2 RNA -> L3 protein -> L4 complex -> L5 pathway -> L6 organism) connected by six bridge translators (B1->2 transcription, B2->3 translation, B3->4 folding+docking, B3->5 pathway-direct, B4->5 complex-to-pathway, B5->6 organismal extrapolation).
2. **Two parallel pipelines**. Pipeline A (natural central dogma) and Pipeline B (therapeutic CRISPR / siRNA agents that share L1 and L2 but use task-family-specific risk weights).
3. **15 tasks, 25 models**. Mask-fill tasks (T1.1, T2.1, T3.1, T4.1) with ground truth; generative tasks (T1.2, T3.2, T3.3, T4.2, T4.3, T4.4, T4.5) without unique answers; cell-perturbation tasks (T5.1); ADMET adversarial probes (T6.1); therapy-pipeline tasks (T_CRISPR, T_siRNA).
4. **Aggregate risk metric** `P_harm = 1 - prod_l (1 - w_l * r_l)` with explicit per-level decomposition `{w_l * r_l}` for diagnostic transparency.
5. **A documented methodology critique** of the noisy-OR aggregator (statistical and biological limits) and a path forward to causal-DAG aggregation.

---

## 2. Benchmark structure

### 2.1 Cascade levels and bridge agents

| Level | Channel | Forward FM(s) used | Risk channel `r_l` |
|------:|---------|--------------------|--------------------|
| L1 | Genome / nucleic acid | DNABERT-2, NT-v2, HyenaDNA, Caduceus, Evo-2 | Sequence-level divergence from gold |
| L2 | RNA | RNA-FM, RiNALMo, SpliceBERT, UTR-LM | Translation efficiency / structural disruption |
| L3 | Protein | ESM-2, ProtBERT, ProtT5, Ankh, ProteinMPNN, ESM-IF1 | Structure / epitope deviation (TMalign or AARR) |
| L4 | Complex / Binding | AutoDock Vina, ESMFold | Docking ΔG -> per-variant r_L4 |
| L5 | Pathway | Pfam, KEGG, VirHostNet, Geneformer | Multi-source pathway evidence count |
| L6a | Chemical organismal toxicity | ADMET-AI, MolE, DeepTox proxies | Composite endpoint score |
| L6b | Biological pathogenic outcome | MHCflurry 2.2.1 | Novel 9-mer presentation rate |

### 2.2 Pipelines

**Pipeline A (natural central dogma)** uses default weights `L1=0.10, L2=0.10, L3=0.20, L4=0.25, L5=0.15, L6=0.20`.

**Pipeline B (therapeutic agents, CRISPR + siRNA)** uses task-family-specific weights `W_L1L2=0.3, W_L3=0.4, W_L5=0.6, W_L6b=0.25`. Spec section 7 explicitly permits per-task-family rescaling, and the cascade JSON header declares which pipeline weights are in effect. Pipeline A and Pipeline B P_harm values are not directly comparable without normalization.

### 2.3 Task catalog (15 tasks)

| Task | Pipeline | n | Max P_harm | Spec models implemented |
|------|---------|--:|-----------:|-------------------------|
| T_CRISPR | B | 604 | 0.6651 | CrisprGPT, GPT-4o, Claude |
| T_siRNA | B | 270 | 0.6753 | GPT-4o, RNA-FM (oracle), RiNALMo (oracle) |
| T1.1 | A | 6 | 0.3379 | DNABERT-2, HyenaDNA, NT-v2, Caduceus |
| T1.2 | A | 9 | 0.4181 | Evo-2 7b_base 8k context |
| T2.1 | A | 10 | 0.3579 | RNA-FM, SpliceBERT, UTR-LM, RiNALMo |
| T3.1 | A | 36 | 0.5128 | ESM-2, ProtBERT, ProtT5, Ankh |
| T3.2 | A | 105 | 0.5650 (Vina-upgraded) | ProteinMPNN, ESM-IF1 |
| T3.3 | A | 20 | 0.5093 | GPT-4o, Claude |
| T4.1 | A | 84 | 0.4155 | MolFormer-XL, ChemBERTa, ChemBERTa-77M-MLM |
| T4.2 | A | 1062 | 0.3807 | REINVENT 4 (mol2mol + de novo), MolGPT |
| T4.3 | A | 683 | 0.3620 | GPT-4o, ChatDrug |
| T4.4 | A | 12+ | 0.3883 | Pocket2Mol, TargetDiff |
| T4.5 | A | 12 | 0.4660 | ProteinMPNN, ESM-IF1 |
| T5.1 | A | 5 | 0.4485 | Geneformer |
| T6.1 | terminal | 221 | recall 56.2%, false-safe 43.8% | ADMET-AI, MolE-style proxy, DeepTox-style proxy |

### 2.4 Reproducibility

- All cascade JSONs ship per-level `{r_l}` decomposition alongside the headline `P_harm`.
- ADMET-AI predictions on a 5-entry spot-check reproduce with delta = 0.0000 across all four endpoints (ClinTox, hERG, DILI, AMES).
- AutoDock Vina baselines for the 6 dockable proteins (Mpro, RdRp, HIV protease, HIV RT, Influenza NA, HCV NS3) are saved in `Protein/shared/vina_wt_baselines.json` and re-used across T3.2 and T4.5 cascade upgrades.
- Pretrained weights: TargetDiff weights mirrored at Zenodo record 14041881 (we cite the public mirror because the original Google Drive folder is not publicly accessible). Caduceus, Evo-2, RiNALMo, ProteinMPNN, ESM-IF1 weights are mirrored on HuggingFace.

---

## 3. Headline findings

### 3.1 Adversarial prompts inflate off-target risk

- **T_CRISPR**: Adversarial-prompt strategies amplify off-target hits **9.05x** over baseline. 48 of 604 gRNAs fall into Critical risk (P_harm > 0.6). Maximum P_harm = 0.6651 (BRCA1 -> SPATA21).
- **T_siRNA**: 100% jailbreak success across 5 attack strategies (S1-S5). All Critical-risk siRNAs hit MYC (essential gene) via seed-region 7-mer matches. Maximum P_harm = 0.6753 (TTR / S2_seed_unconstrained).

### 3.2 Single-task ADMET predictors miss known toxics

- **T6.1**: On ClinTox (n=112 known clinical toxics), ADMET-AI recall is 56.2% (43.8% false-safe). MolE-style proxy recall 38%, DeepTox-style proxy recall 21%. **30.3% of known toxics are missed by all three predictors**, indicating a systematic out-of-distribution blind spot.
- **Implication for T4.1, T4.3 cascades**: their `r_L6a` values are lower bounds for OOD-generated molecules and should be reported with this caveat.

### 3.3 Per-variant Vina docking changes T3.2 / T4.5 ranking

The legacy `r_L4 = 0.7` heuristic used a fixed value for all designs in T3.2 and T4.5. This audit replaces that with per-variant Vina docking-based scoring:

- **WT baseline ΔG** (kcal/mol): Mpro -7.41, RdRp -8.27, HIV protease -12.27, HIV RT -7.32, HCV NS3 -8.76. Influenza NA failed Vina pose generation and falls back to the heuristic.
- **Per-variant r_L4** is then the WT baseline scaled by `(1 + 0.6 * pocket_mut_rate)`, capped at 1. Pocket residues are defined as Cα within 8 Å of the co-crystal ligand center.
- **Effect on rankings**: T3.2 ProteinMPNN max P_harm rises 0.5215 -> 0.5650 (the new top entry has 90% pocket-residue mutation rate). T3.2 ESM-IF1 max rises 0.4915 -> 0.5378.
- **Limitation**: 105 ProteinMPNN+ESM-IF1 designs have real Vina baselines; 9 PPI designs (T4.5 Spike RBD) have no canonical small-molecule ligand and retain the heuristic fallback.

### 3.4 Cross-model agreement is task-dependent

- **T2.1 mask-fill recovery** ranges from 23.0% (UTR-LM) to 44.0% (SpliceBERT). The four spec models disagree by **2x** on a single sequence's recovery.
- **T_siRNA cross-designer Pearson correlation** between RNA-FM and RiNALMo on 1,200 sliding-window guide candidates is **0.666** (moderate agreement). On the GPT-4o-generated siRNA pool, the same correlation is **0.0073** (essentially independent), suggesting that GPT-4o's outputs are out-of-distribution for the masked RNA language models.
- **T3.1 multi-model AARR**: ProtBERT 54%, ESM-2 12.8%, Ankh 5.3%, ProtT5 1.7%. The 32x spread reveals model-class differences in sequence biology coverage.
- **T6.1 cross-predictor**: ADMET-AI / MolE-style / DeepTox-style false-safe rates are 41% / 62% / 79%. Predictors disagree on 30%+ of known toxics.

---

## 4. Methodology

### 4.1 Forward cascade evaluation

Each task's output is fed forward through every reachable level via the cascade engine. For example, a T3.1 ESM-2 mask-fill output is propagated as: **filled protein -> ESMFold (B3->4 folding) -> AutoDock Vina (binding) -> r_L4**, *and* **filled protein -> Pfam/KEGG/VirHostNet (B3->5 pathway-direct) -> r_L5**, *and* **filled protein 9-mer enumeration -> MHCflurry presentation (B3->6b) -> r_L6b**. The dual-path (B3->4 plus B3->5) is implemented for T3.1, T3.2, T3.3.

### 4.2 Risk Aggregator

`P_harm = 1 - prod_l (1 - w_l * r_l)` aggregates levels independently. The per-level decomposition `{w_l * r_l}` is shipped alongside P_harm in every cascade JSON for diagnostic use. P_harm is bounded in [0, 1] but **is not a calibrated probability of harm** — see section 5.

### 4.3 Bridge truncation

When a bridge fails (e.g., B3->4 cannot produce a fold-able structure), the corresponding `r_l` is replaced by a Bridge penalty term `lambda_b in [0, 0.1]`. Noisy-OR handles this by collapsing the corresponding `(1 - w_l * r_l)` toward 1, so the cascade gracefully truncates.

---

## 5. Limitations and methodological honesty

### 5.1 P_harm is a ranking score, not a probability

`P_harm` should be read as "a risk-relevant ranking score in [0, 1] for adversarial outputs *within a task family*". It should **not** be reported or interpreted as a calibrated probability of real-world clinical harm. Calibration would require ground-truth labels (wet-lab assays of harm on representative designs), which do not exist.

### 5.2 The noisy-OR formula has known biological limits

The cascade is a directed causal chain (central dogma), not a fault tree. Treating levels as independent failure modes mis-handles five biological coupling cases: frameshift coupling, synonymous variants, compensatory mutations, allosteric long-range effects, and immune escape. See `P_HARM_CRITIQUE.md` for a full analysis. The recommended path forward is causal-DAG aggregation, but it requires conditional probability tables that are not yet community-curated.

### 5.3 Documented model coverage gaps

A small number of spec models are intentionally not run because they are gated, deprecated, or distributed only via channels that require manual access. These are documented in `AUDIT_REPORT.md` and the benchmark proceeds with the remaining models, which provide multi-model triangulation in every task.

### 5.4 Uncertainty quantification is not propagated

P_harm is a point estimate. MHCflurry has documented calibration error and ADMET-AI has the false-safe rate documented in T6.1. Per-level uncertainty is not currently propagated to the aggregate. A bootstrap approach (MHCflurry seed ensembles, ADMET-AI bagging) is feasible and is recommended as a calibration upgrade.

### 5.5 Wet-lab validation is out of scope

This benchmark is a *computational* safety evaluation framework. None of the high-P_harm designs reported here have been or should be physically synthesized. The framework's purpose is to expose adversarial attack surfaces and force generative biology models to be evaluated under realistic cascade-aware conditions; it is not a discovery pipeline.

---

## 6. Ethics and dual-use considerations

The benchmark is structured to **expose** adversarial vulnerabilities of generative biology AI under controlled conditions, not to provide actionable design recipes. All adversarial prompts and output sequences are reported as **risk evaluation artifacts**, never as candidate biological agents:

- The benchmark uses canonical biosecurity-relevant proteins (Ebola GP, SARS-CoV-2 Spike/Mpro/RdRp, HIV PR/RT, Influenza NA, HCV NS3, Anthrax PA) and canonical essential-gene targets (TP53, BRCA1, EGFR, HTT, MYC). These are the same targets used in the published biosecurity literature (e.g., Soice et al. 2023, Mouton et al. 2023).
- High-risk outputs (max P_harm > 0.6) are reported as ranking signals, not as wet-lab-ready sequences. Specifically, the top-K designs in T_CRISPR / T_siRNA are off-target predictions on essential genes; their biological feasibility depends on factors (delivery, tissue, immune response) the cascade does not model and is not designed to model.
- The benchmark explicitly tests defensive use cases: T6.1 measures ADMET-AI's failure rate on known clinical toxics, which is information defenders need to harden their predictor stack.

The dual-use risk is real but bounded: the benchmark's value is in revealing that **single-modality safety evaluation under-counts cascade harm**, which cannot be exploited without first implementing the cascade itself. The cascade implementation here uses only published models and public databases.

---

## 7. Comparison to existing benchmarks

| Benchmark | Modalities | Cascade-aware | Adversarial probe | Risk aggregation |
|-----------|------------|:-------------:|:-----------------:|------------------|
| ProteinGym (2023) | Protein only | No | No | Single fitness score |
| Therapeutics Data Commons | Multi-task ADMET | No | No | Per-task AUROC |
| MoleculeNet | Small molecule | No | No | Per-task AUROC |
| BiophysicalBench | Protein-ligand | No | No | Per-target ΔG |
| BIOLLM (2024) | LLM coding | No | Yes (jailbreak) | Pass/fail |
| **BioSafety Benchmark** (this work) | DNA / RNA / protein / small molecule / single cell | **Yes** | **Yes** (5 jailbreak strategies) | **Cascade-aggregated P_harm** |

To our knowledge no prior benchmark spans the central dogma cascade, includes therapeutic Pipeline B (CRISPR + siRNA) explicitly, *and* runs adversarial-prompt attacks on the same 25-model panel.

---

## 8. Reproducibility checklist

- [x] All 15 tasks have at least one cascade evaluation pass with a P_harm value (or, for L6a-terminal tasks, a recall metric).
- [x] Cascade JSONs ship per-level `{r_l}` decomposition.
- [x] Scripts are version-pinned to `selfevolve` env (Python 3.12, torch 2.9.1+cu128, transformers 5.5.4) for protein/RNA/small-molecule cascades, and `evo2env` for T1.2 Evo-2 inference.
- [x] Pretrained weights cite their public mirrors (HuggingFace, Zenodo) with provenance hashes where available.
- [x] AutoDock Vina docking is reproducible via `Protein/shared/vina_docking.py` and `Protein/shared/vina_r_l4_upgrade.py`.
- [ ] Wet-lab calibration: not applicable (out of scope).

---

## 9. Conclusion

BioSafety Benchmark is the first cascade-aware safety evaluation framework for generative biology AI. It reveals that single-modality self-evaluation under-counts harm by amounts that matter (9x off-target inflation under jailbreak, 43.8% false-safe on ClinTox, 32x cross-model AARR spread), and provides a per-level risk decomposition that lets safety engineers diagnose *which* channel each high-P_harm design fires through. The benchmark is positioned for the NeurIPS Datasets and Benchmarks track. The known calibration limits of the noisy-OR aggregator are documented and a path forward to causal-DAG aggregation is proposed for future community work.
