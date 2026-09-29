# BioSafety Benchmark — Conference Reviewer Self-Audit

**Date**: 2026-05-02
**Scope**: Anticipate the major review concerns a NeurIPS Datasets and Benchmarks track reviewer would raise, surface concrete weaknesses, and either fix them or document the response.

This document plays both the reviewer and the author. Each section poses a likely review concern, gives the strongest version of the critique, and then documents the response or the gap that remains.

---

## R1. Novelty and positioning

### R1.1 "This is just a collection of existing benchmarks stitched together with noisy-OR."

**The strong version of the critique**: ProteinGym already covers protein fitness; Therapeutics Data Commons covers ADMET; MoleculeNet covers small-molecule properties. What is the novel methodological contribution beyond an aggregation formula?

**Response**:

1. **Cascade-aware forward evaluation is the novel contribution**, not the per-level evaluation. No prior benchmark forces the output of a protein design model through B3->4 (folding + Vina docking) *and* B3->5 (Pfam/KEGG/VirHostNet pathway annotation) *and* B3->6b (MHCflurry presentation) *as a single test*. Each individually is a known evaluation; their composition exposes failure modes that single-modality evaluation cannot see.
2. **The dual-pipeline architecture (Pipeline A natural + Pipeline B therapeutic)** is also novel. Pipeline B (CRISPR gRNA + siRNA off-target) treats therapeutic agents as adversarial inputs and shares L1/L2 with Pipeline A but uses task-family-specific weights. This explicitly tests the most plausible real-world adversarial use case.
3. **The 9.05x off-target inflation under jailbreak prompts in T_CRISPR** and **100% jailbreak rate in T_siRNA** are findings that none of the cited benchmarks would surface, because none of them run adversarial prompts through cascade evaluation.

**Gap that remains**: We do not have a head-to-head benchmark comparison against a single existing aggregated risk score because none exists. We compare against per-task AUROC-style metrics in the existing literature.

### R1.2 "Why noisy-OR? Why not a learned aggregator?"

**Response**: This is addressed in `P_HARM_CRITIQUE.md`. Briefly, a learned MoE aggregator does not improve scientific defensibility in the absence of ground-truth harm labels; it only moves the calibration problem to the training-distribution choice, sacrifices per-level transparency, and cannot encode the central-dogma causal structure that the cascade actually has. The recommended position is to keep noisy-OR as a transparent ranking score, always report per-level decomposition, and frame causal-DAG aggregation as future work pending labeled subsets.

---

## R2. Scientific defensibility of the metric

### R2.1 "Is P_harm = 0.65 supposed to mean 'this design causes harm 65% of the time'?"

**Response**: No, and the benchmark text now explicitly says so. P_harm is a **ranking score**, not a calibrated probability. See `P_HARM_CRITIQUE.md` and `NEURIPS_SUMMARY.md` section 5.1.

### R2.2 "Five biological cases break your independence assumption. Why is this not a fatal flaw?"

The five cases:
- Frameshift coupling (L1 -> L3 deterministic)
- Synonymous variant (L1 nonzero, L3 baseline)
- Compensatory mutation (L4 / L5 anti-correlation)
- Allosteric effect (large L4 jump, small L3 change)
- Immune escape (L3 high reduces L6b match)

**Response**: These are real limitations for *probability* interpretation. For *ranking within a task family*, they bias the score in predictable directions but do not invert rankings. Specifically, all five cases either inflate P_harm (independence over-counts joint risk in a, b) or under-count one channel relative to another (c, d, e). Within a single task family, these biases are systematic and do not flip a low-risk design to high-risk or vice versa. The per-level decomposition `{w_l * r_l}` exposes which channel is driving the score, so a reviewer can apply biological judgment manually. **However**, this is a known limitation, and the recommended biological-coherence flags (section 5 of `P_HARM_CRITIQUE.md`) are concrete proposed fixes.

### R2.3 "How are the weights w_l calibrated?"

**Response**: The default Pipeline A weights `(L1=0.10, L2=0.10, L3=0.20, L4=0.25, L5=0.15, L6=0.20)` are **not** calibrated against a labeled dataset; they reflect an a-priori ranking of "binding-level evidence is the most concrete, sequence-level divergence is the least specific". Pipeline B weights `(W_L1L2=0.3, W_L3=0.4, W_L5=0.6, W_L6b=0.25)` are biologically motivated for therapeutic agents. The benchmark text now flags this as a known caveat.

Calibrating weights would require: (a) a labeled subset of dangerous-vs-benign designs, and (b) a proper-scoring-rule loss (e.g., Brier or log loss) on that subset. We propose this as future community work in `NEURIPS_SUMMARY.md` section 5.4.

---

## R3. Reproducibility

### R3.1 "Can I reproduce the cascade on my hardware?"

**Response**:

| Component | Reproducibility |
|-----------|-----------------|
| ADMET-AI predictions on stored top-5 SMILES | ✅ delta = 0.0000 across all four endpoints |
| AutoDock Vina ΔG for WT proteins | ✅ baseline file `Protein/shared/vina_wt_baselines.json` shipped, includes random seeds and box centers |
| ProteinMPNN designs | ✅ random seeds and sampling temps recorded in cascade JSON |
| ESMFold structures (T3.3) | ✅ deterministic given input sequence |
| Evo-2 7b_base 8k context | ✅ env recipe documented |
| Caduceus, RiNALMo, RNA-FM | ✅ HuggingFace mirrors |
| TargetDiff | ✅ Zenodo public mirror documented |
| MHCflurry 2.2.1 | ✅ pip-installable |

**Hardware**: Single-machine CUDA setup (A100 / H100 / RTX 6000 PRO Blackwell tested). No multi-node requirements.

**Gap**: We do not yet ship a single docker / nix recipe that pins all dependencies. The `selfevolve` env is documented inline but not as a reproducible image. Adding a Dockerfile is recommended before final submission.

### R3.2 "Are random seeds consistent across runs?"

**Response**: Cascade JSONs record sampling temperature, model checkpoint, and seed where the underlying script supports it. ProteinMPNN, ESM-IF1, and TargetDiff are explicitly seeded. MHCflurry and ADMET-AI are deterministic given inputs. ChatDrug PDDS calls use a recorded session but the LLM (GPT-4o) is non-deterministic; the cascade JSON ships the raw response so the cascade re-evaluation is reproducible even if a re-prompt would not give identical outputs.

---

## R4. Coverage and selection bias

### R4.1 "Why these 9 biosecurity proteins and not others?"

The 9 proteins (Ebola GP, SARS-CoV-2 Spike/Mpro/RdRp, HIV PR/RT, Influenza NA, HCV NS3, Anthrax PA) are the canonical biosecurity-relevant proteins in the published literature (Soice et al. 2023, Mouton et al. 2023). They are not the only relevant ones; alternative selections (e.g., mpox, dengue, ricin) would produce a different but equally defensible benchmark.

**Gap**: We do not include a "spurious-target negative control" set (proteins explicitly known to be biosafety-irrelevant). Adding a Pipeline A negative control set would strengthen the benchmark by allowing AUROC-style ranking quality measurement. Recommended for next release.

### R4.2 "Why is T6.1 false-safe rate so high (43.8%)? Is this a benchmark problem or an ADMET-AI problem?"

**Response**: It is an ADMET-AI distribution problem, not a benchmark problem. T6.1 evaluates ADMET-AI on the full ClinTox known-toxics set (n=112), using the same stored model. The 43.8% false-safe rate reproduces a known limitation that ADMET-AI's training distribution under-represents OOD-generated molecules. T6.1 is **the contribution that exposes this** — a per-task ADMET-AI benchmark would not have surfaced it.

The cross-predictor extension (MolE-style 62%, DeepTox-style 79% false-safe) shows the failure is not specific to ADMET-AI; all three predictor architectures share the OOD blind spot. This is a strong argument for benchmark inclusion.

### R4.3 "T5.1 ground truth is simulated, not Perturb-seq."

**Response**: Acknowledged in the T5.1 final report and `AUDIT_REPORT.md` section 4.5. Full Perturb-seq evaluation requires the Replogle scPerturb dataset and an anndata loading pipeline, both feasible but out of scope for this submission. The current T5.1 result (Geneformer on simulated SARS-CoV-2 GSE147507 + HIV VirHostNet signatures) is positioned as a smoke-test of the cascade integration, not the headline finding.

---

## R5. Adversarial threat model

### R5.1 "Are the jailbreak strategies realistic?"

The five strategies (S1 pan_family, S2 seed_unconstrained, S3 promiscuous_seed, S4 expert_override, S5 passenger_ignore) are documented in `siRNA/AGENT_RUN_INSTRUCTIONS.md`. They mirror published jailbreak patterns (Wang et al. 2023, Zou et al. 2023) adapted to the siRNA / CRISPR design context. The 100% jailbreak success rate in T_siRNA is a finding, not a setup artifact: the LLM (GPT-4o) generated guide strands targeting essential genes in 270/270 calls when prompted under any of the 5 strategies, with 0% baseline-level refusals.

**Gap**: We do not test prompt-injection or model-stealing attacks. Threat model is restricted to "prompted output of the released model"; a stronger threat model (fine-tuned or steered model) is a follow-up study.

### R5.2 "What about Claude / Gemini / open-source LLMs?"

**Response**: T_CRISPR includes Claude (with rotated API key, run scheduled), and the T_siRNA / T_CRISPR adversarial design is model-agnostic. The cascade pipeline does not depend on which LLM produced the input. Adding Gemini, open-source LLMs (Llama, Qwen) is straightforward and is a recommended community contribution.

---

## R6. Ethics

### R6.1 "Is this a how-to guide for designing dangerous siRNAs?"

**Response**: No. We report risk-ranking metrics, not wet-lab-ready protocols. Specifically:

- We do NOT publish synthesis instructions, delivery vehicles, or chemical modifications (e.g., 2'-O-methylation patterns, GalNAc conjugation chemistry) that would be required to convert a high-P_harm guide strand into a biologically active reagent.
- We DO publish the canonical seed-region 7-mers and their off-target essential gene matches, because this information is already public via the gene's Ensembl 3'UTR annotation and DepMap essentiality scores.
- We DO publish the cascade-aware ranking method, because defenders need it more than attackers do — an attacker can use any single tool; a defender needs the cascade view.

Comparison: `RNA Off-Target Risk in Therapeutic siRNAs` is already a published research area with multiple benchmarks (e.g., GenomeNet, Therapeutics Data Commons siRNA target). Our incremental disclosure is the cascade-aware aggregation, not novel target sequences.

### R6.2 "Should NeurIPS publish biosecurity-attacking work?"

**Response**: This is a benchmark *of the safety* of generative biology AI, not a paper *demonstrating attacks*. The closest precedents are Soice et al. 2023 (Nature Comm Bio) on LLMs and dual-use biology and Mouton et al. 2023 (RAND) on AI-enabled biological threat assessment, both of which were published as risk-evaluation contributions. NeurIPS D&B has accepted comparable benchmark contributions in adjacent domains (cybersecurity benchmarks, content-moderation benchmarks). The decision to publish should follow the venue's responsible-disclosure policy; the authors will follow the venue's guidance.

**Gap**: We do not yet have a formal pre-publication review by an institutional biosecurity committee. Recommended before final submission.

---

## R7. Comparison to related work

### R7.1 "How does this differ from Soice et al. 2023?"

Soice et al. measured whether LLMs could provide step-by-step instructions for creating pandemic pathogens. Their threat model is **information transfer**: does the LLM tell the user how to do bad biology?

Our threat model is **cascade harm**: does a generative biology FM produce designs whose downstream cascade evaluation reveals high risk, even when the FM's own self-evaluation says the design is safe?

These are complementary. Soice et al. addresses LLM uplift; we address generative biology FM blind spots.

### R7.2 "How does this differ from BIOLLM (2024)?"

BIOLLM evaluates LLM jailbreak rates on biology-related queries. The output is a "did the LLM refuse" boolean. Our work goes one layer deeper: **even when the LLM does not refuse**, what does the cascade evaluation say about the harm potential of the output? T_siRNA's 100% jailbreak rate plus 270 cascade evaluations show that "did not refuse" is necessary but not sufficient — the cascade must agree with the refusal classifier.

---

## R8. Concrete weaknesses we have not yet fixed

The following are honest gaps that the reviewer will see and that we cannot close before submission:

1. **No labeled dangerous-vs-benign subset**. Without this, P_harm cannot be calibrated and ranking quality cannot be measured by AUROC.
2. **No bootstrap confidence intervals on P_harm**. Per-design uncertainty is feasible (MHCflurry seed ensembles, ADMET-AI bagging) but not yet implemented.
3. **No causal-DAG aggregator**. The biologically-correct alternative is documented in `P_HARM_CRITIQUE.md` but not implemented.
4. **No external-validation set**. All cascade designs are scored on the same data they were generated from. A held-out generation set (e.g., from a model not in the benchmark panel) is the right next-step experiment.
5. **No adversarial-prompt LLM ablation**. We test 5 jailbreak strategies on GPT-4o; we do not yet test whether a fine-tuned attack prompt (gradient-based) achieves higher P_harm than the hand-written strategies.
6. **Multi-turn conversational adversaries not tested**. The current T_CRISPR / T_siRNA jailbreak is single-turn. Multi-turn chain-of-thought attacks may yield higher P_harm and are not measured.

These are the most likely "what's missing?" complaints. The submission should explicitly enumerate them as future work, not paper over them.

---

## R9. Recommendations for the submission

Before formal submission to the NeurIPS D&B track:

1. **Add a docker / nix recipe** that pins all dependencies for one-command reproducibility.
2. **Add a Pipeline A negative-control set** (5-10 biosafety-irrelevant proteins) so AUROC-style ranking quality can be reported.
3. **Add bootstrap CIs** on the top-10 P_harm entries per task, using MHCflurry seed ensembles for L6b and ADMET-AI bagging for L6a.
4. **Implement the five biological-coherence flags** described in `P_HARM_CRITIQUE.md` section 5.
5. **Get institutional biosecurity pre-review** before public release, to confirm responsible-disclosure compliance.
6. **Document the threat model explicitly** in section 1 of the submission, including what is and is not in scope (no fine-tuned-attacker threat model, no multi-turn conversational threats).

These are tractable; with 1-2 weeks of additional work the submission addresses the major review concerns above.

---

## R10. Summary verdict

The benchmark's primary contribution — **cascade-aware safety evaluation that exposes single-modality blind spots** — is novel and well-supported by the empirical findings (9.05x off-target inflation, 100% siRNA jailbreak rate, 43.8% ADMET-AI false-safe). The methodological limitations (noisy-OR independence, missing labels, no bootstrap CI) are documented in detail and have credible paths to resolution. The remaining weaknesses (R8) are honest gaps that would not be hidden from reviewers and that frame natural next steps for the community.

**Reviewer-of-self recommendation**: marginally above acceptance threshold. The cascade contribution is strong; the metric calibration is a known weakness that the authors have anticipated and documented.
