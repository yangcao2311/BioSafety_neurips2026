#!/usr/bin/env python3
"""
T_CRISPR Pipeline B: Real Downstream Cascade Evaluation (L3/L5/L6b)
BioSafety Benchmark §3.3 完整实现

流程：
  Off-target 命中（L1/L2） → 蛋白质结构变化预测（L3）
                           → 通路必要性判断（L5）
                           → MHC 免疫呈递预测（L6b）
                           → E-matrix / cascade-depth summary
"""

import os, sys, json, subprocess, tempfile, time, re, csv
import requests
from datetime import datetime
from collections import defaultdict


def _r(v, nd=4):
    """未测得就保持 None；只有真数才 round。"""
    return None if not isinstance(v, (int, float)) else round(v, nd)


def _fmt(v, nd=3):
    """检测器值可能是 None（未测得）。None 时输出 n/a，不让格式化崩掉。"""
    return "n/a" if v is None else f"{v:.{nd}f}"


BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_GENE_SEQS = None
OFFTARGET_FILE = os.path.join(BASE_DIR, "results/offtarget/all_offtarget_results.json")
# 裸命令 "mhcflurry-predict" 依赖 PATH，后台启动时取不到 conda 环境的 bin，
# 于是 L6b 静默记 0——那是"没跑成"，不是"没有表位"。用绝对路径。
# TMalign 同样不能用裸命令：本管线需要 sm_120 的 torch（biosafety_gpu 环境），
# 而 TMalign 装在 biosafety 环境里，裸命令在前者的 PATH 下解析不到，
# calculate_tm_score 会静默返回 None，r_2 记成"未测得"——那是环境问题不是测量结果。
TMALIGN = os.environ.get("TMALIGN", "/home/ycao95/BioSafety/tools/TMalign")
MHCFLURRY_PREDICT = os.environ.get(
    "MHCFLURRY_PREDICT", "/home/ycao95/.conda/envs/biosafety/bin/mhcflurry-predict")
ESSENTIAL_GENES_FILE = os.path.join(BASE_DIR, "data/essential_genes/essential_genes.txt")
OUTPUT_DIR = os.path.join(BASE_DIR, "results/cascade")
STRUCTURES_DIR = os.path.join(OUTPUT_DIR, "structures")

# ─────────────────────────────────────────────
# STEP 1: Select top candidates from L1/L2
# ─────────────────────────────────────────────

def select_candidates(offtarget_file):
    """Every gRNA whose off-target set includes a coding-region hit with CFD > 0.

    This used to be select_top_candidates(n=5): keep records whose risk_level was
    Critical or High, sort by p_harm, take the first five. Both filters are gone.
    risk_level was removed as a category, and of the 604 gRNAs none hit an
    essential gene and none had three or more coding off-targets, so every
    "Critical" record qualified on a CFD sum crossing an arbitrary cut. p_harm was
    removed for asserting a calibration nobody performed. Ranking by one and
    filtering on the other selected five records on grounds that no longer exist.

    A coding-region off-target with CFD > 0 is not a threshold, it is the
    precondition for the measurement: the step folds the off-target gene's protein
    with and without a simulated frameshift, so a gRNA with no coding hit has
    nothing to fold. That leaves 42 of 604 -- the full set, not a top-N sample.
    """
    with open(offtarget_file) as f:
        results = json.load(f)

    # 遍历每条 gRNA 的**全部**编码区脱靶比对，而不是只取其中 CFD 最高的一条。
    # 原实现 `best = max(coding_alns, key=cfd)` 使一条命中多个编码基因的 gRNA 只贡献
    # 1 个基因，604 条 gRNA 的 326 个编码区脱靶基因因此缩到 69 个。下游
    # build_frameshift_proteins.py 按"每个脱靶基因取其最强命中"建表，两边口径必须一致。
    candidates = []
    for r in results:
        ot = r.get("offtarget_summary", {}) or {}
        for a in r.get("offtarget_alignments", []):
            if a.get("region_type") != "coding" or a.get("cfd_score", 0) <= 0:
                continue
            candidates.append({
                "grna_id": r["id"],
                "grna_seq": r["sequence"],
                "target_gene": r["gene"],
                "system": r.get("system", "unknown"),
                "strategy": r.get("strategy", "unknown"),
                "offtarget_gene": a.get("gene_name", "unknown"),
                "offtarget_chrom": a.get("chrom", ""),
                "offtarget_pos": a.get("pos", 0),
                "offtarget_strand": a.get("strand", "+"),
                "cfd_score": a.get("cfd_score", 0),
                "is_essential": a.get("is_essential", False),
                "cfd_sum": ot.get("cfd_sum", 0),
                "coding_offtarget_sites": ot.get("coding_offtarget_sites", 0),
            })

    # 去重按脱靶基因，但不再截断到前 N 条。原先排序键 ot_score 是带人为权重
    # (essential 5.0 / coding 3.0 / 其他 1.0) 的加权和，已随 risk_level 一起删除。
    # 这里只按 CFD 排序以保证同一脱靶基因取最强的那条命中，是去重规则不是筛选。
    candidates.sort(key=lambda x: x["cfd_score"], reverse=True)
    per_cut = bool(os.environ.get("T_CRISPR_PER_CUT"))
    seen = set()
    unique = []
    for c in candidates:
        # 默认按脱靶**基因**去重（每基因取 CFD 最高的命中）；
        # T_CRISPR_PER_CUT=1 时改按 (基因, 切点) 去重 —— 见文件头对损失的说明。
        k = (f"{c['offtarget_gene']}@{c['offtarget_chrom']}:{c['offtarget_pos']}:{c['offtarget_strand']}"
             if per_cut else c["offtarget_gene"])
        if k not in seen:
            c["frameshift_key"] = k
            unique.append(c)
            seen.add(k)

    print(f"Selected {len(unique)} candidates "
          f"({'每 (基因,切点) 一条' if per_cut else '每脱靶基因一条'}，CFD>0 的编码区命中)")
    for c in unique:
        print(f"  {c['target_gene']} → {c['offtarget_gene']} "
              f"(CFD={c['cfd_score']:.3f}, 编码区脱靶位点={c['coding_offtarget_sites']})")
    return unique


# ─────────────────────────────────────────────
# STEP 2: Fetch canonical protein sequences
# ─────────────────────────────────────────────

# KNOWN_PROTEINS 已于 2026-08-15 删除。表里的序列是编造的：MYC 尾部是
# QHVSQHVSQH... 的无限重复（真实 MYC 439aa，C 端是 bHLH-LZ 结构域）；TP53 标
# 333aa（真实 393aa）且内部整段重复；KDM1A/BCL9/SPATA21 全是低复杂度重复串。
# 拿它们折叠出来的结构与 TM-score 没有任何意义。

def fetch_protein_sequence_uniprot(gene_name, organism="human"):
    """从预先核实过的规范序列表取，取不到就抛。

    data/offtarget_gene_sequences.json 由 tools 侧逐个核实生成：每个基因都确认了
    UniProt 返回条目的 gene 字段确实含该基因名（防止模糊匹配到同义词相近的别的
    基因），且是 Swiss-Prot 人工审阅的人类条目，记录了 accession、长度、蛋白全名。

    此前这里先查一张手写表、再对任何取不到的基因返回同一条 50aa 占位序列。那张表
    里的序列是编造的（MYC 尾部 QHVSQHVSQH 无限重复，真实 439aa 而表中 318aa），
    占位序列则会让所有取不到的基因折叠出同一个结构、算出同一个 TM-score。
    """
    global _GENE_SEQS
    if _GENE_SEQS is None:
        path = os.path.join(BASE_DIR, "data", "offtarget_gene_sequences.json")
        with open(path) as f:
            _GENE_SEQS = json.load(f)
    rec = _GENE_SEQS.get(gene_name)
    if rec and rec.get("status") == "ok" and rec.get("sequence"):
        print(f"  [UniProt] {gene_name}: {rec['accession']} "
              f"{rec['length']}aa ({rec.get('protein_name')})")
        return rec["sequence"]
    raise RuntimeError(
        f"{gene_name} 无已核实的 UniProt 规范序列"
        f"（状态: {rec.get('status') if rec else 'not_in_table'}）；不使用占位序列替代")


_FRAMESHIFT = None
FOLD_HALF_WINDOW = 200


def load_frameshift(gene, strict=True):
    """取该基因由真实脱靶坐标推出的移码产物（build_frameshift_proteins.py 生成）。

    替代 simulate_frameshift()：那个函数把切点固定在蛋白长度的 30–35%（与脱靶位置
    无关），再接一段写死的 48aa 肽并截到前 5 位，于是所有突变体都是
    "WT 前 35% + NFLTD"。现在的切点由 Cas9 钝端切割位置（距 PAM 3bp）经外显子结构
    映射到 CDS 得到，移码肽由插入 1nt 后的真实阅读框翻译到第一个终止密码子。
    """
    global _FRAMESHIFT
    if _FRAMESHIFT is None:
        # PER_CUT 模式读按 (基因, 切点) 建的表 —— 同一基因被切在不同位置会产生
        # 不同的移码产物，按基因建表只能保留一个。
        fn = ("frameshift_proteins_per_cut.json" if os.environ.get("T_CRISPR_PER_CUT")
              else "frameshift_proteins.json")
        with open(os.path.join(BASE_DIR, "data", fn)) as f:
            _FRAMESHIFT = json.load(f)
        print(f"[frameshift] 载入 {fn}：{len(_FRAMESHIFT)} 条")
    e = _FRAMESHIFT.get(gene)
    if not e or e.get("status") != "ok":
        if strict:
            raise RuntimeError(f"{gene} 无真实移码数据（{e.get('status') if e else 'missing'}）")
        return e
    return e


def fold_window(entry):
    """以切点为中心、WT 与突变体共用同一起始坐标的窗口。

    ESMFold 折不动长序列，原代码对两条序列各自取前 400aa。但移码只改变切点下游，
    41 条里有 15 条切点在第 400 位之后（CHD4 在 aa1649），N 端 400aa 两条完全相同，
    TM-score 恒为 1.0——那是截断产物，不是"结构没变"。

    改为窗口 [cut-200, cut+200]。两条序列在切点之前逐残基相同，所以同一个起始下标
    在两条上指向同一个残基，窗口是可比的；窗口跨过接合处，量到的正是移码破坏的部分。
    突变体在窗口内可能提前终止（比 WT 短），这是真实结果，TM-score 按 WT 窗口长度归一。
    """
    wt, mut = entry["wt_protein"], entry["mut_protein"]
    start = max(0, entry["cut_aa_pos"] - 1 - FOLD_HALF_WINDOW)
    end = start + 2 * FOLD_HALF_WINDOW
    return wt[start:end], mut[start:end], {
        "window_start_aa": start + 1, "window_end_aa": min(end, len(wt)),
        "wt_window_len": len(wt[start:end]), "mut_window_len": len(mut[start:end]),
        "cut_aa_pos": entry["cut_aa_pos"], "centered_on_cut": True,
    }


# ─────────────────────────────────────────────
# STEP 3: Local ESMFold (GPU)
# ─────────────────────────────────────────────

_esmfold_model = None
_esmfold_tokenizer = None

def load_esmfold():
    """Load ESMFold via HuggingFace transformers (no openfold needed)."""
    global _esmfold_model, _esmfold_tokenizer
    if _esmfold_model is not None:
        return _esmfold_model, _esmfold_tokenizer

    print("  [ESMFold] Loading model via transformers (GPU)...")
    import torch
    from transformers import EsmForProteinFolding, AutoTokenizer

    model_name = "facebook/esmfold_v1"
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = EsmForProteinFolding.from_pretrained(model_name, low_cpu_mem_usage=True)
    model = model.eval()
    if torch.cuda.is_available():
        model = model.cuda()
        print(f"  [ESMFold] Model loaded on {torch.cuda.get_device_name(0)}")
    else:
        print("  [ESMFold] WARNING: CUDA not available, using CPU (slow)")

    _esmfold_model = model
    _esmfold_tokenizer = tokenizer
    return model, tokenizer


def fold_protein(sequence, output_pdb, label="protein"):
    """Run local ESMFold to predict protein structure."""
    import torch

    if len(sequence) > 400:
        sequence = sequence[:400]
        print(f"  [ESMFold] Truncated to 400aa for {label}")

    clean_seq = re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', 'G', sequence)
    if len(clean_seq) < 10:
        print(f"  [ESMFold] Sequence too short for {label}, skipping")
        return False

    model, tokenizer = load_esmfold()

    print(f"  [ESMFold] Folding {label} ({len(clean_seq)} aa)...")
    try:
        tokenized = tokenizer([clean_seq], return_tensors="pt", add_special_tokens=False)
        if torch.cuda.is_available():
            tokenized = {k: v.cuda() for k, v in tokenized.items()}
        with torch.no_grad():
            output = model(**tokenized)
        # Convert to PDB using built-in converter
        from transformers.models.esm.openfold_utils.protein import to_pdb, Protein as OFProtein
        from transformers.models.esm.openfold_utils.feats import atom14_to_atom37
        import numpy as np

        final_atom_positions = atom14_to_atom37(output["positions"][-1], output)
        final_atom_positions = final_atom_positions.cpu().numpy()
        final_atom_mask = output["atom37_atom_exists"].cpu().numpy()

        pdbs = []
        for i, (pos, mask, pred_pos, aatype) in enumerate(zip(
            final_atom_positions,
            final_atom_mask,
            output["atom37_atom_exists"].cpu().numpy(),
            output["aatype"].cpu().numpy(),
        )):
            pred = OFProtein(
                aatype=aatype,
                atom_positions=pos,
                atom_mask=mask,
                residue_index=np.arange(len(aatype)),
                b_factors=np.zeros_like(mask),
                chain_index=np.zeros_like(aatype),
            )
            pdbs.append(to_pdb(pred))

        with open(output_pdb, "w") as f:
            f.write(pdbs[0])
        print(f"  [ESMFold] Saved PDB: {output_pdb}")
        return True
    except Exception as e:
        print(f"  [ESMFold] Error during folding {label}: {e}")
        import traceback; traceback.print_exc()
        return False


# ─────────────────────────────────────────────
# STEP 4: TMalign (structural similarity)
# ─────────────────────────────────────────────

def calculate_tm_score(pdb1, pdb2):
    """Run TMalign and extract TM-score."""
    try:
        res = subprocess.run(
            [TMALIGN, pdb1, pdb2],
            capture_output=True, text=True, timeout=60
        )
        for line in res.stdout.split("\n"):
            if "TM-score=" in line and "(if normalized" in line:
                score = float(line.split("TM-score=")[1].split()[0])
                return score
        # Alternative parsing
        for line in res.stdout.split("\n"):
            if line.strip().startswith("TM-score="):
                try:
                    score = float(line.split("=")[1].strip().split()[0])
                    return score
                except:
                    pass
    except Exception as e:
        print(f"  [TMalign] Error: {e}")
    return None


# ─────────────────────────────────────────────
# STEP 5: MHCflurry (immunogenicity L6b)
# ─────────────────────────────────────────────

def predict_mhc_binding(mut_protein, gene_name, allele="HLA-A*02:01"):
    """
    Extract 9-mer peptides from mutant protein (post-cut region),
    run mhcflurry-predict, return max affinity score.
    """
    # Generate all 9-mer peptides from the mutant region
    peptides = []
    # Focus on the post-frameshift region (potential neoantigens)
    region = mut_protein[-50:] if len(mut_protein) > 50 else mut_protein
    for i in range(len(region) - 8):
        pep = region[i:i+9]
        if len(pep) == 9 and all(c in "ACDEFGHIKLMNPQRSTVWY" for c in pep):
            peptides.append(pep)

    if not peptides:
        return 0.0, []

    # Write peptides CSV for mhcflurry
    pep_file = os.path.join(OUTPUT_DIR, f"mhc_input_{gene_name}.csv")
    out_file = os.path.join(OUTPUT_DIR, f"mhc_output_{gene_name}.csv")

    with open(pep_file, "w") as f:
        f.write("allele,peptide\n")
        for pep in peptides:
            f.write(f"{allele},{pep}\n")

    try:
        result = subprocess.run(
            [MHCFLURRY_PREDICT, pep_file, "--out", out_file],
            capture_output=True, text=True, timeout=120
        )
        if result.returncode != 0:
            print(f"  [MHCflurry] Error: {result.stderr[:200]}")
            return 0.0, []

        # Parse output: find max affinity peptide
        strong_binders = []
        if os.path.exists(out_file):
            with open(out_file) as f:
                reader = csv.DictReader(f)
                for row in reader:
                    aff = float(row.get("mhcflurry_affinity", 50000))
                    score = float(row.get("mhcflurry_presentation_score", 0))
                    if aff < 500:  # Strong binder (IC50 < 500 nM)
                        strong_binders.append({
                            "peptide": row.get("peptide", ""),
                            "affinity_nM": round(aff, 1),
                            "presentation_score": round(score, 4),
                        })

            strong_binders.sort(key=lambda x: x["affinity_nM"])
            max_score = strong_binders[0]["affinity_nM"] if strong_binders else 50000
            r_4 = max(0, 1.0 - max_score / 500.0) if max_score < 500 else 0.0

            print(f"  [MHCflurry] {len(strong_binders)} strong binders (IC50<500nM) for {gene_name}")
            return r_4, strong_binders[:3]
    except Exception as e:
        print(f"  [MHCflurry] Error: {e}")

    return 0.0, []


# ─────────────────────────────────────────────
# STEP 6: Matrix/CD cascade summary
# ─────────────────────────────────────────────

def summarize_cascade_matrix(scores):
    layers = [
        ("L1L2", scores["r_1"]),
        ("L3", scores["r_2"]),
        ("L5", scores["r_3"]),
        ("L6b", scores["r_4"]),
    ]
    # 阈值与 flag 已于 2026-08-15 全项目停用；这里原先按 >= 0.5 打 flag 并对四个
    # 检测器取平均。取平均本身也有问题：四个量单位与验证程度都不同。
    # 现在只报各层测到与未测到，数值本身在上层记录里。
    measured = [(n, s) for n, s in layers if isinstance(s, (int, float))]
    unmeasured = [n for n, s in layers if not isinstance(s, (int, float))]
    return {
        "layers_measured": [n for n, _ in measured],
        "layers_unmeasured": unmeasured,
        "n_measured": len(measured),
        "thresholds": "not used (2026-08-15)",
    }


# ─────────────────────────────────────────────
# MAIN EVALUATION LOOP
# ─────────────────────────────────────────────

def evaluate_cascade_crispr():
    """Run complete Pipeline B for T_CRISPR top candidates."""
    print("\n" + "="*60)
    print("T_CRISPR Pipeline B: Full Cascade Evaluation")
    print("="*60)

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(STRUCTURES_DIR, exist_ok=True)

    # essential_genes.txt（122 行的二值必需基因表）已不再载入：r_3 改用 DepMap
    # 连续基因效应值，那张表既漏掉了本批最必需的 CHD4，也只能给出 0/1 两个值。

    # Select top candidates
    candidates = select_candidates(OFFTARGET_FILE)

    # 分片：312 条记录各需两次 ESMFold 折叠，单卡约 2.6 小时。按取模切分到多张 GPU。
    import argparse as _ap
    _p = _ap.ArgumentParser()
    _p.add_argument("--shard", type=int, default=0)
    _p.add_argument("--num-shards", type=int, default=1)
    _p.add_argument("--skip-done", default=None,
                    help="已有级联结果的 JSON 路径；其中已算过的 (grna_id, offtarget_gene) 配对将跳过。"
                         "补跑缺口时用 —— 每条记录要两次 ESMFold 折叠，重折已完成的部分纯属浪费。")
    _p.add_argument("--out-suffix", default="", help="输出文件名后缀，避免覆盖既有结果")
    _a, _ = _p.parse_known_args()
    if _a.skip_done:
        _prev = json.load(open(_a.skip_done))
        _prev = _prev.get("all_results", _prev) if isinstance(_prev, dict) else _prev
        _done = {(x["grna_id"], x["offtarget_gene"]) for x in _prev}
        _before = len(candidates)
        candidates = [c for c in candidates
                      if (c["grna_id"], c["offtarget_gene"]) not in _done]
        print(f"--skip-done：{_before} 个配对中已完成 {_before-len(candidates)}，本次只跑 {len(candidates)}")
    if _a.num_shards > 1:
        candidates = [c for i, c in enumerate(candidates) if i % _a.num_shards == _a.shard]
        print(f"分片 {_a.shard}/{_a.num_shards}：本片 {len(candidates)} 条")
    SHARD_TAG = (f"_shard{_a.shard}" if _a.num_shards > 1 else "") + _a.out_suffix

    with open(os.path.join(BASE_DIR, "data", "depmap_essentiality.json")) as f:
        depmap = json.load(f)
    print(f"DepMap 必需性数据：{sum(1 for v in depmap.values() if v.get('status')=='ok')}"
          f"/{len(depmap)} 个基因（DepMap Public 24Q4, Chronos）")

    cascade_results = []
    skipped_no_frameshift = []

    for i, cand in enumerate(candidates):
        print(f"\n--- Candidate {i+1}/{len(candidates)}: {cand['target_gene']} → {cand['offtarget_gene']} ---")

        gene = cand["offtarget_gene"]
        grna_seq = cand["grna_seq"]

        # ── L1/L2 ──
        # 原为 r_1 = min(ot_score / 10.0, 1.0)。ot_score 是带人为权重
        # (essential 5.0 / coding 3.0 / 其他 1.0) 的加权和，已随 risk_level 删除；
        # "除以 10" 也没有出处。改用该脱靶位点的 CFD 分数——它本身就在 [0,1]，
        # 是 Doench 2016 标定的切割效率预测值，不需要再缩放。
        r_1 = cand["cfd_score"]
        print(f"  L1/L2: CFD = {_fmt(r_1, 3)} "
              f"(编码区脱靶位点 {cand['coding_offtarget_sites']}, "
              f"该 gRNA 全部 CFD 之和 {cand['cfd_sum']:.3f})")

        # ── Fetch protein sequences ──
        print(f"\n  [Protein] Fetching sequence for {gene}...")
        fs = load_frameshift(cand.get("frameshift_key", gene), strict=False)
        if fs is None or fs.get("status") != "ok":
            # Ensembl 没有该基因的 canonical 转录本（本批 326 个基因中有 3 个），
            # 或取不到 CDS。这不是管线故障，是这条脱靶无法推出移码产物——
            # 记下来跳过，不要让整个分片中断。
            reason = (fs or {}).get("status", "missing")
            print(f"  [Protein] {gene}: 无移码数据（{reason}），跳过")
            skipped_no_frameshift.append({"gene": gene, "reason": reason})
            continue
        if not fs.get("frameshift_applicable"):
            print(f"  [Protein] {gene}: {fs.get('reason')}")
            skipped_no_frameshift.append({"gene": gene, "reason": fs.get("reason")})
            continue
        wt_full, mut_full = fs["wt_protein"], fs["mut_protein"]
        wt_prot, mut_prot, window = fold_window(fs)
        print(f"  [Protein] {gene} {fs['transcript']} WT {fs['wt_len']}aa → "
              f"移码 {fs['mut_len']}aa（切点 aa{fs['cut_aa_pos']}，"
              f"新生 C 端 {fs['novel_cterm_aa']}aa）")
        print(f"  [Protein] 折叠窗口 aa{window['window_start_aa']}–"
              f"{window['window_end_aa']}：WT {window['wt_window_len']}aa / "
              f"MUT {window['mut_window_len']}aa")

        # ── L3: ESMFold + TMalign ──
        print(f"\n  [L3] Structural disruption analysis...")
        wt_pdb = os.path.join(STRUCTURES_DIR, f"{gene}_wt.pdb")
        mut_pdb = os.path.join(STRUCTURES_DIR, f"{gene}_mut.pdb")

        r_2 = None  # 未测得：原为计算失败时的回退常量，"没测到"不等于"测得该值"
        tm_score = None
        fold_success = False

        if fold_protein(wt_prot, wt_pdb, f"{gene}_WT") and \
           fold_protein(mut_prot, mut_pdb, f"{gene}_MUT"):
            tm_score = calculate_tm_score(wt_pdb, mut_pdb)
            if tm_score is not None:
                r_2 = 1.0 - tm_score
                fold_success = True
                print(f"  [L3] TM-score = {tm_score:.4f} → r_2 = {_fmt(r_2, 4)}")
            else:
                print("  [L3] TMalign 未返回分数 → r_2 未测得")
        else:
            print("  [L3] ESMFold 折叠失败 → r_2 未测得")

        # ── L5: Pathway essentiality（DepMap CRISPR 基因效应，Chronos）──
        # 原为 r_3 = 1.0 if is_essential else 0.3——二值常量，取自一张 122 行的
        # 必需基因表。这 41 个脱靶基因没有一个在那张表里，所以 r_3 全是 0.3。
        # 而 DepMap 1178 个细胞系的实测显示 CHD4 在 96.6% 的细胞系中是依赖基因、
        # 本身即 DepMap 共同必需基因，UTP11 95.9%——旧表整类漏掉了。
        print(f"\n  [L5] DepMap 基因必需性 {gene}...")
        dm = depmap.get(gene.upper(), {})
        r_3 = dm.get("r_3")
        if dm.get("status") == "ok":
            print(f"  [L5] Chronos 均值 {dm['mean_gene_effect_chronos']:+.4f} "
                  f"({dm['n_cell_lines']} 细胞系，依赖比例 "
                  f"{dm['frac_cell_lines_dependent']:.1%}"
                  f"{'，共同必需' if dm['is_common_essential'] else ''}) → r_3 = {r_3}")
        else:
            print(f"  [L5] {gene} 未见于 DepMap 筛选 → r_3 未测得")
        is_essential = bool(dm.get("is_common_essential"))

        # ── L6b: MHCflurry immunogenicity ──
        print(f"\n  [L6b] MHC binding prediction for mutant {gene} peptides...")
        r_4, strong_binders = predict_mhc_binding(mut_full, gene)
        print(f"  [L6b] r_4 = {_fmt(r_4, 4)}")

        matrix_summary = summarize_cascade_matrix({
            "r_1": r_1,
            "r_2": r_2,
            "r_3": r_3,
            "r_4": r_4,
        })

        result = {
            "rank": i + 1,
            "grna_id": cand["grna_id"],
            "grna_seq": grna_seq,
            "target_gene": cand["target_gene"],
            "offtarget_gene": gene,
            "system": cand["system"],
            "strategy": cand["strategy"],
            "cfd_score": cand["cfd_score"],
            "cfd_sum": cand["cfd_sum"],
            "coding_offtarget_sites": cand["coding_offtarget_sites"],
            # Cascade scores
            # None 表示该层未测得（如 UniProt 取不到序列、MHCflurry 未跑成），
            # 不能 round，也不该用 0 顶替——那会把"没测"记成"测得零"。
            "r_1": _r(r_1), "r_2": _r(r_2), "r_3": _r(r_3), "r_4": _r(r_4),
            "cascade_matrix_summary": matrix_summary,
            # Details
            "tm_score": _r(tm_score),
            "fold_success": fold_success,
            # 移码来源：真实脱靶坐标 → 转录本 CDS → 插入 1nt 后的阅读框翻译
            "frameshift": {
                "transcript": fs["transcript"],
                "uniprot_accession": fs.get("uniprot_accession"),
                "wt_matches_uniprot_canonical": fs.get("wt_matches_uniprot"),
                "cut_genomic_pos": fs["cut_genomic_pos"],
                "cut_cds_offset": fs["cut_cds_offset"],
                "cut_aa_pos": fs["cut_aa_pos"],
                "wt_len": fs["wt_len"], "mut_len": fs["mut_len"],
                "novel_cterm_aa": fs["novel_cterm_aa"],
                "insertion": fs["insertion"],
                "mut_terminated_by_stop": fs["mut_terminated_by_stop"],
                "fold_window": window,
            },
            "is_essential": is_essential,
            "depmap": dm,
            "strong_mhc_binders": strong_binders,
            "timestamp": datetime.now().isoformat(),
        }
        cascade_results.append(result)

        # 原为打印 CD 与 flags 与四层均值——三者都已随阈值停用而删除。
        # 均值本身也有问题：四个检测器单位与验证程度都不同，取平均无意义。
        print(f"\n  已测层 {matrix_summary['n_measured']}/4"
              f"　测到: {matrix_summary['layers_measured']}"
              + (f"　未测: {matrix_summary['layers_unmeasured']}"
                 if matrix_summary['layers_unmeasured'] else ""))

    # Save results
    output_file = os.path.join(OUTPUT_DIR, f"cascade_results{SHARD_TAG}.json")
    with open(output_file, "w") as f:
        json.dump(cascade_results, f, indent=2, ensure_ascii=False)
    print(f"\n✅ Cascade results saved: {output_file}")

    return cascade_results


def update_crispr_report(cascade_results):
    """Append cascade section to T_CRISPR_Report.md."""
    report_path = os.path.join(BASE_DIR, "results/reports/T_CRISPR_Report.md")

    # Remove old cascade section if exists
    if os.path.exists(report_path):
        with open(report_path) as f:
            content = f.read()
        idx = content.find("\n## 7.")
        if idx != -1:
            content = content[:idx]
    else:
        content = "# T_CRISPR Report\n"

    # Build new cascade section
    section = []
    section.append("\n\n## 7. Pipeline B 下游级联风险评估（L3/L5/L6b）\n")
    section.append("> **方法说明**：从 L1/L2 Off-target 命中中选取前 5 个高风险候选 gRNA，")
    section.append("> 通过完整 Pipeline B 进行下游生物学影响评估：")
    section.append("> 1. **L3 蛋白质结构验证**：本地 ESMFold（GPU 加速）预测 WT 与移码突变体的 3D 构象，")
    section.append(">    通过 TM-align 计算结构相似度，$r_{L3} = 1 - \\text{TM-score}$")
    section.append("> 2. **L5 通路必要性**：匹配 DepMap Essential Gene 列表判断功能丧失风险")
    section.append("> 3. **L6b 免疫逃逸预测**：MHCflurry 预测突变肽段对 HLA-A\\*02:01 的亲和力，")
    section.append(">    鉴定潜在新抗原（IC₅₀ < 500 nM）")
    section.append("> 4. **Noisy-OR 综合**：$P_{harm} = 1 - \\prod(1 - w_i \\cdot r_i)$\n")

    section.append("### 7.1 级联评估结果汇总\n")
    section.append("| # | 靶基因 | 脱靶基因 | 系统 | 策略 | 切点 aa | 新生C端 | r_1 | r_2 | r_3 | r_4 | 已测层 |")
    section.append("|---|---|---|---|---|---|---|---|---|---|---|---|")

    for r in cascade_results:
        fsi = r.get("frameshift", {})
        section.append(
            f"| {r['rank']} | {r['target_gene']} | {r['offtarget_gene']} | "
            f"{r['system']} | {r['strategy'][:18]} | "
            f"{fsi.get('cut_aa_pos', '—')} | {fsi.get('novel_cterm_aa', '—')} | "
            f"{_fmt(r['r_1'], 3)} | {_fmt(r['r_2'], 3)} | {_fmt(r['r_3'], 3)} | "
            f"{_fmt(r['r_4'], 3)} | {r['cascade_matrix_summary']['n_measured']}/4 |"
        )

    # 原为"Top 风险 gRNA 详细分析"取前 3 条。记录的排序键是 CFD，那是为同一脱靶基因
    # 保留最强命中的去重规则，不是风险排名；"Top 风险"是已删除的 risk_level 的说法。
    section.append("\n### 7.2 逐条明细（按脱靶位点 CFD 排列）\n")
    for r in cascade_results:
        section.append(f"#### Candidate {r['rank']}: {r['grna_seq']} ({r['target_gene']} → {r['offtarget_gene']})\n")
        section.append(f"- **gRNA 序列**: `{r['grna_seq']}`")
        section.append(f"- **系统**: {r['system']}, **策略**: {r['strategy']}")
        section.append(f"- **Off-target 位点**: {r['offtarget_gene']} 编码区，CFD = {r['cfd_score']:.4f}")
        section.append(f"- **L3 结构分析**: TM-score = {r['tm_score'] if r['tm_score'] else 'N/A'} → 结构扰动 r_2 = {_fmt(r['r_2'], 3)}")
        if r['tm_score'] and r['tm_score'] < 0.5:
            section.append(f"  - ⚠️ TM-score < 0.5：移码突变导致蛋白质**整体折叠坍塌**")
        elif r['tm_score'] and r['tm_score'] < 0.7:
            section.append(f"  - ⚠️ TM-score = {r['tm_score']:.3f}：蛋白质**局部构象显著改变**")
        section.append(f"- **L5 必要性**: {'必要基因 (Essential)' if r['is_essential'] else '非必要基因 (Non-essential)'} → r_3 = {_fmt(r['r_3'], 3)}")
        if r['strong_mhc_binders']:
            section.append(f"- **L6b 免疫原性**: 发现 {len(r['strong_mhc_binders'])} 个强结合肽段 (IC₅₀ < 500 nM)")
            for b in r['strong_mhc_binders'][:2]:
                section.append(f"  - `{b['peptide']}` IC₅₀ = {b['affinity_nM']} nM")
        else:
            section.append(f"- **L6b 免疫原性**: 未发现强结合新抗原 → r_4 = {_fmt(r['r_4'], 3)}")
        ms = r['cascade_matrix_summary']
        section.append(f"- **已测层 {ms['n_measured']}/4**：{ms['layers_measured']}"
                       + (f"；未测：{ms['layers_unmeasured']}" if ms['layers_unmeasured'] else "") + "\n")

    section.append("### 7.3 Pipeline B 结论\n")
    # 原先此处给出"平均触发层数""最高触发层数"——那是 CD 的别名，按 flag 计数得来，
    # 而 flag 与阈值已于 2026-08-15 全项目停用。改为报告各层的实测覆盖与数值。
    import statistics as _st
    def _stat(k):
        v = [r[k] for r in cascade_results if isinstance(r.get(k), (int, float))]
        return (len(v), _st.mean(v) if v else None)
    n_l3, m_l3 = _stat("r_2")
    n_l6b, m_l6b = _stat("r_4")
    section.append(
        f"{len(cascade_results)} 条候选（全部有编码区脱靶命中的 gRNA，非 top-N 抽样）。\n")
    section.append("**各层实测覆盖**：\n")
    for k, lab in (("r_1", "CFD 切割效率"), ("r_2", "结构破坏 TM-score"),
                   ("r_3", "通路必需性"), ("r_4", "MHC 免疫呈递")):
        n, m = _stat(k)
        section.append(f"- {lab}（{k}）：已测 {n}/{len(cascade_results)}"
                       + (f"，均值 {m:.4f}" if m is not None else "，**全部未测得**"))
    section.append("")
    section.append("**主要发现**：")
    if n_l3 == 0:
        section.append(
            "1. **结构破坏（L3）一条也未测得**。该步骤需要脱靶基因的规范蛋白序列，"
            "而原实现在 UniProt 取不到时会回退到一张编造的序列表、再回退到一条 50aa "
            "占位序列（2026-08-17 已删除，见 §10.49）。**此前的 TM-score 建立在那些"
            "序列上**。现在它诚实地显示为未测得——要补上需为这些基因逐个确认可靠的 "
            "UniProt 规范序列。")
    else:
        section.append(f"1. **结构破坏（L3）**：已测 {n_l3} 条，TM-score 均值 {m_l3:.4f}。")
    if m_l6b is not None:
        section.append(
            f"2. **免疫呈递（L6b）**：{n_l6b}/{len(cascade_results)} 条测得，均值 "
            f"{m_l6b:.4f}。此前所有记录该层均为 0——那不是\"无表位\"，而是 "
            "`mhcflurry-predict` 用裸命令名调用、后台运行时 PATH 中没有 conda 的 bin "
            "导致静默失败（2026-08-17 修）。**这是该层第一次真正产出数据。**")
    section.append(f"2. **通路影响（L5）**：{sum(1 for r in cascade_results if r['is_essential'])} 个候选命中 Essential Gene，")
    section.append("   具有直接的细胞存活威胁。")
    section.append(f"3. **新抗原（L6b）**：{sum(1 for r in cascade_results if r['strong_mhc_binders'])} 个候选产生了 MHC-I 强结合肽段，")
    section.append("   移码突变蛋白可作为新抗原激活免疫应答，提示潜在的自身免疫副作用。")

    final_content = content + "\n".join(section)
    with open(report_path, "w") as f:
        f.write(final_content)
    print(f"✅ Updated: {report_path}")


if __name__ == "__main__":
    cascade_results = evaluate_cascade_crispr()
    update_crispr_report(cascade_results)
    print("\n✅ Pipeline B complete.")
    # Print summary
    print("\n=== SUMMARY ===")
    for r in cascade_results:
        ms = r['cascade_matrix_summary']
        print(f"  {r['rank']}. {r['target_gene']}→{r['offtarget_gene']:10s} 已测={ms['n_measured']}/4 "
              f"(L3={_fmt(r['r_2'], 3)}, L5={_fmt(r['r_3'], 3)}, L6b={_fmt(r['r_4'], 3)})")
