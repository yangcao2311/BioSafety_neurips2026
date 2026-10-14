#!/usr/bin/env python3
"""T_siRNA：**专用 siRNA 设计模型** OligoFormer 作为设计器，与通用 RNA-LM 并列。

设计口径
--------
唯一的自变量是**用哪个模型排序**。候选池、脱靶通道、检测器全部与
`rna_fm_designer.py` **同函数复用**（直接 import，不重写），因此三组 top-K 的差异
只能来自设计器本身。

为什么要加这一档
----------------
RNA-FM / RiNALMo 给出的是**通用语言模型似然（NLL）**，与"这条 siRNA 能不能沉默
目标基因"没有任何关系 —— 它们只是判别式 oracle，被借来当排序器用。
OligoFormer（Bai et al., Bioinformatics 2024）是在**真实 knockdown 数据**上训练的
siRNA 效力预测器，它的排序反映的是实际沉默效力。

于是本层出现一条与 T_CRISPR 完全对称的轴：
    T_CRISPR  CrisprGPT（专用工具） vs GPT-4o（通用 LLM）
    T_siRNA   OligoFormer（专用设计器） vs RNA-FM / RiNALMo（通用 RNA-LM）

**被检验的假设**：按沉默效力选出的 siRNA，种子区互补性更强，
脱靶负担可能**反而更高** —— 「设计得越好 = 越危险」。
这个假设 RNA-FM 测不出来，因为它根本不看沉默效力。

⚠ 三条必须与结果一起带出的口径
-------------------------------
1. **OligoFormer 也不是生成模型**，同样是判别式 oracle。加它解决的不是
   "被测的不是生成能力"，而是"排序信号与任务无关"。
2. **方向相反**：RNA-FM/RiNALMo 的 NLL **越低越好**，OligoFormer 的 efficacy
   **越高越好**。选择时按各自方向取 top-K，`designer_score` 字段带 `direction`。
3. **长度口径**：本管线候选是 21 nt 引导链（19 nt 双链区 + 2 nt 3′ 突出端），
   OligoFormer 只吃 19 nt。取 `guide[0:19]` 作为双链核心送入，
   种子区（引导链 2–8 位）完整包含在内。此换算写进每条记录的 `of_core_19nt` 字段。

输出
----
  results/oligoformer_designer/of_scores_<gene>.json   每基因全部候选的 efficacy
  results/oligoformer_designer/cascade_results.json    OligoFormer top-K 的级联记录
  results/oligoformer_designer/summary.json            与两个通用 RNA-LM 的三方对比
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

BASE = Path(__file__).resolve().parent.parent
SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))

import rna_fm_designer as RFD  # noqa: E402  —— 候选池 / 通道 / 检测器全部复用

OF_ROOT = Path("/home/ycao95/BioSafety/Task/tools/OligoFormer")
OF_PY = Path("/home/ycao95/.conda/envs/biosafety/bin/python")
OUT_DIR = BASE / "results" / "oligoformer_designer"

DESIGNER_LABEL = "OligoFormer"
DESIGNER_DIRECTION = "higher_is_better"   # efficacy：越高 = 预测沉默效力越强


# ============================================================
# OligoFormer 打分
# ============================================================
def score_with_oligoformer(gene: str, mrna: str, cands: list[dict],
                           cuda: int) -> dict[str, float]:
    """把本基因的候选送 OligoFormer，返回 {candidate_id: efficacy}。

    ⚠ 不使用 OligoFormer 自带的 -off / -tox 脱靶与毒性模块 ——
      本项目的脱靶判据必须与 RNA-FM 那两档**逐字相同**，否则三组不可比。
    """
    cores, keep = [], []
    for c in cands:
        core = c["guide"][:19]          # 21 nt 引导链的 19 nt 双链核心
        if len(core) != 19 or set(core) - set("AUCG"):
            continue                    # 含 N 或长度异常的直接丢弃，不补齐
        c["of_core_19nt"] = core
        cores.append(core)
        keep.append(c)
    if not cores:
        return {}

    work = Path(tempfile.mkdtemp(prefix=f"of_{gene}_"))
    tag = f"OFDES{gene}"                # 名字里不能有下划线以外的怪字符
    (work / "mrna.fa").write_text(f">{tag}\n{mrna}\n")
    (work / "sirna.fa").write_text("\n".join(f">{i}\n{s}" for i, s in enumerate(cores)) + "\n")

    outdir = work / "result"
    outdir.mkdir(exist_ok=True)
    cmd = [str(OF_PY), "scripts/main.py", "--infer", "1",
           "-i1", str(work / "mrna.fa"), "-i2", str(work / "sirna.fa"),
           "--cuda", str(cuda), "--output_dir", str(outdir) + "/"]
    env = {**__import__("os").environ,
           "PATH": f"{OF_PY.parent}:{__import__('os').environ.get('PATH','')}"}
    print(f"  [{gene}] OligoFormer 打分 {len(cores)} 条 ...", flush=True)
    r = subprocess.run(cmd, cwd=OF_ROOT, env=env, capture_output=True, text=True)
    res = outdir / f"{tag}.txt"
    if r.returncode != 0 or not res.exists():
        # 失败就抛，不返回默认分 —— 未测得不等于测得该值
        raise RuntimeError(
            f"{gene}: OligoFormer 推理失败（returncode={r.returncode}）\n"
            f"stdout tail:\n{r.stdout[-1500:]}\nstderr tail:\n{r.stderr[-1500:]}")

    import pandas as pd
    df = pd.read_csv(res, sep="\t")
    if len(df) != len(cores):
        raise RuntimeError(f"{gene}: OligoFormer 返回 {len(df)} 行，送入 {len(cores)} 条，行数对不上")
    # 输出保持送入顺序（infer.py 按 siRNA 列表顺序构表）
    out = {}
    for c, eff in zip(keep, df["efficacy"].tolist()):
        out[c["id"]] = float(eff)
    return out


# ============================================================
# 主流程
# ============================================================
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidate-cap", type=int, default=200,
                    help="每基因候选数上限（必须与 rna_fm_designer 那次运行相同，否则池子不同）")
    ap.add_argument("--top-k", type=int, default=5)
    ap.add_argument("--tissue", default="Liver")
    ap.add_argument("--cuda", type=int, default=0)
    ap.add_argument("--genes", nargs="*", default=None, help="只跑这些基因（冒烟用）")
    ap.add_argument("--skip-channels", action="store_true",
                    help="只打分不跑脱靶通道（冒烟用）")
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    targets = {g: a for g, a in RFD.TARGET_REFSEQS.items()
               if args.genes is None or g in args.genes}
    if not targets:
        sys.exit(f"没有匹配的基因。可选：{list(RFD.TARGET_REFSEQS)}")

    # ---- 1. 候选池：与 rna_fm_designer 同函数、同参数 ----
    all_candidates, mrna_by_gene = [], {}
    for gene, accession in targets.items():
        mrna = RFD.parse_fasta(RFD.fetch_mrna(gene, accession))
        if len(mrna) < 40:
            raise RuntimeError(f"{gene} mRNA 只有 {len(mrna)} nt，不足以滑窗。")
        mrna_by_gene[gene] = mrna
        cands = RFD.generate_guide_candidates(mrna, args.candidate_cap)
        for c in cands:
            c.update({"gene": gene, "accession": accession,
                      "id": f"des_{gene}_{c['target_pos']}"})
        all_candidates.extend(cands)
        print(f"  {gene} ({accession}): mRNA={len(mrna)} nt，候选 {len(cands)} 条")

    # ---- 2. OligoFormer 打分 ----
    print(f"\n=== {DESIGNER_LABEL} 打分（efficacy，越高越好）===")
    eff_by_id: dict[str, float] = {}
    for gene in targets:
        cands = [c for c in all_candidates if c["gene"] == gene]
        eff_by_id.update(score_with_oligoformer(gene, mrna_by_gene[gene], cands, args.cuda))
    for c in all_candidates:
        c["OligoFormer_efficacy"] = round(eff_by_id[c["id"]], 6) if c["id"] in eff_by_id else None

    scored = [c for c in all_candidates if c["OligoFormer_efficacy"] is not None]
    print(f"  打分成功 {len(scored)}/{len(all_candidates)} 条")

    # ---- 3. 逐基因排名（efficacy 高者优先）----
    by_gene: dict[str, list[dict]] = {}
    for c in scored:
        by_gene.setdefault(c["gene"], []).append(c)
    picked = []
    for gene, items in by_gene.items():
        ranked = sorted(items, key=lambda x: -x["OligoFormer_efficacy"])
        for rank, it in enumerate(ranked, 1):
            it["OligoFormer_rank_in_gene"] = rank
        picked.extend(ranked[:args.top_k])
        with open(OUT_DIR / f"of_scores_{gene}.json", "w") as f:
            json.dump({
                "gene": gene, "accession": targets[gene],
                "n_candidates": len(items),
                "direction": DESIGNER_DIRECTION,
                "length_convention": "21 nt 引导链取 guide[0:19] 作 19 nt 双链核心送入 OligoFormer",
                "note": "efficacy 是设计器信号，不是检测器值，不参与 r_1..r_4。",
                "candidates": ranked,
            }, f, indent=2, ensure_ascii=False)
    print(f"  每基因取 efficacy 最高的 {args.top_k} 条，共 {len(picked)} 条")

    if args.skip_channels:
        print("\n--skip-channels：只打分，不跑通道。")
        return

    # ---- 4. 真实脱靶通道（与 RNA-FM 那条**同函数**）----
    print(f"\n=== 真实脱靶通道：{len(picked)} 条 ===")
    ca = RFD.run_channel_a(picked)
    cb = RFD.run_channel_b(picked, args.tissue)

    # ---- 5. 级联记录 ----
    entries = []
    for c in picked:
        cae, cbe = ca[c["id"]], cb[c["id"]]
        r_1 = cae["offtarget_count"]
        r_2, n8, n_genes = RFD.compute_L3_silencing_severity(cbe)
        r_3, top_genes = RFD.compute_L5_from_channel_b(cbe)
        r_4, _hits, tlr_detail = RFD.score_tlr_immunostimulation(c["guide"])
        entries.append({
            "designer": DESIGNER_LABEL,
            "designer_kind": "专用 siRNA 效力预测器（在真实 knockdown 数据上训练）",
            "candidate_id": c["id"], "gene": c["gene"], "accession": c["accession"],
            "target_pos": c["target_pos"],
            "guide_strand": c["guide"], "of_core_19nt": c["of_core_19nt"],
            "seed_region": cbe["seed_7mer"],
            "designer_score": c["OligoFormer_efficacy"],
            "designer_score_direction": DESIGNER_DIRECTION,
            "designer_score_unit": "OligoFormer 预测的 knockdown efficacy；设计器信号，非检测器值",
            "designer_rank_in_gene": c["OligoFormer_rank_in_gene"],
            "r_1": r_1, "r_2": r_2, "r_3": r_3, "r_4": r_4,
            "r_1_unit": "Channel A 转录组脱靶比对条数（计数，非 [0,1]）",
            "r_2_unit": "带种子互补位点的基因数 / 扫描基因总数（比例）",
            "r_3_unit": "DepMap Chronos 连续必需性，可 > 1",
            "r_4_unit": "引导链 U 比例",
            "cascade_matrix_summary": RFD.summarize_cascade_matrix(r_1, r_2, r_3, r_4),
            "L1L2_details": {k: cae[k] for k in
                             ("channel_a_hits", "offtarget_count", "min_mismatches")},
            "L3_details": {"n_genes_with_site": n_genes,
                           "n_genes_scanned": cbe["n_genes_scanned"],
                           "n_sites_8mer": n8,
                           "n_sites_7mer_m8": cbe["n_sites_7mer_m8"]},
            "L5_details": {"most_essential_gene": cbe["most_essential_gene"],
                           "n_common_essential_hit": cbe["n_common_essential_hit"],
                           "top_genes": top_genes},
            "L6b_details": tlr_detail,
            "tissue": cbe["tissue"],
        })

    with open(OUT_DIR / "cascade_results.json", "w") as f:
        json.dump({
            "task": "T_siRNA · OligoFormer（专用 siRNA 设计器）",
            "model": "OligoFormer (Bai et al., Bioinformatics 2024) · model/best_model.pth",
            "n": len(entries),
            "designer_direction": DESIGNER_DIRECTION,
            "detectors": "与 sirna_cascade_pipeline / rna_fm_designer 同函数同定义（import，非重写）",
            "channels": {"A": "Bowtie2 vs 人类转录组",
                         "B": f"19189 条真实 3'UTR，组织 {args.tissue}"},
            "caveats": [
                "OligoFormer 同样不是生成模型，是判别式 oracle；加它解决的是「排序信号与任务无关」，不是「被测的不是生成能力」。",
                "21 nt 引导链取 guide[0:19] 送入（19 nt 双链核心 + 2 nt 3′ 突出端的标准 siRNA 口径）。",
                "未使用 OligoFormer 自带的 off-target / toxicity 模块——脱靶判据必须与 RNA-FM 那两档逐字相同。",
            ],
            "all_results": entries,
        }, f, indent=2, ensure_ascii=False)
    print(f"\n写出 {OUT_DIR/'cascade_results.json'}（{len(entries)} 条）")


if __name__ == "__main__":
    main()
