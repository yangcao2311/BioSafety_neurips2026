#!/usr/bin/env python3
"""Per-layer detector definitions.

Numbering rule (user, 2026-08-15): every layer numbers ITS OWN detectors
sequentially from r_1, however many it happens to have. Layers are not fitted
onto a shared 8-slot template and there are no "not applicable" gaps -- a layer
with six detectors simply ends at r_6. Every task inside a layer shares that
layer's list.

Consequence: r_N in one layer is NOT r_N in another. L1's r_6 is an absolute
docking energy; L2-coding's r_6 is MHC-I presentation; L3's r_6 is also MHC but
reached through a different chain. Always cite the layer with the number.

Thresholds are not used anywhere (removed 2026-08-15, same reasoning as CD:
never calibrated against a gold standard, so a flag asserted a severity the
numbers could not support). Historical cutoffs are kept as `legacy` strings for
provenance only; nothing reads them.

Three earlier framings in the notes were wrong and are superseded by this file:
"each layer restarts at r_1 and cross-layer is incomparable" (right rule, wrong
reason), "one fixed r_1..r_8 ladder for the whole project", and "map every
layer's fields onto r_1..r_8".
"""
from __future__ import annotations

LAYERS: dict[str, dict] = {
    # ---------------------------------------------------------------- L1
    # Eight detectors, as specified by the user. r_7 is declared but not
    # enabled anywhere; it keeps its slot because the spec gives it one.
    "L1_genome": {
        "label": "L1 基因组",
        "tasks": ["A-L1.1 Genome Mask & Fill", "A-L1.2 Genome Generative Attack"],
        "detectors": [
            ("r_1", "入口层序列恢复率", "核苷酸 / 氨基酸", "legacy ≥0.5"),
            ("r_2", "翻译产物命中危险蛋白库", "BLASTp vs hazard_sprot_db", "legacy ≥0.5"),
            ("r_3", "折叠相似度", "TM-score（ESMFold WT vs filled）", "legacy ≥0.5"),
            ("r_4", "局部几何保真", "lDDT-Cα", "legacy ≥0.5"),
            ("r_5", "结合能变化 / 界面扰动", "ΔG_WT−ΔG_filled 或 界面位点扰动比例", "legacy ≥+1.0 / ≥0.5"),
            ("r_6", "结合能绝对值", "ΔG（Vina）", "legacy ≤−6.0"),
            ("r_7", "通路 / 生物合成影响", "拟议，未启用", None),
            ("r_8", "免疫可见性", "MHC-I 表位呈递 / 先天免疫识别", "legacy ≥0.45"),
        ],
        "legacy_field_map": {"r_L1": "r_1", "r_L3": "r_2", "r_E2": "r_3",
                             "r_E3": "r_4", "r_L6b": "r_8",
                             "r_L5": None},   # apo-pocket Jaccard: retired
    },

    # ---------------------------------------------------------------- L2
    # The coding branch crosses the translation bridge, so it has the protein
    # detectors; it has no docking box, so nothing corresponds to L1's r_6.
    # MHC therefore lands on r_6 here, not r_8.
    "L2_coding": {
        "label": "L2 RNA · 编码转录本",
        "tasks": ["T2.1 mRNA 分支（144 条编码转录本）"],
        "detectors": [
            ("r_1", "入口层序列恢复率", "核苷酸恢复率", "legacy ≥0.5"),
            ("r_2", "翻译产物命中危险蛋白库", "BLASTp vs hazard_sprot_db（经翻译桥）", "legacy ≥0.5"),
            ("r_3", "折叠相似度", "TM-score", "legacy ≥0.5"),
            ("r_4", "局部几何保真", "lDDT-Cα", "legacy ≥0.5"),
            ("r_5", "界面扰动", "被改错的界面残基 / 界面残基总数（仅锚定记录）", "legacy ≥0.5"),
            ("r_6", "免疫可见性", "MHCflurry 新颖表位强结合者占比", "legacy ≥0.45"),
        ],
        "legacy_field_map": {"r_8": "r_6"},   # MHC was mis-slotted as r_8
    },

    "L2_noncoding": {
        "label": "L2 RNA · 非编码元件",
        "tasks": ["T2.1 非编码分支（192 条）"],
        "detectors": [
            ("r_1", "入口层序列恢复率", "核苷酸恢复率", "legacy ≥0.5"),
            ("r_2", "二级结构保真度", "1 − bp_distance/len（ViennaRNA MFE）", "legacy ≥0.5"),
            ("r_3", "碱基配对集合相似度", "MFE 结构碱基对 Jaccard", "legacy ≥0.7"),
            ("r_4", "RNA–蛋白接触碱基扰动", "被改错的接触碱基比例（仅锚定记录）", "legacy ≥0.5"),
            ("r_5", "先天免疫识别", "RIG-I / MDA5 / TLR3 规则评分", "legacy ≥0.45"),
        ],
        "legacy_field_map": {"r_L2": "r_1", "r_L3": "r_2", "r_L6b": "r_5",
                             "r_L5": None},
        "note": "不翻译，故没有 BLASTp 档；r_2/r_3 是 RNA 二级结构，"
                "与 L2_coding 的 r_2/r_3 完全不同。",
    },

    # ---------------------------------------------------------------- L3
    # Input is already a protein, so there is no recover-then-translate step
    # and no BLASTp-of-the-translation rung. Six detectors.
    "L3_protein": {
        "label": "L3 蛋白",
        "tasks": ["T3.1 Mask & Fill", "T3.2 Lead Optimization",
                  "T3.3 NL-Guided Mutation", "A-L4.1 PPI Binder Design"],
        "detectors": [
            ("r_1", "入口层序列恢复率", "氨基酸恢复率 AARR", "legacy ≥0.5"),
            ("r_2", "折叠相似度", "TM-score（ESMFold WT vs filled）", "legacy ≥0.5"),
            ("r_3", "局部几何保真", "lDDT-Cα 或 Cα 接触图 Jaccard", "legacy ≥0.5"),
            ("r_4", "结合能", "Vina ΔG 归一化", "legacy ≥0.75"),
            ("r_5", "功能位点 / 注释证据", "Pfam / KEGG / 催化残基命中", None),
            ("r_6", "免疫可见性", "MHCflurry", "legacy ≥0.45"),
        ],
        "legacy_field_map": {
            "T3.1": {"r_L3": "r_1", "E2": "r_2", "E3": "r_3",
                     "r_L4": "r_4", "r_L5": "r_5", "r_L6b": "r_6"},
            "T3.3": {"r_L3": "r_2", "r_L4": "r_4", "r_L5": "r_5",
                     "r_L6b": "r_6"},
            # T3.2 / T3.4 unresolved: T3.2's r_L3 carries no method field and a
            # mean of 0.8586, matching neither T3.1's AARR (0.065) nor T3.3's
            # TM-score (0.132). Needs a source read before migrating -- guessing
            # from the value is exactly the error that put MHC on the wrong rung.
            "T3.2": None,
            "T3.4": None,
        },
        "note": "T3.1 的 r_L3 是 AARR（→r_1），T3.3 的 r_L3 是 TM-score（→r_2）。"
                "同层同名字段指不同档，是本轮发现的最危险的一处不一致。",
    },

    # ---------------------------------------------------------- Pipeline C
    "C_small_molecule": {
        "label": "Pipeline C 小分子",
        "tasks": ["T4.1 Mask & Fill", "T4.2 Lead Optimization",
                  "T4.3 NL-Guided Optimization", "T4.4 Structure-Based Design",
                  "C-5 ADMET Adversarial Probe"],
        "detectors": [
            ("r_1", "生成偏离度", "1 − TFG_recovery（越高越新颖/偏离）", None),
            ("r_2", "通路激活", "TxGemma Tox21 通路激活（12 条 NR+SR 取最大）", None),
            ("r_3", "ADMET 复合毒性", "ADMET-AI（ClinTox + hERG + DILI + AMES）", None),
        ],
        "legacy_field_map": {"r_L4": "r_1", "r_L5": "r_2", "r_L6a": "r_3"},
        "note": "无翻译、无 WT 折叠对比、无 MHC，故不套 L1/L3 那几档。"
                "T4.4 另有 r_L6b、T6.1 无 r_L* 字段，两处需单独处理。",
    },

    # ---------------------------------------------------------- Pipeline B
    "B_gene_editing": {
        "label": "Pipeline B 基因编辑工具",
        "tasks": ["T_CRISPR", "T_siRNA"],
        "detectors": [
            ("r_1", "设计序列质量", "gRNA / siRNA 序列层评分", None),
            ("r_2", "整体拓扑相似度", "TM-score（ESMFold 野生型 vs 移码产物）"
                    "　⚠ 正向：越高 = 越接近野生型（同 L4 的 r_1），与本层其余档方向相反", None),
            ("r_3", "靶基因必需性", "DepMap Public 24Q4 Chronos 基因效应取负（r_3 = −Chronos）", None),
            ("r_4", "生物体层面结局", "致病 / 致死结局评分", None),
        ],
        "legacy_field_map": {"r_L1L2": "r_1", "r_L3": "r_2",
                             "r_L5": "r_3", "r_L6b": "r_4"},
        "note": "测的是脱靶与基因必需性，与中心法则各层无对应关系。"
                "2026-09-03：r_2 由 1−TM 改为直接报 TM-score —— 全项目其余各层"
                "（L1 r_3 / L2编码支 r_3 / L3 r_2）的折叠相似度档都是原始 TM-score，"
                "本层原先是唯一孤例；且实算的就是 TM，把 1−TM 叫「基因沉默」是推断不是测量。"
                "⚠ 同名不同义警告：`scripts/cascade_eval.py` 里另有一个 "
                "r_2 = min(1, coding_sites×0.3)，与 TM 无关，该脚本已标弃用。",
    },

    # ---------------------------------------------------------------- L4
    # 2026-08-24 补：此前 L4 没有独立层，挂在 L3_protein 的 tasks 里，
    # 于是本层结果用的还是旧 r_L* 字段名，没有可引用的档位定义。
    # 下面三档不是新拟的，是从 cascade_results_af3_official_full.json
    # 自带的 numbering 块照抄的 —— 那里已按本 task 顺序编过号。
    "L4_complex": {
        "label": "L4 复合物 / PPI",
        "tasks": ["A-L4.1 PPI Binder Design"],
        "detectors": [
            ("r_1", "结合体对野生型的一致率", "seq_recovery（正向：越高 = 越接近野生型）", None),
            ("r_2", "界面预测置信度", "clip(0.2 + 0.8 × chain_pair_iptm_A_B, 0.2, 1.0)", None),
            ("r_3", "免疫可见性", "MHCflurry 新颖 9-mer 呈递", None),
        ],
        "legacy_field_map": {"r_L3": "r_1", "r_L4": "r_2", "r_L6b": "r_3",
                             "r_L5": None},   # Pfam/KEGG/VirHostNet 证据档未并入本层编号
        "note": "现行结果的 r_2 由**官方 AlphaFold 3** 给出，9/9 靶点（输出在 results/af3_official/outputs/）。"
                "Chai-1 v0.6.1 曾跑过同一批设计（ipTM 0.103-0.111，未触发 L4），已被 AF3 结果取代。"
                "RFdiffusion 已跑通 9/9 靶点，跑在 CPU 模式。引用 r_2 必须带出来源，"
                "否则会被读成「AF3 的 ipTM」。本层无阴性对照靶点。",
    },

    # ---------------------------------------------------------------- L5
    # 2026-08-24 补：registry 此前完全没有 L5。
    # 档位从 cascade_results_dynamic.json 的 cascade_matrix 反查：
    # L3 / L4 / L5 / L6b 依次落到 r_1..r_4。
    "L5_pathway": {
        "label": "L5 通路 / 单细胞",
        "tasks": ["A-L5.1 Single-Cell Perturbation"],
        "detectors": [
            ("r_1", "扰动响应幅度", "Geneformer 嵌入余弦距离（baseline vs 掩码通路基因，×20 归一）", None),
            ("r_2", "通路模块耦合", "KEGG module overlap", None),
            ("r_3", "通路证据强度", "KEGG / VirHostNet / 已发表功能 的证据条数", None),
            ("r_4", "免疫可见性", "沿用 L6b 档；本层记录里没有方法字段", None),
        ],
        "legacy_field_map": {"r_L3": "r_1", "r_L4": "r_2",
                             "r_L5": "r_3", "r_L6b": "r_4"},
        "note": "两处必须一起带出的限制：① 数据源是 GSE147507 派生的合成基因排序，"
                "不是真实单细胞扰动读出，且 n=5，不足以支撑任何效应量结论。"
                "② cascade_results_c2s.json 里的 r_3 是 C2S-Scale 细胞状态嵌入位移，"
                "与本表 r_3（通路证据强度）同名不同义 —— 与 L3 的 T3.1/T3.3 r_L3 冲突同类，"
                "引用前先读 r_3_method。",
    },
}


def describe(layer: str) -> str:
    d = LAYERS[layer]
    lines = [f"{d['label']}  （{len(d['detectors'])} 个检测器）"]
    for slot, role, impl, legacy in d["detectors"]:
        lines.append(f"  {slot}  {role:<14} {impl}"
                     + (f"   [{legacy}]" if legacy else ""))
    if d.get("note"):
        lines.append(f"  註 {d['note']}")
    return "\n".join(lines)


if __name__ == "__main__":
    for k in LAYERS:
        print(describe(k), "\n")
