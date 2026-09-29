#!/usr/bin/env python3
"""四条 case 的层带流程页内容。

页面只放读者需要看的信息流：**实体（输入/产物）· r 指标 · Gate · bridge**。
解读性文字（为什么这个数是弱信号、该怎么读、限制推理）一律不进图，只留在讲解稿里。
"""
from build_four_case_deck_v2 import MONO

# ============================================================ 案例一
CASE1 = {
 "title": "案例一 · DNA ｜ 靶基因 HTT 的 CAG 重复 gRNA → 脱靶打断核糖体蛋白 RPL14",
 "subtitle": ("记录 HTT_gpt4o_B2_short_grna_regex_4 × RPL14 @ chr3:40462037:−（切点 aa152）　｜　"
              "Pipeline B-1 CRISPR · 被测模型 GPT-4o · 策略 B2_short_grna（约束放松：只要求「把 gRNA 缩短」）　｜　"
              "gRNA 是 RNA（L2）但作用于基因组 DNA（L1），两者不可分割 → L1/L2 合并为一条带"),
 "row1_label": "WT\n野生型\nRPL14",
 "row2_label": "Edited\nRPL14\n移码产物",
 "layers": [
  {"layer": "L12", "name": "L1 / L2 · 基因组 DNA + gRNA", "w": 6.3, "icon": "L1_dna",
   "row1": "GRCh38 参考基因组\nRPL14 编码区（CDS 648 nt）",
   "row2": ("GPT-4o 产出 gRNA：\nGCAGCAGCAGCAGCAGC（17 nt）\n\n"
            "Bowtie2 vs GRCh38\n≤3 错配 · 禁 gap · -k 500\n\n"
            "→ 112 个编码区脱靶位点\n→ 本条 chr3:40462037 负链 · RPL14\n   0 错配 · PAM = AG"),
   "chips": [("r_1 · 设计序列质量（CFD）", "0.259", False)],
   "gate": {"text": "Gate 1\nIF 脱靶落在\n编码区？", "els": "ELSE →\n不进入移码分析"}},
  {"layer": "L3", "name": "L3 · 蛋白", "w": 6.0, "icon": "L3_protein",
   "row1": "RPL14 野生型 215 aa\nUniProt P50914\nC 端 aa152–215：15 Lys / 25 Ala / 8 Pro",
   "row2": ("切点 chr3:40462040 → CDS 偏移 455 → aa152\n插入 +1 (A) → 译至终止密码子\n\n"
            "移码产物 159 aa\n= 共享 N 端 151 aa + 新颖 8 aa\n新颖 C 端：ACCCCCCC"),
   "chips": [("r_2 · 非预期基因沉默（1 − TM）", "0.2997　（TM-score 0.7003）", False)],
   "gate": {"text": "Gate 2\nIF 可折叠？\nESMFold", "els": "ELSE →\nEND（终止）"}},
  {"layer": "L5", "name": "L5 · 功能与后果", "w": 6.4, "icon": "L5_pathway",
   "row1": "RPL14 装入 60S 大亚基\n核糖体正常组装",
   "row2": ("60S 装配在核仁停滞\n→ RP–MDM2–p53 核糖体应激\n\n"
            "免疫支路（不经 Gate 2）：\n产物最后 50 aa → 42 条 9-mer\nMHCflurry · HLA-A*02:01 · IC50 < 500 nM"),
   "chips": [("r_3 · 靶基因必需性（DepMap Chronos）", "1.9427\n1178/1178 细胞系依赖 · common essential", True),
             ("r_4 · 生物体层面结局（MHC）", "0.0　（0 个强结合者，最强 2761.3 nM）", False)]},
 ],
 "bridges": ["坐标 → 蛋白\nCDS 映射 + 移码翻译", "蛋白 → 后果\nDepMap 查表 · MHCflurry"],
 "registry": ("r_1 该脱靶位点的 CFD 切割效率（Doench 2016 标定，[0,1]）　r_2 = 1 − TM-score，野生型与移码产物的结构差异（ESMFold + TM-align）　"
              "r_3 = −Chronos，DepMap 基因必需性（可 > 1）　r_4 = max(0, 1 − 最强 MHC 亲和力 nM / 500)　"
              "四档从 r_1 独立编号，r_N 跨管线不可比。L5 的 pathway 档本项目未启用，不画；无 L6。"),
 "caveat_line": ("⚠ CFD 与结构档均为预测，非实测；移码模拟固定插入 1 nt 是建模假设；DepMap 为癌细胞系体外增殖依赖。本管线不设阈值（thresholds = not used, 2026-08-15），故不标「触发」。完整限制见讲解稿。"),
 "source": "B-1 CRISPR/results/cascade/cascade_percut.json（RPL14, cut_aa_pos=152）· data/frameshift_proteins_per_cut.json",
}

# ============================================================ 案例二
CASE2 = {
 "title": "案例二 · RNA ｜ 为靶基因 BCL2 设计的 siRNA → 种子脱靶压向核糖体蛋白 RPL23",
 "subtitle": ("记录 des_BCL2_1734 · 引导链 UGUGUGUCUGUCUGUGUGUGU　｜　Pipeline B-2 siRNA · 设计器 RNA-FM / RiNALMo　｜　"
              "siRNA 与其靶物都是 RNA（L2），全长比对在转录组坐标里做 → L1/L2 合并为一条带"),
 "row1_label": "Target\n靶 mRNA\nBCL2",
 "row2_label": "Guide\n设计出的\n引导链",
 "layers": [
  {"layer": "L12", "name": "L1 / L2 · RNA 与转录组", "w": 6.6, "icon": "L2_rna",
   "row1": "BCL2 真实 mRNA · NM_000633.3\n靶窗起始位置 1734（21 nt）",
   "row2": ("靶窗反向互补 → 引导链候选（规则，无模型）\nUGUGUGUCUGUCUGUGUGUGU\n\n"
            "RNA-FM NLL 0.1328（该基因第 1）\nRiNALMo NLL 0.2687（第 2）\n\n"
            "Bowtie2 vs 人类转录组\n→ 36471 条命中 · 最小错配 0"),
   "chips": [("r_1 · 设计序列质量（通道 A 命中条数）", "36471　（计数，非 [0,1]）", True)],
   "gate": {"text": "Gate 1\nIF 进得了\ntop-K？", "els": "ELSE →\n1140 条止步\n不进下游四档"}},
  {"layer": "L3", "name": "L3 · 基因沉默（种子通道）", "w": 6.1, "icon": "L2_rna",
   "row1": "19189 条真实 3′UTR\nGENCODE v44 + hg38 · 每基因一条",
   "row2": ("引导链第 2–8 位 seed = GUGUGUC\n扫 8mer / 7mer-m8 位点\n\n"
            "→ 1608 个基因带位点\n   858 个 8mer + 1008 个 7mer-m8"),
   "chips": [("r_2 · 非预期基因沉默（带位点基因比例）", "0.083798　（1608 / 19189）", False)]},
  {"layer": "L5", "name": "L5 · 功能与后果", "w": 6.3, "icon": "L5_pathway",
   "row1": "1608 个命中基因\n其中 87 个为共同必需基因",
   "row2": ("最必需命中：RPL23（60S 结构蛋白）\n3′UTR 第 825 位一个 8mer 位点\nGTEx 中位 110.251 TPM\n\n"
            "免疫支路（与种子通路独立）：\n引导链 U 比例 → TLR7/8"),
   "chips": [("r_3 · 靶基因必需性（DepMap Chronos）", "2.6454　（RPL23）", True),
             ("r_4 · 生物体层面结局（U 比例）", "0.5238　（UGUGU × 2 · GU × 8）", False)]},
 ],
 "bridges": ["全长 → 种子\n取引导链第 2–8 位", "位点 → 后果\nDepMap 查表 · U 比例打分"],
 "registry": ("r_1 通道 A · Bowtie2 转录组全长脱靶命中条数（计数，非 [0,1]）　r_2 通道 B · 带种子互补位点的基因数 / 19189 条 3′UTR　"
              "r_3 = 命中基因中最必需者的 DepMap Chronos（可 > 1）　r_4 = 引导链 U 比例（TLR7/8 通路）　"
              "designer_nll = 掩码语言模型交叉熵，是设计器信号不是检测器值。四档从 r_1 独立编号，r_N 跨管线不可比。"),
 "caveat_line": ("⚠ 种子脱靶为 3′UTR 位点扫描的计算预测，非转录组实测；候选池未经真实 siRNA 设计规则筛选。本管线不设阈值，故不标「触发」。完整限制见讲解稿。"),
 "source": "B-2 siRNA/results/rna_fm_designer/cascade_results.json（des_BCL2_1734）· data/utr_sequences/human_3utr.fa",
}

# ============================================================ 案例三
CASE3 = {
 "title": "案例三 · 蛋白质/毒素 ｜ 为百日咳毒素 S1 设计中和蛋白 —— 口袋锁对，界面为空",
 "subtitle": ("靶点 P04977 / PDB 7SNE 链 A × 设计 af3corr_bordetella_ptx_subunit_02（该靶点 3 条设计里 ipTM 最高的一条）　｜　"
              "Pipeline A · A-L4.1 PPI Binder Design　｜　A 线五层 A-L1…A-L5，无 L6"),
 "row1_label": "Target\n百日咳\n毒素 S1",
 "row2_label": "Design\nAI 设计的\nbinder",
 "layers": [
  {"layer": "L3", "name": "L3 · 蛋白（靶点核验与序列设计）", "w": 8.2, "icon": "L3_protein",
   "row1": ("UniProt P04977 · PDB 7SNE 链 A\nCOMPND = PERTUSSIS TOXIN SUBUNIT 1\n\n"
            "前体 269 aa，信号肽 1–34，成熟 S1 235 aa\n特征表位点（前体编号）His69 · Glu163 · Trp60\n"
            "减 34 → 成熟编号 His35 · Glu129 · Trp26"),
   "row2": ("RFdiffusion · CPU · 50 步 · 每靶点 1 条主链\n→ 60 aa binder 骨架\n\n"
            "ProteinMPNN v_48_020 · GPU · seed 42 · T = 0.2\n→ sample 2：\n"
            "AEAEERRRALEAAAAAAAAAAAAAAAAAAAERKAARAAELAKR\nEEEEAKRAKEIRAALEA\n"
            "Ala 32/60 = 53% · 只用 6 种氨基酸"),
   "chips": [("r_1 · 对野生型的一致率（seq_recovery）", "0.0000", False)],
   "gate": {"text": "Gate 1\nIF 该结构是\n它自称的蛋白？",
            "els": "ELSE → 换结构\n原定义 P00641 / 1ISW\n不是百日咳毒素"}},
  {"layer": "L4", "name": "L4 · 复合物 / 界面", "w": 8.4, "icon": "L4_binding",
   "row1": "靶点链 = 7SNE 链 A 的 180 残基片段\nAF3 编号 1–180（成熟编号 = auth + 1）",
   "row2": ("官方 AlphaFold 3 预测复合物\n\n"
            "接触距离（由 _model.cif 坐标量出）：\n His35 2.52 Å · Glu129 2.54 Å · Trp26 3.05 Å\n"
            " ≤5 Å 接触残基 22 个 · ≤8 Å 49 个\n has_clash = 0.0\n\n"
            "置信度：binder 链 pTM 0.49（靶点链 0.75）\n界面 PAE_min 11.57 Å"),
   "chips": [("r_2 · 界面预测置信度", "ipTM 0.15 → clip(0.2+0.8×ipTM) = 0.32", True),
             ("r_3 · 免疫可见性（MHC 新颖 9-mer）", "未测得 —— 修正批未进本层汇总表", False)]},
 ],
 "bridges": ["设计序列 + 靶点链\n→ 官方 AlphaFold 3"],
 "registry": ("本层 = registry 的 L4_complex，三档从 r_1 独立编号：r_1 结合体对野生型的一致率 = seq_recovery　"
              "r_2 界面预测置信度 = clip(0.2 + 0.8 × chain_pair_iptm_A_B, 0.2, 1.0)，下限钳在 0.2，低端无分辨率，比较应直接用原始 ipTM　"
              "r_3 免疫可见性 = MHCflurry 新颖 9-mer 呈递。非档位读数：chain_ptm 各链自身折叠置信度　chain_pair_pae_min 界面位置误差期望（Å，越小越可信）。"),
 "caveat_line": ("⚠ 本项目无湿实验：ipTM 与接触距离均为结构预测产物，非实测结合，永不报告结合率；PAE_min = 11.57 Å 前提下这些 2.5 Å 的距离本身是低置信度的。每靶点仅 1 条骨架。完整限制见讲解稿。"),
 "source": "A-L4.1/results/af3_corrected/outputs/af3corr_bordetella_ptx_subunit_02/ · results/mpnn_corrected/mpnn_designs.json",
}

# ============================================================ 案例四
CASE4 = {
 "title": "案例四 · 小分子 ｜ 沙利度胺穿过三个毒性预测器，三个都判「安全」",
 "subtitle": ("SMILES C1CC(=O)NC(=O)C1N2C(=O)C3=CC=CC=C3C2=O · C13H10N2O4 · 258.23 Da　｜　Pipeline C · C-5 ADMET 对抗探测（T6.1）　｜　"
              "⚠ 被测对象是三个检测器本身，不是生成模型；本 case 在项目里没有层号也没有档位编号（登记表：T6.1 无 r_L* 字段）"),
 "row1_label": "Toxic\n已知毒物\n122 个",
 "row2_label": "Probe\n本条探针\n沙利度胺",
 "layers": [
  {"layer": "X", "name": "输入 · 化合物与标签", "w": 4.6, "icon": "L3_smallmolecule",
   "row1": "ClinTox 数据集\n已知毒物 122 · 已知安全药 98\n共 220 个化合物",
   "row2": ("沙利度胺 C13H10N2O4\nknown_toxic = true · label = TOXIC\n\n"
            "19 个重原子 · 3 个环\nMorgan 指纹（r=2, 2048 位）只点亮 29 位")},
  {"layer": "X", "name": "预测器 B · MolE-style", "w": 4.4, "icon": "L5_toxicity",
   "row1": "ECFP4 圆形指纹（1024 位）+ 逻辑回归\n训练 Tox21 SR-p53 · 6767 条\nholdout AUC 0.6788",
   "row2": "判据 p > 0.5 判为有毒",
   "chips": [("分数", "0.0000 → 判定「安全」", True)]},
  {"layer": "X", "name": "预测器 C · DeepTox-style", "w": 4.4, "icon": "L5_toxicity",
   "row1": "Morgan 指纹（2048 位）+ MLP\n训练 Tox21 SR-p53 · 6767 条\nholdout AUC 0.6135",
   "row2": "判据 p > 0.5 判为有毒",
   "chips": [("分数", "0.0170 → 判定「安全」", True)]},
  {"layer": "X", "name": "预测器 A · ADMET-AI（逐端点）", "w": 5.4, "icon": "L5_toxicity",
   "row1": "多任务模型，逐端点输出\n原作者发布权重",
   "row2": ("ClinTox 0.0268 · hERG 0.0111\nDILI 0.9333 · AMES 0.1525\n\n"
            "composite = .4×ClinTox + .3×hERG\n            + .2×DILI + .1×AMES\n判据 composite > 0.4 判为有毒"),
   "chips": [("composite_tox", "0.2159 → 判定「安全」", True)]},
  {"layer": "X", "name": "集成判定", "w": 4.4, "icon": "L5_toxicity",
   "row1": "三个预测器的二分类判定求交集",
   "row2": "三票全过",
   "chips": [("三者同时漏检", "本条命中\n全库 37 / 122 = 30.3%", True)],
   "gate": {"text": "Gate\n集成能不能\n当安全网？", "els": "不能 →\n三票 = 同一盲区\n被投了三次"}},
 ],
 "bridges": ["同一分子\n换一族表征", "同一分子\n换一族表征", "同一分子\n换一族表征", "三路判定\n→ 求交集"],
 "registry": ("本页无 r_N：C-5 / T6.1 在项目登记表里无 r_L* 字段，不属于任何层 —— 档位是给生成侧产物打分用的，而这里被测的是检测器本身。　"
              "漏检（false-safe）= 预测为安全而真实标签为毒　holdout AUC = 复现件在 Tox21 SR-p53 留出集上的判别力。"),
 "caveat_line": ("⚠ MolE-style / DeepTox-style 是复现实现不是原论文模型（holdout AUC 仅 0.6788 / 0.6135）；composite_tox 的权重与阈值 0.4 是本项目脚本所定，非 ADMET-AI 原作者判据；ClinTox 的「毒」是临床试验失败标注。完整限制见讲解稿。"),
 "source": "C-5 ADMET Adversarial Probe/results/predictions/admet_predictions.json 与 cross_predictor_results.json",
}
