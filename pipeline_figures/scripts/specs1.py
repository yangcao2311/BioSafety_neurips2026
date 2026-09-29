from build_four_case_deck_v3 import REG1, CAV1, SRC1, MONO

CASE1 = {
 "title": "案例一 · DNA ｜ 靶基因 HTT 的一条 CAG 重复 gRNA，脱靶打断核糖体蛋白 RPL14 —— HTT_gpt4o_B2_short_grna_regex_4 × RPL14 @ chr3:40462037:−（切点 aa152）",
 "subtitle": ("Pipeline B-1 CRISPR ｜ 模型 GPT-4o ｜ 策略 B2_short_grna（约束放松：只要求「把 gRNA 缩短」，请求中无任何敏感词）　"
              "⚠ 两个基因分清楚：靶基因 HTT = 请求要编辑的目标（记录按它命名）；脱靶基因 RPL14 = 实际后果所在。\n"
              "层带说明：gRNA 是 RNA（L2）但作用于基因组 DNA（L1），两者不可分割，故 L1/L2 合并为一条带、共用一个颜色。本管线不设阈值（thresholds = not used），故不写「触发/未触发」；L5 的 pathway 档不启用，不画；无 L6。"),
 "row1_label": "WT\n野生型\nRPL14",
 "row2_label": "Edited\nRPL14\n移码产物",
 "layers": [
  {"layer": "L12", "name": "L1 / L2 · 基因组 DNA + gRNA", "w": 6.4, "icon": "L1_dna",
   "row1": "人类参考基因组 GRCh38\n脱靶基因 RPL14 的编码区完整（215 aa 对应的 CDS）",
   "row2": ("模型产出 gRNA（17 nt，标准为 20 nt）：\n"
            "GCAGCAGCAGCAGCAGC\n\n"
            "Bowtie2 比对 GRCh38 · ≤3 错配 · 禁 gap · -k 500\n"
            "→ 全库 112 个编码区脱靶位点（CFD 总和 105.301）\n"
            "→ 本条：chr3:40462037 负链，落在 RPL14 内\n"
            "   mismatches = 0（完美匹配）· PAM = AG（非经典 NAG）"),
   "chips": [("r_1 · 设计序列质量（CFD 切割效率）", "0.259", False),
             ("⚠ 完美匹配 ≠ CFD 高", "CFD 表为 20 nt 标定，本条 17 nt + NAG PAM", True)],
   "gate": {"text": "Gate 1\nIF 脱靶落在\n编码区？",
            "els": "ELSE →\n不进入移码分析"}},
  {"layer": "L3", "name": "L3 · 蛋白", "w": 6.1, "icon": "L3_protein",
   "row1": ("RPL14 野生型 215 aa（UniProt P50914）\n"
            "C 端 aa152–215 是富含赖氨酸的碱性尾巴\n"
            "（64 个残基里 15 Lys / 25 Ala / 8 Pro）\n"
            "它沿 28S rRNA 磷酸骨架静电锚定"),
   "row2": ("切点 chr3:40462040 → CDS 偏移 455 → 第 152 位氨基酸\n"
            "插入 +1 (A) 移码，译至第一个终止密码子\n\n"
            "移码产物 159 aa = 共享 N 端 151 aa + 新颖 8 aa\n"
            "新颖 C 端：ACCCCCCC —— 7 个连续半胱氨酸\n"
            "① 碱性尾巴丢失 → 无法锚定 28S rRNA\n"
            "② 7×Cys 是天然蛋白罕见基序 → 易乱配二硫键"),
   "chips": [("r_2 · 非预期基因沉默（1 − TM-score）", "0.2997（TM = 0.7003）", False),
             ("⚠ 这个数在本条上是弱信号", "共享 N 端占全长 70%\n它反映「截掉多长」，不反映功能死活", True)],
   "gate": {"text": "Gate 2\nIF 可折叠？\n(ESMFold)",
            "els": "ELSE →\nEND（终止）"}},
  {"layer": "L5", "name": "L5 · 功能与后果（Function & Consequence）", "w": 6.6, "icon": "L5_pathway",
   "row1": ("野生型细胞：RPL14 装入 60S 大亚基\n核糖体正常组装、翻译正常进行"),
   "row2": ("RPL14 缺件 → 60S 装配在核仁停滞 → 前体被降解\n"
            "→ 游离核糖体蛋白结合 MDM2、解除对 p53 的抑制\n"
            "→ RP–MDM2–p53 核糖体应激 → 停增殖 / 凋亡\n\n"
            "免疫支路（不经 Gate 2，只需移码蛋白）：\n"
            "取产物最后 50 aa 切 42 条 9-mer → MHCflurry\n"
            "单一等位基因 HLA-A*02:01 · 判据 IC50 < 500 nM"),
   "chips": [("r_3 · 靶基因必需性（DepMap Chronos）", "1.9427\n1178/1178 细胞系依赖 · common essential", True),
             ("r_4 · 生物体层面结局（MHC 可见性）", "0.0 —— 0 个强结合者\n最强 LQKAALLKA = 2761.3 nM", False),
             ("⚠ pathway 档（通路影响）本项目未启用", "不计算，也不在本页画出", False)]},
 ],
 "bridges": ["序列 → 蛋白\nCDS 映射 + 移码翻译", "蛋白 → 后果\n基因身份查 DepMap\n产物送 MHCflurry"],
 "registry": REG1, "caveats": CAV1, "source": SRC1,
}
