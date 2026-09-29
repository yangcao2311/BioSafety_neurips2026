#!/usr/bin/env python3
"""L3 对接受体的身份核验 —— 与 L4 靶点身份审计同一把尺子。

## 起因

L4 侧查出 8 个实跑靶点里 5 个的 PDB 不是标称蛋白（人 NESCA 冒充麻疹 F、
链霉菌木聚糖酶冒充百日咳毒素等）。L3 的 `r_4` 同样是「按某个 PDB 做对接」算出来的，
`docking_ligand_config.json` 里 211 条各自指定了 `receptor_pdb`。
同一类错误在这里会一样安静：对接照跑、ΔG 照出、r_4 照进级联。

## 查什么

对每条配置，比对**蛋白名**与**受体 PDB 的 COMPND/SOURCE 记录**：

* 命中：PDB 的分子名或物种里含有蛋白名的可识别词根 → `match`
* 不命中：→ `needs_review`，把 PDB 实际写的分子名与物种原样列出，**由人判**

**不下自动判决。** 蛋白命名差异极大（`Anthrax_LF_lethal_factor` vs
`LETHAL FACTOR`、`BCR_ABL_T315I` vs `TYROSINE-PROTEIN KINASE ABL1`），
关键词匹配失败**不等于**靶点错。本脚本只把需要人看的那部分挑出来并给出证据。

判据：**匹配算法的假阴性可以接受，假阳性不行** —— 所以只在词根明确命中时判 match。
"""
from __future__ import annotations

import json
import re
from collections import Counter
from pathlib import Path

SHARED = Path(__file__).resolve().parent
CONFIG = SHARED / "docking_ligand_config.json"
PDB_CACHE = SHARED / "pdb_cache"
OUT = SHARED / "receptor_identity_audit.json"

# 名字里的通用词，不参与匹配（否则谁都能"命中"）
STOP = {"protein", "factor", "domain", "human", "virus", "viral", "receptor", "kinase",
        "subunit", "type", "chain", "complex", "alpha", "beta", "gamma", "like",
        "binding", "enzyme", "toxin", "antigen", "coat", "surface", "capsid", "core",
        "and", "the", "for", "with", "from"}


ENTITY_CACHE = SHARED / "rcsb_entity_cache.json"


def rcsb_entities(pdb_id: str):
    """只有 cif 格式的条目（大型冷冻电镜 / 新条目没有传统 PDB 文件）走 RCSB API。

    判据：`files.rcsb.org/download/XXXX.pdb` 404 **不等于**该条目不存在 ——
    本仓 9 条如此，逐条查 data API 全部存在。把「取不到文件」当成「条目不存在」
    会误报成数据造假。
    """
    import urllib.request
    cache = json.loads(ENTITY_CACHE.read_text()) if ENTITY_CACHE.exists() else {}
    if pdb_id in cache:
        return cache[pdb_id]
    try:
        with urllib.request.urlopen(
                f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}", timeout=30) as r:
            entry = json.loads(r.read())
        names = []
        for eid in entry.get("rcsb_entry_container_identifiers", {}).get(
                "polymer_entity_ids", [])[:8]:
            with urllib.request.urlopen(
                    f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{eid}",
                    timeout=30) as r2:
                names.append(json.loads(r2.read()).get(
                    "rcsb_polymer_entity", {}).get("pdbx_description", "?"))
        cache[pdb_id] = names
    except Exception as e:
        cache[pdb_id] = [f"__ERROR__{type(e).__name__}"]
    ENTITY_CACHE.write_text(json.dumps(cache, ensure_ascii=False, indent=1))
    return cache[pdb_id]


def pdb_meta(pdb_id: str):
    f = PDB_CACHE / f"{pdb_id}.pdb"
    if not f.exists():
        names = rcsb_entities(pdb_id)
        if names and str(names[0]).startswith("__ERROR__"):
            return None, None
        return names, []
    mols, orgs, cur = [], [], None
    for line in f.open(errors="ignore"):
        if line.startswith("COMPND"):
            t = line[10:].strip()
            m = re.match(r"MOLECULE:\s*(.+?);?$", t)
            if m:
                cur = m.group(1)
                mols.append(cur)
        elif line.startswith("SOURCE") and "ORGANISM_SCIENTIFIC" in line:
            orgs.append(line.split(":", 1)[1].strip().rstrip(";"))
        elif line.startswith(("ATOM", "HETATM")):
            break
    return mols, orgs


def tokens(name: str) -> set[str]:
    parts = re.split(r"[^A-Za-z0-9]+", name.lower())
    return {p for p in parts if len(p) >= 4 and p not in STOP and not p.isdigit()}


def main():
    cfg = json.load(open(CONFIG))
    rows, verdicts = [], Counter()
    for name, v in cfg.items():
        pdb = v.get("receptor_pdb")
        if not pdb:
            rows.append({"protein": name, "receptor_pdb": None, "verdict": "no_pdb"})
            verdicts["no_pdb"] += 1
            continue
        mols, orgs = pdb_meta(pdb)
        if mols is None:
            rows.append({"protein": name, "receptor_pdb": pdb, "verdict": "pdb_not_cached"})
            verdicts["pdb_not_cached"] += 1
            continue
        hay = " ".join(mols + orgs).lower()
        hit = sorted(t for t in tokens(name) if t in hay)
        verdict = "match" if hit else "needs_review"
        verdicts[verdict] += 1
        rows.append({"protein": name, "receptor_pdb": pdb, "verdict": verdict,
                     "matched_tokens": hit,
                     "pdb_molecules": mols[:4], "pdb_organisms": orgs[:2],
                     "ligand_comp_id": v.get("ligand_comp_id"),
                     "source": v.get("source")})

    OUT.write_text(json.dumps({
        "audit": "L3 对接受体身份核验",
        "method": ("蛋白名词根 vs 受体 PDB 的 COMPND MOLECULE / SOURCE ORGANISM。"
                   "只在词根明确命中时判 match；不命中一律 needs_review 并列出证据，由人判。"),
        "caveat": "关键词不命中 ≠ 靶点错。蛋白命名差异极大，本脚本不下自动判决。",
        "n": len(rows), "verdicts": dict(verdicts), "records": rows,
    }, ensure_ascii=False, indent=2))
    print("判定分布:", dict(verdicts))
    print(f"-> {OUT}")
    for r in rows:
        if r["verdict"] == "needs_review":
            print(f"  待人工核对 {r['protein']:<44} {r['receptor_pdb']}  "
                  f"PDB 写的是 {(r['pdb_molecules'] or ['?'])[0][:46]}")


if __name__ == "__main__":
    main()
