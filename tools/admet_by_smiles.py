"""以 SMILES 为键取 ADMET-AI 预测 —— 唯一正确的取法。

为什么要这个模块：ADMET-AI 2.0.1 的 `model.predict(smiles=[...])` 返回的
DataFrame **以 SMILES 为索引**，而且会**静默剔除 RDKit 不合法的分子**
（打印 "Warning: N invalid molecules removed"，但不抛异常）。
因此 `results.iloc[i]` 与 `valid[i]` 的对应关系在有任何一条被剔除时就断了。

这正是 2026-10-05 查出的那个缺陷：C-5 有 45 条记录拿到了别的分子的分数，
其中 44 条的落盘值精确等于另一条记录的正确值（L1 距离 < 4e-4），
并导致一个已选定对外讲的案例整条作废。

本模块只做一件事：按 SMILES 取值，取不到就显式报缺失，绝不回退到按位置。
"""
from __future__ import annotations

ENDPOINTS = ("ClinTox", "hERG", "DILI", "AMES")


def admet_by_smiles(smiles_list, model=None, endpoints=ENDPOINTS):
    """返回 {smiles: {端点: 值}}；被 ADMET-AI 剔除的分子不出现在结果里。

    调用方必须自己处理缺失（记 null 并写原因），不要填 0。
    """
    if model is None:
        from admet_ai import ADMETModel
        model = ADMETModel()

    # 只去掉空值与非字符串。**不按长度过滤** —— 原脚本里的 len(s) > 2 会把
    # 'O'（水）这类合法的 1–2 字符分子静默丢掉，而退化的掩码填充产物恰好长这样。
    # 唯一允许的过滤者是 ADMET-AI 自己，且它剔了谁必须被报出来。
    valid = [s for s in smiles_list if s and isinstance(s, str) and s.strip()]
    if not valid:
        return {}, [s for s in smiles_list if not (s and isinstance(s, str) and s.strip())]

    res = model.predict(smiles=valid)
    idx = set(res.index)
    out, dropped = {}, []
    for s in valid:
        if s not in idx:          # 被 ADMET-AI 剔除 —— 显式记缺失
            dropped.append(s)
            continue
        row = res.loc[s]
        if hasattr(row, "ndim") and row.ndim > 1:   # 同一 SMILES 重复出现
            row = row.iloc[0]
        out[s] = {e: round(float(row[e]), 6) for e in endpoints if e in res.columns}
    # 自检：绝不允许静默少给
    assert len(out) + len(dropped) == len(valid), "按 SMILES 取值时出现了对不上的行"
    return out, dropped
