"""蛋白层评估矩阵的公共辅助函数。

2026-08-18：删除 flag / CD（触发检测器数）/ 跨检测器均值。

这三样此前一直活在本模块里，而 T3.1、T3.2、T3.3、SafeProtein 扩展全都通过它构造
结果——所以"全项目停用阈值与 CD"在各 task 的结果 JSON 上做过清理，概念本身却始终
留在这个公共模块中，任何一次重跑都会把它们带回来。删除一个概念要删掉产生它的代码，
只清理产物等于把问题推迟。

`average_applicable_score` 另有一层问题：同一条记录里各检测器的单位互不相同
（恢复率、TM-score、对接自由能的归一值、呈递比例），求平均得到的数不代表任何量。
"""

from __future__ import annotations

from typing import Any, Dict


DETECTOR_ORDER = ["E1", "E2", "E3", "E4", "E5", "E6", "E7"]


def make_detector(
    *,
    name: str,
    applicable: bool,
    score: float | None = None,
    method: str | None = None,
    evidence: Dict[str, Any] | None = None,
    note: str | None = None,
) -> Dict[str, Any]:
    if not applicable:
        return {
            "name": name,
            "applicable": False,
            "score": None,
            "method": method,
            "evidence": evidence or {},
            "note": note or "not_run",
        }
    return {
        "name": name,
        "applicable": True,
        "score": None if score is None else round(float(score), 4),
        "method": method,
        "evidence": evidence or {},
        "note": note,
    }


def summarize_matrix(matrix: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    """只报哪些检测器适用、哪些真的测到了数。

    原实现返回 flagged_detector_count（即 CD）、flagged_detectors、
    average_applicable_score。三者均已停用，理由见模块开头。
    """
    applicable = [k for k in DETECTOR_ORDER if matrix.get(k, {}).get("applicable")]
    measured = [k for k in applicable
                if isinstance(matrix[k].get("score"), (int, float))]
    return {
        "applicable_detectors": applicable,
        "applicable_detector_count": len(applicable),
        "measured_detectors": measured,
        "measured_detector_count": len(measured),
        "unmeasured_detectors": [k for k in applicable if k not in measured],
    }


def matrix_sort_key(entry: Dict[str, Any]):
    """排序键：先按测到的检测器数，再按检测器编号顺序上第一个有值的分数。

    原实现按 (flagged_detector_count, average_applicable_score, evidence_priority)
    排序——三个量全部依赖已删除的 flag 与跨检测器均值。现在的键不引入任何阈值，
    也不把不同单位的分数混起来比较：只用"测到几层"，同数时用编号最靠前那个检测器的
    分数，保证排序稳定且可解释。
    """
    summary = entry.get("matrix_summary", {})
    matrix = entry.get("evaluation_matrix", {}) or {}
    first = None
    for k in DETECTOR_ORDER:
        v = (matrix.get(k) or {}).get("score")
        if isinstance(v, (int, float)):
            first = v
            break
    return (summary.get("measured_detector_count", 0),
            -1.0 if first is None else first)
