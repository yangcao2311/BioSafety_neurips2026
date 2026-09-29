#!/usr/bin/env python3
"""T_CRISPR Report Generator"""
import json, os
from datetime import datetime
from collections import defaultdict

def generate_report(offtarget_dir, output_dir):
    results_file = os.path.join(offtarget_dir, "all_offtarget_results.json")
    comparison_file = os.path.join(offtarget_dir, "offtarget_comparison.json")
    with open(results_file) as f: results = json.load(f)
    with open(comparison_file) as f: comparison = json.load(f)
    
    report = []
    report.append("# T_CRISPR 实验报告：CRISPR gRNA Off-Target Safety Evaluation")
    report.append(f"\n> 生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    report.append("")
    
    report.append("## 1. 实验概览")
    report.append(f"\n- 总评估 gRNA 数量: {len(results)}")
    bl_count = sum(1 for r in results if r.get("strategy") in ("baseline","unknown"))
    jb_count = len(results) - bl_count
    report.append(f"- Baseline gRNA: {bl_count}")
    report.append(f"- Jailbreak gRNA: {jb_count}")
    systems = set(r.get("system","unknown") for r in results)
    report.append(f"- 被测系统: {', '.join(systems)}")
    genes = set(r.get("gene","unknown") for r in results)
    report.append(f"- 靶基因: {', '.join(genes)}")
    
    report.append("\n## 2. Baseline vs Jailbreak 对比")
    comp = comparison.get("comparison",{})
    b = comparison.get("baseline",{})
    j = comparison.get("jailbreak",{})
    if comp:
        report.append(f"\n| 指标 | Baseline | Jailbreak | 变化 |")
        report.append(f"|------|----------|-----------|------|")
        report.append(f"| 平均 OT Score | {b.get('avg_offtarget_score','N/A')} | {j.get('avg_offtarget_score','N/A')} | {comp.get('score_inflation_ratio','N/A')}x |")
        report.append(f"| 最大 OT Score | {b.get('max_offtarget_score','N/A')} | {j.get('max_offtarget_score','N/A')} | — |")
    
    report.append("\n## 3. Risk Level 分布")
    for label in ["baseline","jailbreak"]:
        dist = comparison.get(label,{}).get("risk_distribution",{})
        if dist:
            report.append(f"\n### {label.capitalize()}")
            report.append(f"\n| Risk Level | Count |")
            report.append(f"|------------|-------|")
            for level in ["Low","Medium","High","Critical"]:
                report.append(f"| {level} | {dist.get(level,0)} |")
    
    report.append("\n## 4. 各攻击策略效果")
    strat_stats = defaultdict(lambda:{"count":0,"scores":[],"risk_levels":[]})
    for r in results:
        strat = r.get("strategy","unknown")
        score = r.get("offtarget_summary",{}).get("offtarget_score",0)
        risk = r.get("offtarget_summary",{}).get("risk_level","Unknown")
        strat_stats[strat]["count"] += 1
        strat_stats[strat]["scores"].append(score)
        strat_stats[strat]["risk_levels"].append(risk)
    report.append(f"\n| 策略 | 数量 | 平均 OT Score | High/Critical 比例 |")
    report.append(f"|------|------|---------------|-------------------|")
    for strat, stats in sorted(strat_stats.items()):
        avg = sum(stats["scores"])/max(len(stats["scores"]),1)
        hc = sum(1 for r in stats["risk_levels"] if r in ("High","Critical"))
        hc_rate = hc/max(stats["count"],1)*100
        report.append(f"| {strat} | {stats['count']} | {avg:.4f} | {hc_rate:.1f}% |")
    
    report.append("\n## 5. 最高风险 gRNA（Top 10）")
    scored = [r for r in results if isinstance(r.get("offtarget_summary"),dict) and "offtarget_score" in r["offtarget_summary"]]
    scored.sort(key=lambda x:x["offtarget_summary"]["offtarget_score"], reverse=True)
    report.append(f"\n| Rank | Gene | System | Strategy | Sequence | OT Score | Risk |")
    report.append(f"|------|------|--------|----------|----------|----------|------|")
    for i, r in enumerate(scored[:10]):
        ot = r["offtarget_summary"]
        report.append(f"| {i+1} | {r['gene']} | {r['system']} | {r['strategy']} | `{r['sequence'][:10]}...` | {ot['offtarget_score']:.2f} | {ot['risk_level']} |")
    
    report.append("\n## 6. 初步结论")
    report.append("\n> 此部分需人工补充分析。以上数据提供了：")
    report.append("> - Baseline vs Jailbreak 的 off-target score 对比")
    report.append("> - 各攻击策略的有效性排名")
    report.append("> - 最高风险 gRNA 的具体信息")
    
    os.makedirs(output_dir, exist_ok=True)
    report_text = "\n".join(report)
    report_path = os.path.join(output_dir, "T_CRISPR_Report.md")
    with open(report_path, "w") as f: f.write(report_text)
    print(f"\n✅ Report generated: {report_path}")
    return report_path

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--offtarget-dir", default="results/offtarget")
    parser.add_argument("--output", default="results/reports")
    args = parser.parse_args()
    generate_report(args.offtarget_dir, args.output)
