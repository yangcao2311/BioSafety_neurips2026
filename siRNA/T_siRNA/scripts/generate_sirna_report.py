#!/usr/bin/env python3
"""T_siRNA Report Generator"""
import json, os
from datetime import datetime
from collections import defaultdict

def generate_report(integrated_dir, output_path):
    with open(os.path.join(integrated_dir,"composite_risk_all.json")) as f:
        results = json.load(f)
    
    report = []
    report.append("# T_siRNA 实验报告：siRNA Off-Target Safety Evaluation")
    report.append(f"\n> 生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    report.append("")
    
    report.append("## 1. 实验概览")
    report.append(f"\n- 总评估 siRNA 数量: {len(results)}")
    bl = [r for r in results if r.get("strategy") in ("baseline","unknown")]
    jb = [r for r in results if r.get("strategy") not in ("baseline","unknown")]
    report.append(f"- Baseline siRNA: {len(bl)}")
    report.append(f"- Jailbreak siRNA: {len(jb)}")
    genes = set(r.get("gene","unknown") for r in results)
    report.append(f"- 靶基因: {', '.join(genes)}")
    
    report.append("\n## 2. 双通道风险分析")
    report.append("\n### Channel A: 全长 Off-Target")
    for label, group in [("Baseline",bl),("Jailbreak",jb)]:
        if not group: continue
        a_scores = [r.get("channel_a_offtargets",0) for r in group]
        report.append(f"\n**{label}**: 平均 OT count = {sum(a_scores)/max(len(a_scores),1):.1f}, "
            f"最大 = {max(a_scores) if a_scores else 0}")
    
    report.append("\n### Channel B: Seed-Mediated Off-Target")
    for label, group in [("Baseline",bl),("Jailbreak",jb)]:
        if not group: continue
        b_scores = [r.get("channel_b_total_risk",0) for r in group]
        report.append(f"\n**{label}**: 平均 risk = {sum(b_scores)/max(len(b_scores),1):.1f}, "
            f"最大 = {max(b_scores) if b_scores else 0:.1f}")
    
    report.append("\n## 3. Risk Level 分布")
    report.append("\n| Category | Low | Medium | High | Critical |")
    report.append("|----------|-----|--------|------|----------|")
    for label, group in [("Baseline",bl),("Jailbreak",jb)]:
        dist = defaultdict(int)
        for r in group: dist[r.get("overall_risk","Low")] += 1
        report.append(f"| {label} | {dist['Low']} | {dist['Medium']} | {dist['High']} | {dist['Critical']} |")
    
    report.append("\n## 4. 各攻击策略效果")
    strat_stats = defaultdict(lambda:{"count":0,"composite":[],"risks":[]})
    for r in results:
        s = r.get("strategy","unknown")
        strat_stats[s]["count"] += 1
        strat_stats[s]["composite"].append(r.get("composite_score",0))
        strat_stats[s]["risks"].append(r.get("overall_risk","Low"))
    report.append("\n| 策略 | 数量 | 平均 Composite | High/Critical % |")
    report.append("|------|------|---------------|-----------------|")
    for s, st in sorted(strat_stats.items()):
        avg = sum(st["composite"])/max(len(st["composite"]),1)
        hc = sum(1 for r in st["risks"] if r in ("High","Critical"))/max(st["count"],1)*100
        report.append(f"| {s} | {st['count']} | {avg:.4f} | {hc:.1f}% |")
    
    report.append("\n## 5. 最高风险 siRNA (Top 10)")
    results.sort(key=lambda x:x.get("composite_score",0), reverse=True)
    report.append("\n| Rank | Gene | Strategy | Composite | Ch-A OT | Ch-B Risk | Overall |")
    report.append("|------|------|----------|-----------|---------|-----------|---------|")
    for i, r in enumerate(results[:10]):
        report.append(f"| {i+1} | {r['gene']} | {r['strategy']} | {r['composite_score']:.4f} | "
            f"{r['channel_a_offtargets']} | {r['channel_b_total_risk']:.1f} | {r['overall_risk']} |")
    
    report.append("\n## 6. 结论")
    report.append("\n> 此报告展示了 siRNA 设计系统在 jailbreak 攻击下的安全性变化。")
    report.append("> 双通道评估（全长 off-target + seed-mediated off-target）提供了全面的风险画像。")
    
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as f: f.write("\n".join(report))
    print(f"✅ Report: {output_path}")

if __name__=="__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--integrated", default="results/integrated")
    p.add_argument("--output", default="results/reports/T_siRNA_Report.md")
    a = p.parse_args()
    generate_report(a.integrated, a.output)
