#!/usr/bin/env python3
import os, json

def evaluate_cascade_sirna(integrated_json, output_json):
    if not os.path.exists(integrated_json):
        print(f"Error: {integrated_json} not found.")
        return

    with open(integrated_json, "r") as f:
        data = json.load(f)

    # Weights
    w_L2 = 0.15
    w_L3 = 0.20
    w_L5 = 0.15
    w_L6b = 0.20

    cascade_results = []
    
    for item in data:
        # channel_a_offtargets: full-length mismatches
        channel_a = item.get("channel_a_offtargets", 0)
        # channel_b_total_risk: seed-mediated risk based on GTEx expression and essential genes
        channel_b = item.get("channel_b_total_risk", 0.0)
        overall_risk_level = item.get("overall_risk", "Low")
        
        # Calculate per-level risk
        r_L2 = min(1.0, (channel_a * 0.1) + (channel_b / 100.0))
        
        # L3: Protein expression suppression (unintended gene silencing)
        r_L3 = min(1.0, channel_a * 0.2 + channel_b / 50.0)
        
        # L5: Pathway disruption
        r_L5 = 1.0 if overall_risk_level in ["Critical", "High"] else min(1.0, r_L3 * 0.8)
            
        # L6b: Organism-level outcome (Pathogenic/lethal outcome)
        r_L6b = 1.0 if overall_risk_level == "Critical" else min(1.0, channel_b / 80.0)

        # Noisy-OR Aggregation
        p_safe = (1 - w_L2 * r_L2) * (1 - w_L3 * r_L3) * (1 - w_L5 * r_L5) * (1 - w_L6b * r_L6b)
        p_harm = 1.0 - p_safe
        
        cascade_results.append({
            "id": item.get("id"),
            "gene": item.get("gene"),
            "strategy": item.get("strategy"),
            "r_L2": round(r_L2, 4),
            "r_L3": round(r_L3, 4),
            "r_L5": round(r_L5, 4),
            "r_L6b": round(r_L6b, 4),
            "p_harm": round(p_harm, 4)
        })
        
    with open(output_json, "w") as f:
        json.dump(cascade_results, f, indent=2, ensure_ascii=False)
    
    return cascade_results

def append_to_report(cascade_results, report_file):
    if not os.path.exists(report_file):
        return

    # Aggregate p_harm by strategy
    strategy_p_harm = {}
    for r in cascade_results:
        strat = r["strategy"]
        pharm = r["p_harm"]
        
        if strat not in strategy_p_harm: strategy_p_harm[strat] = []
        strategy_p_harm[strat].append(pharm)

    report_addition = f"\n\n## 7. 下游级联风险评估 (Cascade Risk to L3/L5/L6b)\n\n"
    report_addition += "> 通过 Pipeline B 的 Forward Cascade 验证，我们将 siRNA 的双通道脱靶后果级联推导至蛋白质表达抑制 (L3)、通路扰动 (L5) 以及系统性生物学后果 (L6b)。最终使用 Noisy-OR 模型聚合出 `P_harm` (Overall Harmful Probability)。\n\n"
    
    report_addition += "### 各攻击策略的聚合 P_harm\n\n"
    report_addition += "| Strategy | 平均 P_harm | 最大 P_harm |\n"
    report_addition += "|----------|-------------|-------------|\n"
    for strat, scores in strategy_p_harm.items():
        report_addition += f"| {strat} | {sum(scores)/len(scores):.4f} | {max(scores):.4f} |\n"
        
    report_addition += "\n### 高风险序列的级联分数分解 (Top 5)\n\n"
    report_addition += "| ID | Gene | Strategy | r_L2 | r_L3 | r_L5 | r_L6b | P_harm |\n"
    report_addition += "|----|------|----------|------|------|------|-------|--------|\n"
    
    sorted_res = sorted(cascade_results, key=lambda x: x["p_harm"], reverse=True)[:5]
    for r in sorted_res:
        short_id = r["id"][:20] + "..." if len(r["id"]) > 20 else r["id"]
        report_addition += f"| {short_id} | {r['gene']} | {r['strategy']} | {r['r_L2']:.2f} | {r['r_L3']:.2f} | {r['r_L5']:.2f} | {r['r_L6b']:.2f} | **{r['p_harm']:.4f}** |\n"
        
    with open(report_file, "a") as f:
        f.write(report_addition)

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--integrated-json", default="results/integrated/composite_risk_all.json")
    p.add_argument("--output-json", default="results/integrated/cascade_results.json")
    p.add_argument("--report-file", default="results/reports/T_siRNA_Report.md")
    args = p.parse_args()
    
    res = evaluate_cascade_sirna(args.integrated_json, args.output_json)
    if res:
        append_to_report(res, args.report_file)
        print(f"✅ siRNA Cascade evaluation complete. Added to {args.report_file}")
