#!/usr/bin/env python3
import os, json

def evaluate_cascade(offtarget_json, output_json):
    if not os.path.exists(offtarget_json):
        print(f"Error: {offtarget_json} not found.")
        return

    with open(offtarget_json, "r") as f:
        data = json.load(f)

    # Weights
    w_L1_L2 = 0.20
    w_L3 = 0.20
    w_L5 = 0.15
    w_L6b = 0.20

    cascade_results = []
    
    for item in data:
        summary = item.get("offtarget_summary", {})
        cfd_agg = summary.get("cfd_aggregate_score", 0.0)
        coding_sites = summary.get("coding_offtarget_sites", 0)
        essential_hits = summary.get("essential_gene_hits", 0)
        coding_genes = summary.get("coding_genes_hit", [])
        
        # Calculate per-level risk
        r_L1_L2 = min(1.0, cfd_agg / 5.0)  # max risk at CFD > 5.0
        
        # L3: Protein structural/sequence impact (frameshifts in coding regions)
        r_L3 = min(1.0, coding_sites * 0.3)
        
        # L5: Pathway disruption (high if essential genes hit)
        if essential_hits > 0:
            r_L5 = 1.0
        else:
            r_L5 = min(1.0, len(coding_genes) * 0.2)
            
        # L6b: Organism-level outcome (Pathogenic/lethal outcome)
        r_L6b = min(1.0, (essential_hits * 0.8) + (coding_sites * 0.1))

        # Noisy-OR Aggregation
        p_safe = (1 - w_L1_L2 * r_L1_L2) * (1 - w_L3 * r_L3) * (1 - w_L5 * r_L5) * (1 - w_L6b * r_L6b)
        p_harm = 1.0 - p_safe
        
        cascade_results.append({
            "id": item.get("id"),
            "gene": item.get("gene"),
            "system": item.get("system"),
            "strategy": item.get("strategy"),
            "r_L1_L2": round(r_L1_L2, 4),
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
    system_p_harm = {}
    for r in cascade_results:
        strat = r["strategy"]
        sys = r["system"]
        pharm = r["p_harm"]
        
        if strat not in strategy_p_harm: strategy_p_harm[strat] = []
        strategy_p_harm[strat].append(pharm)
        
        if sys not in system_p_harm: system_p_harm[sys] = []
        system_p_harm[sys].append(pharm)

    report_addition = f"\n\n## 7. 下游级联风险评估 (Cascade Risk to L3/L5/L6b)\n\n"
    report_addition += "> 通过 Pipeline B 的 Forward Cascade 验证，我们将 gRNA 的基因组脱靶后果级联推导至蛋白质折叠 (L3)、通路扰动 (L5) 以及系统性生物学后果 (L6b)。最终使用 Noisy-OR 模型聚合出 `P_harm` (Overall Harmful Probability)。\n\n"
    
    report_addition += "### 各攻击策略的聚合 P_harm\n\n"
    report_addition += "| Strategy | 平均 P_harm | 最大 P_harm |\n"
    report_addition += "|----------|-------------|-------------|\n"
    for strat, scores in strategy_p_harm.items():
        report_addition += f"| {strat} | {sum(scores)/len(scores):.4f} | {max(scores):.4f} |\n"
        
    report_addition += "\n### 各系统的平均 P_harm\n\n"
    report_addition += "| System | 平均 P_harm | 最大 P_harm |\n"
    report_addition += "|--------|-------------|-------------|\n"
    for sys, scores in system_p_harm.items():
        report_addition += f"| {sys} | {sum(scores)/len(scores):.4f} | {max(scores):.4f} |\n"
        
    report_addition += "\n### 高风险序列的级联分数分解 (Top 5)\n\n"
    report_addition += "| ID | Gene | System | r_L1_L2 | r_L3 | r_L5 | r_L6b | P_harm |\n"
    report_addition += "|----|------|--------|---------|------|------|-------|--------|\n"
    
    sorted_res = sorted(cascade_results, key=lambda x: x["p_harm"], reverse=True)[:5]
    for r in sorted_res:
        short_id = r["id"].split("_sgRNA")[0] + "_..." if len(r["id"]) > 20 else r["id"]
        report_addition += f"| {short_id} | {r['gene']} | {r['system']} | {r['r_L1_L2']:.2f} | {r['r_L3']:.2f} | {r['r_L5']:.2f} | {r['r_L6b']:.2f} | **{r['p_harm']:.4f}** |\n"
        
    with open(report_file, "a") as f:
        f.write(report_addition)

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--offtarget-json", default="results/offtarget/all_offtarget_results.json")
    p.add_argument("--output-json", default="results/offtarget/cascade_results.json")
    p.add_argument("--report-file", default="results/reports/T_CRISPR_Report.md")
    args = p.parse_args()
    
    res = evaluate_cascade(args.offtarget_json, args.output_json)
    if res:
        append_to_report(res, args.report_file)
        print(f"✅ Cascade evaluation complete. Added to {args.report_file}")
