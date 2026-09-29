#!/usr/bin/env python3
"""Integrate Channel A + Channel B risks"""
import os, json

def integrate(channel_a_dir, channel_b_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    a_file = os.path.join(channel_a_dir, "channelA_results.json")
    b_file = os.path.join(channel_b_dir, "channelB_results.json")
    with open(a_file) as f: a_results = json.load(f)
    with open(b_file) as f: b_results = json.load(f)
    
    a_map = {r["id"]:r for r in a_results}
    b_map = {r["id"]:r for r in b_results}
    all_ids = set(list(a_map.keys()) + list(b_map.keys()))
    
    integrated = []
    for sid in all_ids:
        a = a_map.get(sid, {})
        b = b_map.get(sid, {})
        a_score = a.get("offtarget_count", 0)
        b_risk = b.get("expression_weighted_risk", {})
        b_score = b_risk.get("total_risk", 0) if isinstance(b_risk, dict) else 0
        
        alpha, beta = 0.4, 0.6
        composite = alpha * a_score + beta * (b_score / 100)
        
        risk_a = a.get("risk_a", "Low")
        risk_b = b_risk.get("risk_level", "Low") if isinstance(b_risk, dict) else "Low"
        levels = ["Low","Medium","High","Critical"]
        overall = max([risk_a, risk_b], key=lambda x: levels.index(x) if x in levels else 0)
        
        integrated.append({"id": sid, "gene": a.get("gene", b.get("gene","unknown")),
            "strategy": a.get("strategy", b.get("strategy","unknown")),
            "channel_a_offtargets": a_score, "channel_a_risk": risk_a,
            "channel_b_total_risk": round(b_score,4) if b_score else 0,
            "channel_b_risk": risk_b,
            "composite_score": round(composite, 4),
            "overall_risk": overall})
    
    with open(os.path.join(output_dir, "composite_risk_all.json"), "w") as f:
        json.dump(integrated, f, indent=2, ensure_ascii=False)
    print(f"✅ Integration complete: {len(integrated)} results")

if __name__=="__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--channel-a", default="results/offtarget_channelA")
    p.add_argument("--channel-b", default="results/offtarget_channelB")
    p.add_argument("--output", default="results/integrated")
    a = p.parse_args()
    integrate(a.channel_a, a.channel_b, a.output)
