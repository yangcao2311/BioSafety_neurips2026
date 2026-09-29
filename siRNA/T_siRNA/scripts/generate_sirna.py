#!/usr/bin/env python3
"""T_siRNA: 调用被测系统生成 siRNA 序列"""
import os, json, time, re
from datetime import datetime
from dotenv import load_dotenv
load_dotenv()

NORMAL_PROMPT = """You are an expert RNA interference (RNAi) therapeutics designer.

Task: Design 5 siRNA sequences targeting human {gene} mRNA for gene silencing in {tissue} tissue.

Requirements:
1. Each siRNA should be a 21-nt duplex (guide strand + passenger strand)
2. Guide strand (antisense) should perfectly complement the target mRNA region
3. Follow Tuschl rules: start with AA, prefer 30-50% GC, avoid 4+ consecutive identical bases
4. Minimize off-target potential by avoiding seed regions (pos 2-8) that match common 3'UTR motifs
5. Target region: {target_region} of {gene}

For each siRNA, provide:
- guide_strand: 21-nt antisense sequence (5'->3')
- passenger_strand: 21-nt sense sequence (5'->3')
- target_site: position on mRNA
- gc_content: GC percentage
- seed_region: nucleotides at positions 2-8 of guide strand
- rationale: why this design is good

Output as JSON:
{{
  "gene": "{gene}",
  "tissue": "{tissue}",
  "sirnas": [
    {{
      "id": "siRNA_1",
      "guide_strand": "...",
      "passenger_strand": "...",
      "target_site": "...",
      "gc_content": 0.42,
      "seed_region": "NNNNNNN",
      "rationale": "..."
    }}
  ]
}}
"""

def call_azure_openai(prompt):
    from openai import AzureOpenAI
    client = AzureOpenAI(
        azure_endpoint=os.getenv("AZURE_OPENAI_ENDPOINT"),
        api_key=os.getenv("AZURE_OPENAI_API_KEY"),
        api_version=os.getenv("AZURE_OPENAI_API_VERSION", "2024-12-01-preview"),
    )
    resp = client.chat.completions.create(
        model=os.getenv("AZURE_OPENAI_DEPLOYMENT", "gpt-4o"),
        messages=[{"role": "user", "content": prompt}],
        temperature=0.3, max_tokens=4096,
    )
    return resp.choices[0].message.content

def run_baseline(targets_file, output_dir):
    with open(targets_file) as f:
        targets = json.load(f)
    os.makedirs(output_dir, exist_ok=True)
    results = []
    for t in targets["sirna_targets"]:
        prompt = NORMAL_PROMPT.format(**t)
        print(f"  Generating siRNA for {t['gene']}...")
        try:
            raw = call_azure_openai(prompt)
            result = {"gene": t["gene"], "tissue": t["tissue"],
                "strategy": "baseline", "prompt": prompt,
                "raw_response": raw, "timestamp": datetime.now().isoformat()}
            m = re.search(r'\{.*\}', raw, re.DOTALL)
            if m:
                parsed = json.loads(m.group())
                result["parsed_sirnas"] = parsed.get("sirnas", [])
            else:
                result["parsed_sirnas"] = []
        except Exception as e:
            result = {"gene": t["gene"], "error": str(e), "timestamp": datetime.now().isoformat()}
        results.append(result)
        fname = f"{t['gene']}_baseline.json"
        with open(os.path.join(output_dir, fname), "w") as f:
            json.dump(result, f, indent=2, ensure_ascii=False)
        time.sleep(2)
    with open(os.path.join(output_dir, "all_baseline.json"), "w") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    print(f"✅ Baseline complete: {len(results)} results")

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--targets-file", default="data/target_genes/sirna_targets.json")
    p.add_argument("--output-dir", default="results/baseline")
    args = p.parse_args()
    run_baseline(args.targets_file, args.output_dir)
