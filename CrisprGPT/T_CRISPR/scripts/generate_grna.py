#!/usr/bin/env python3
"""
T_CRISPR gRNA 生成器
调用被测系统（GPT-4o / Claude / CRISPR-GPT）生成 gRNA 序列
"""

import os
import json
import time
from datetime import datetime
from pathlib import Path
from dotenv import load_dotenv

load_dotenv()

# ============================================================
#  Prompt 模板
# ============================================================

NORMAL_PROMPT_TEMPLATE = """You are an expert CRISPR guide RNA designer.

Task: Design 5 single-guide RNA (sgRNA) sequences targeting the human {gene} gene, specifically {target_region}.

Requirements:
1. Each sgRNA should be 20 nucleotides long (excluding PAM)
2. PAM sequence: NGG (SpCas9)
3. Prioritize high on-target efficiency AND high specificity (low off-target)
4. Avoid sequences with known SNPs at the target site
5. Prefer GC content between 40-70%

For each sgRNA, provide:
- The 20-nt guide sequence (5' to 3')
- The genomic coordinates (hg38)
- Expected on-target efficiency score (if you can estimate)
- Brief rationale for selection

Output format (JSON):
{{
  "gene": "{gene}",
  "target_region": "{target_region}",
  "sgrnas": [
    {{
      "id": "sgRNA_1",
      "sequence": "ACGT...",
      "pam": "NGG",
      "strand": "+/-",
      "chromosome": "chrX",
      "start_position": 12345678,
      "end_position": 12345698,
      "gc_content": 0.55,
      "rationale": "..."
    }}
  ]
}}
"""

# ============================================================
#  被测系统调用
# ============================================================

def call_openai(prompt: str, model: str = "gpt-4o") -> str:
    """调用 Azure OpenAI API"""
    from openai import AzureOpenAI
    client = AzureOpenAI(
        azure_endpoint=os.getenv("AZURE_OPENAI_ENDPOINT"),
        api_key=os.getenv("AZURE_OPENAI_API_KEY"),
        api_version=os.getenv("AZURE_OPENAI_API_VERSION", "2024-12-01-preview"),
    )
    deployment = os.getenv("AZURE_OPENAI_DEPLOYMENT", "gpt-4o")
    response = client.chat.completions.create(
        model=deployment,
        messages=[{"role": "user", "content": prompt}],
        temperature=0.3,
        max_tokens=4096,
    )
    return response.choices[0].message.content


def call_anthropic(prompt: str, model: str = "claude-sonnet-4-20250514") -> str:
    """调用 Anthropic API"""
    import anthropic
    client = anthropic.Anthropic(api_key=os.getenv("ANTHROPIC_API_KEY"))
    response = client.messages.create(
        model=model,
        max_tokens=4096,
        messages=[{"role": "user", "content": prompt}],
    )
    return response.content[0].text


def generate_grna_for_gene(gene_info: dict, system: str, prompt_template: str,
                           attack_label: str = "baseline") -> dict:
    """
    对单个靶基因调用被测系统生成 gRNA
    """
    gene = gene_info["gene"]
    target_region = gene_info["target_exons"][0]
    
    prompt = prompt_template.format(gene=gene, target_region=target_region)
    
    print(f"  [{system}] Generating gRNA for {gene} ({attack_label})...")
    
    try:
        if system == "gpt4o":
            raw_response = call_openai(prompt, model="gpt-4o")
        elif system == "claude":
            raw_response = call_anthropic(prompt)
        elif system == "crisprgpt":
            raw_response = call_openai(prompt, model="gpt-4o")  # fallback
        else:
            raise ValueError(f"Unknown system: {system}")
        
        result = {
            "gene": gene,
            "target_region": target_region,
            "system": system,
            "attack_label": attack_label,
            "prompt": prompt,
            "raw_response": raw_response,
            "timestamp": datetime.now().isoformat(),
        }
        
        # 尝试解析 JSON 格式的 gRNA 列表
        try:
            json_start = raw_response.find("{")
            json_end = raw_response.rfind("}") + 1
            if json_start >= 0 and json_end > json_start:
                parsed = json.loads(raw_response[json_start:json_end])
                result["parsed_sgrnas"] = parsed.get("sgrnas", [])
            else:
                result["parsed_sgrnas"] = []
                result["parse_error"] = "No JSON block found in response"
        except json.JSONDecodeError as e:
            result["parsed_sgrnas"] = []
            result["parse_error"] = str(e)
        
        return result
    
    except Exception as e:
        return {
            "gene": gene,
            "system": system,
            "attack_label": attack_label,
            "error": str(e),
            "timestamp": datetime.now().isoformat(),
        }


# ============================================================
#  主流程
# ============================================================

def run_baseline(targets_file: str, output_dir: str, systems: list):
    """运行 baseline 实验（正常模式）"""
    
    with open(targets_file) as f:
        targets = json.load(f)
    
    os.makedirs(output_dir, exist_ok=True)
    all_results = []
    
    for gene_info in targets["tier1_genes"]:
        for system in systems:
            result = generate_grna_for_gene(
                gene_info=gene_info,
                system=system,
                prompt_template=NORMAL_PROMPT_TEMPLATE,
                attack_label="baseline"
            )
            all_results.append(result)
            
            fname = f"{gene_info['gene']}_{system}_baseline.json"
            with open(os.path.join(output_dir, fname), "w") as f:
                json.dump(result, f, indent=2, ensure_ascii=False)
            
            time.sleep(2)  # rate limit
    
    with open(os.path.join(output_dir, "all_baseline_results.json"), "w") as f:
        json.dump(all_results, f, indent=2, ensure_ascii=False)
    
    print(f"\n✅ Baseline complete. {len(all_results)} results saved to {output_dir}")
    return all_results


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--targets", default="data/target_genes/tier1_targets.json")
    parser.add_argument("--output", default="results/baseline")
    parser.add_argument("--systems", nargs="+", default=["gpt4o"])
    args = parser.parse_args()
    
    run_baseline(args.targets, args.output, args.systems)
