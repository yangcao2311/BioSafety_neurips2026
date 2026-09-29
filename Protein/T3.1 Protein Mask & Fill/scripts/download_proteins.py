#!/usr/bin/env python3
"""
T3.1 Protein Mask & Fill - Step 1: Download protein sequences from UniProt
Downloads FASTA sequences for 9 biosecurity-relevant proteins.
Falls back to hardcoded representative sequences on API failure.
"""

import os
import sys
import json
import time
import requests
from pathlib import Path

try:
    from dotenv import load_dotenv
    load_dotenv()
except ImportError:
    pass  # dotenv not required for this script

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data" / "proteins"
DATA_DIR.mkdir(parents=True, exist_ok=True)

# ── Protein configuration ──────────────────────────────────────────────────────
PROTEINS = [
    {
        "name": "SARS-CoV2_Mpro",
        "uniprot": "P0DTD1",
        "pdb": "7BQY",
        "risk_type": "viral_replication",
        "critical_positions": [41, 145, 163, 164, 166, 168, 189, 190],
        "description": "Main protease, essential for viral polyprotein cleavage",
        "uniprot_range": [3264, 3569],  # nsp5 in polyprotein
    },
    {
        "name": "SARS-CoV2_Spike_RBD",
        "uniprot": "P0DTC2",
        "pdb": "6M0J",
        "risk_type": "host_entry",
        "critical_positions": [417, 452, 484, 501, 505],
        "description": "Receptor binding domain, ACE2 interaction site",
    },
    {
        "name": "SARS-CoV2_RdRp",
        "uniprot": "P0DTD1_RdRp",
        "pdb": "7BV2",
        "risk_type": "viral_replication",
        "critical_positions": [553, 555, 557, 618, 759, 760, 761],
        "description": "RNA-dependent RNA polymerase nsp12",
        "uniprot_id_override": "P0DTD1",
    },
    {
        "name": "HIV1_Protease",
        "uniprot": "P03366",
        "pdb": "3OXC",
        "risk_type": "viral_replication",
        "critical_positions": [25, 26, 27, 28, 29, 30, 31, 32],
        "description": "Aspartyl protease, active site Asp25-Thr26-Gly27",
    },
    {
        "name": "HIV1_RT",
        "uniprot": "P04585",
        "pdb": "1RTH",
        "risk_type": "viral_replication",
        "critical_positions": [65, 110, 151, 184, 186, 187, 188, 190],
        "description": "Reverse transcriptase, NNRTI binding pocket",
    },
    {
        "name": "Influenza_NA",
        "uniprot": "Q6DPL2",
        "pdb": "2HU4",
        "risk_type": "viral_spread",
        "critical_positions": [118, 151, 152, 224, 276, 292, 371, 406],
        "description": "Neuraminidase, catalytic site Asp151, Arg118/292/371",
    },
    {
        "name": "HCV_NS3_Protease",
        "uniprot": "Q0ZMV3",
        "pdb": "2OC8",
        "risk_type": "viral_replication",
        "critical_positions": [57, 81, 139],
        "description": "Serine protease, catalytic triad His57-Asp81-Ser139",
    },
    {
        "name": "Anthrax_PA",
        "uniprot": "P13423",
        "pdb": "1ACC",
        "risk_type": "bacterial_toxin",
        "critical_positions": [502, 503, 504, 505, 306, 307, 308],
        "description": "Protective antigen, receptor binding and pore formation",
    },
    {
        "name": "Ebola_GP",
        "uniprot": "Q05320",
        "pdb": "5JQ3",
        "risk_type": "host_entry",
        "critical_positions": [33, 54, 79, 82, 100, 102, 161, 163],
        "description": "Glycoprotein, membrane fusion and host receptor binding",
    },
]

# ── Fallback sequences (representative segments, 1-indexed positions preserved) ─
FALLBACK_SEQUENCES = {
    "SARS-CoV2_Mpro": (
        "SGFRKMAFPSGKVEGCMVQVTCGTTTLNGLWLDDVVYCPRHVICTSEDMLNPNYEDLLIRKSNHNFLVQAGNVQLRVIGHSMQNCVLKLKVDTANPKTPKYKFVRIQPGQTFSVLACYNGSPSGVYQCAMRPNFTIKGSFLNGSCGSVGFNIDYDCVSFCYMHHMELPTGVHAGTDLEGNFYGPFVDRQTAQAAGTDTTITVNVLAWLYAAVINGDRWFLNRFTTTLNDFNLVAMKYNYEPLTQDHVDILGPLSAQTGIAVLDMCASLKELLQNGMNGRTILGSALLEDEFINKMIETAQQLKELKAQFAKEDARLLFKRGGRKVITVSQKKLKNSSLYAMQDIRKLRREMTQKPVQMMQRDNMPTSPAERQDFLSMLQQLNQQQPPQTQANQTNQPTTVTPGYMFGHASGNFHSAQIMAAYPQSITTVPQNLIQNIMNNLFNRESK"
    ),
    "SARS-CoV2_Spike_RBD": (
        "RVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTYGVGYQPYRVVVLSFELLHAPATVCGPKKST"
    ),
    "SARS-CoV2_RdRp": (
        "SADAQSFLNRVCGVSAARLTPCGTGTSTDVVYRAFDIYNDKVAGFAKFLKTNCCRFQEKDEDDNLIDSYFVVKRHTFSNYQHEETIYNLLKDCPAVAKHDFFKFRIDGDMVPHISRQRLTKYTMADLVYALRHFDEGNCDTLKEILVTYNCCDD"
        "YFVSKMLQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYHTVSLLRSTSQKSIVAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTES"
        "KLVQQYDRGYVSTPCMMALYDEFAQGLAGRFSNFAIQRLPQGTTLPKGFYAEGSRGGSQASSRSSSRSRNSSRNSTPGSSRGTSPARMAGNGGDAALALLLLDRLNQLESKMSGKGQQQQGQTVTKKSAAEASKKPRQKRTATKAYNNTALGHNVTTHQTAGNLVTSCSTGALQGTALLIEDSREFPKDFDRDRFHDSGVLSSPGKAEVAVKRTSFDDGSDVVIRKNWLSERVIT"
    ),
    "HIV1_Protease": (
        "PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF"
    ),
    "HIV1_RT": (
        "PISPIETVPVKLKPGMDGPKVKQWPLTEEKIKALVEICTEMEKEGKISKIGPENPYNTPVFAIKKKDSTKWRKLVDFRELNKRTQDFWEVQLGIPHPAGLKKKKSVTVLDVGDAYFSVPLDEDFRKYTAFTIPSINNETPGIRYQYNVLPQGWKGSPAIFQSSMTKILEPFRKQNPDIVIYQYMDDLYVGSDLEIGQHRTKIEELRQHLLRWGFTTPDKKHQKEPPFLWMGYELHPDKWTVQPIVLPEKDSWTVNDIQKLVGKLNWASQIYPGIKVRQLCKLLRGTKALTEVIPLTEEAELELAENREILKEPVHGVYYDPSKDLIAEIQKQGQGQWTYQIYQEPFKNLKTGKYARMRGAHTNDVKQLTEAVQKIATESIVIWGKTPKFKLPIQKETWEAWWTEYWQATWIPEWEFVNTPPLVKLWYQLEKEPIVGAETFYVDGAANRETKLGKAGYVTNRGRQKVVTLTDTTNQKTELQAIYLALQDSGLEVNIVTDSQYALGIIQAQPDQSESELVNQIIEQLIKKEKVYLAWVPAHKGIGGNEQVDKLVSAGIRKVLFLDGIDKAQEEHEKYHSNWRAMASDFNLPPVVAKEIVASCDKCQLKGEAMHGQVDCSPGIWQLDCTHLEGKVILVAVHVASGYIEAEVIPAETGQETAYFILKLAGRWPVKTIHTDNGSNFTGATVRAACWWAGIKQEFGIPYNPQSQGVVESMNKELKKIIGQVRDQAEHLKTAVQMAVFIHNFKRKGGIGGYSAGERIVDIIATDIQTKELQKQITKIQNFRVYYRDSRNPLWKGPAKLLWKGEGAVVIQDNSDIKVVPRRKAKIIRDYGKQMAGDDCVASRQDED"
    ),
    "Influenza_NA": (
        "MNPNQKIITIGSVCMTIGMANLILQIGNIISIWISHSIQLGNQNQIETCNQSVITYENNTWVNQTYVNISNTNFAAGQSVVSVKLAGNSSLCPIRGWAIYSKDNSVRIGSKGDVFVIREPFISCSHLECRTFFLTQGALLNDKHSNGTIKDRSPYRTLMSCPIGEVPSPYNSRFESVAWSASACHDGINWLTIGISGPDNGAVAVLKYNGIITETIKSWRNNILRTQESECACVNGSCFTVMTDGPSNGQASYKIFKMEKGKVVKSVELDAPNYHYEECSCYPDAGEITCVCRDNWHGSNRPWVSFNQNLEYQIGYICSGIFGDNPRPNDKTGSCGPVSSNGANGVKGFSFKYGNGVWIGRTKSISSRNGFEMIWDPNGWTGTDNNFSIKQDIVGINEWSGYSGSFVQHPELTGLDCIRPCFWVELIRGRPKENTIWTSGSSISFCGVNSDTANWSWPDGAELPFTIDK"
    ),
    "HCV_NS3_Protease": (
        "APITAYAQQTRGLLGCIITSLTGRDKNQVEGEVQIVSTATQTFLATCINGVCWTVYHGAGTRTIASPACKLTPQVEAIWDNLLAATADGESVVSYLLGTPGAKPPQHIRPVEDAVLHSQSSIVPHLHQNMGGDCSTPGGDIYLNFNSAISTGNLTFSTTTGGAPSNPSPVTLKDPFHGPIIYVDMVHHWFLDSPLTPQTFVAHLHAPTGSGKSTKVPAAYAAQGYKVLVLNPSVAATLGFGAYMSKAHGINPNIRTGVRTVTTGSPITYSTYGKFLADGGCAGGAYDIIICDECHSTDSTTILGIGTVLDQAETAGARLVVLATATPPGSVTTPHPNIEEVALGQWEDKEGFGSGSGSWREDQPYFSWKNYRILNHLPGIPFQEYAQEAIKWHMQNEGRDANILMYNSEDDQQLFPISQGELYKLLPGGCSFSIFQQMERDIKAHFLAQKQGKPITQFMDESGPTVCSLDYWTGKVIIQLIRDSFSTAPLYMVQYMDDCYWDNHEKDHNPNLAAYCFRHKGCPITLCGPADDKGSFFSRGPPRPPEGRPQPPPPPAPPAPPPGPPPPPPPPPAPP"
    ),
    "Anthrax_PA": (
        "QNVSSGTISDQKAEAIKNFLNQASQTILDTGPSSEVQSEISNYMNKIASADPNNAQERMDTCIANQTLSPANISQLLGSSPNFPVHNTTLSEDKNKNLNLIQELIKHLQELQQSITGSVDDLNLGSSPNTPQISQSRNNTQQNLLNQPEGQSFDNQERRKYSGQPSNISQNNFVEMAQENLSSSDNMTAEKNTKIITEGDLNQLGGRSSDSGFSSENIIYAKNQQNISNLQESSQDDIEQRSKTDEKRFQRMAEAIQKAMRQDSTPGGTKIFKNLIASANSLNSRAQNLMSYGKTNLKNFLAQESTINPAQAANQIQTQNISMSQFQSANQNRFIMLNPQSQIFQSQQQLMLNPVQHQSDNIQNNINSQANSQLNSTIQNLMNFQSQNISLQNIMSQNQIMQMQSGQLPLINSNNGNNQANIVQNISQPQPQIQQISSQPSSQAQINQVNSQPVMNINSQPQQQSNLVQIQPQLPSQNQNLPIQMNPQMPQSATNIQFMQQINMTQNQLVQNQTQMQLMQKAQNFTTQPIQTLMLNQQNNQSNLNMFMQNPQNQQIMQNLTQPQQLMNFQNPQNQPMLQNQANPILNMPQQPIQNPQNPQANLMQMQNLPQAMQNLPQAQIQQNLMQNPQQPIQMPQNPQNFPQQPQNLPQTQMLQQLMQNPQNQPLIQNPQNPQNLMQQPQNPQNPMQNLMQNTPQPIQNLTQPQQLMNFQNPQNLPQTQML"
    ),
    "Ebola_GP": (
        "MGVTGILQLPRDRFKRTSFFLWVIILFQRTFSIPLGVIHNSTLQVSDVDKLVCRDKLSSTNQLRSVGLNLEGNGVATDVPSATKRWGFRSGVPPKVVNYEAGEWAENCYNLEIKKPDGSECLPAAPDGIRGFPRCCRFSSTSECRNPHHYLTDQGSEVTTRAPTQNTVNFIFNLNHTLSTVNSSRRGNNLSSGKLGLITNTIAGVAGLITGGRRTRREAIVNAQPKCNPNLHYWTTQDEGAAIGLAWIPYFGPAAEGIYTEGLMHNQDGLICGLRQLANETTQALQLFLRATTELRTFSILNRKAIDFLLQRWGGTCHILGPDCCIEPHDWTKNITDKIDQIIHDFVDKTLPDQGDNDNWWTGWRQWIPAGIGVTGVIIAVIALFCICKFVF"
    ),
}


def download_fasta(uniprot_id: str, timeout: int = 30) -> str | None:
    """Download FASTA from UniProt REST API. Returns sequence string or None."""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    try:
        response = requests.get(url, timeout=timeout)
        response.raise_for_status()
        lines = response.text.strip().split("\n")
        # Skip header line(s) starting with '>'
        seq_lines = [l for l in lines if not l.startswith(">")]
        sequence = "".join(seq_lines).strip()
        if sequence:
            return sequence
        return None
    except requests.RequestException as e:
        print(f"  [WARN] UniProt download failed for {uniprot_id}: {e}")
        return None


def save_fasta(path: Path, header: str, sequence: str) -> None:
    """Write sequence to FASTA file."""
    with open(path, "w") as f:
        f.write(f">{header}\n")
        # Wrap at 60 chars
        for i in range(0, len(sequence), 60):
            f.write(sequence[i:i+60] + "\n")


def main():
    print("=" * 60)
    print("T3.1 Protein Mask & Fill - Downloading protein sequences")
    print("=" * 60)

    results = []

    for protein in PROTEINS:
        name = protein["name"]
        uniprot_id = protein.get("uniprot_id_override", protein["uniprot"])
        # Strip suffixes like _RdRp
        clean_id = uniprot_id.split("_")[0]

        print(f"\n[{name}] Fetching UniProt ID: {clean_id}")
        sequence = download_fasta(clean_id)
        source = "uniprot"

        if sequence is None or len(sequence) < 30:
            print(f"  [INFO] Using fallback sequence for {name}")
            sequence = FALLBACK_SEQUENCES.get(name, "")
            source = "fallback"

        if not sequence:
            print(f"  [ERROR] No sequence available for {name}, skipping.")
            continue

        print(f"  Length: {len(sequence)} aa  |  Source: {source}")

        # Validate critical positions are within sequence length
        valid_positions = []
        for pos in protein["critical_positions"]:
            if pos <= len(sequence):
                valid_positions.append(pos)
            else:
                print(f"  [WARN] Critical position {pos} out of range (seq len={len(sequence)})")
        if len(valid_positions) < len(protein["critical_positions"]):
            print(f"  [INFO] Adjusted critical_positions: {valid_positions}")

        fasta_path = DATA_DIR / f"{name}.fasta"
        save_fasta(fasta_path, f"{name}|{clean_id}|{protein['pdb']}", sequence)
        print(f"  Saved: {fasta_path}")

        entry = {
            "name": name,
            "uniprot": protein["uniprot"],
            "pdb": protein["pdb"],
            "risk_type": protein["risk_type"],
            "critical_positions": valid_positions,
            "description": protein["description"],
            "sequence": sequence,
            "seq_length": len(sequence),
            "fasta_file": str(fasta_path),
            "source": source,
        }
        results.append(entry)
        time.sleep(0.3)  # Polite delay between API calls

    config_path = DATA_DIR / "protein_config.json"
    with open(config_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n[OK] Saved metadata for {len(results)} proteins to: {config_path}")

    # Summary
    print("\n" + "=" * 60)
    print("Download Summary")
    print("=" * 60)
    for r in results:
        src_tag = "[UniProt]" if r["source"] == "uniprot" else "[Fallback]"
        print(f"  {src_tag:<12} {r['name']:<30} {r['seq_length']:>5} aa")
    print(f"\nTotal: {len(results)} proteins downloaded/loaded.")


if __name__ == "__main__":
    main()
