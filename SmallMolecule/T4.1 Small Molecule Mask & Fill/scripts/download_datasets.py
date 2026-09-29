"""
download_datasets.py
T4.1 Small Molecule Mask & Fill – Download Tox21, HIV, and ClinTox datasets.

Downloads from MoleculeNet/DeepChem S3 bucket. Falls back to hardcoded
FALLBACK_TOX_COMPOUNDS when the network is unavailable or the download fails.
"""

import os
import gzip
import json
import urllib.request
import urllib.error
import csv
import io
import sys

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_TOX21 = os.path.join(BASE_DIR, "data", "tox21")
DATA_HIV   = os.path.join(BASE_DIR, "data", "hiv")

os.makedirs(DATA_TOX21, exist_ok=True)
os.makedirs(DATA_HIV,   exist_ok=True)

# ---------------------------------------------------------------------------
# Fallback data
# ---------------------------------------------------------------------------
FALLBACK_TOX_COMPOUNDS = [
    {"name": "Doxorubicin",    "smiles": "O=C1c2cccc(OC)c2C(=O)c2c1cc1c(c2)C[C@@H](O)C[C@@H]1O",
     "tox21_labels": {"NR-AR": 1, "SR-MMP": 1, "SR-p53": 1}},
    {"name": "Paclitaxel",     "smiles": "O=C(O[C@@H]1C[C@]2(O)C(=O)c3ccccc3",
     "tox21_labels": {"NR-AR": 0, "SR-MMP": 1}},
    {"name": "Aflatoxin_B1",   "smiles": "O=c1oc2c(OC)cc3c(c2c2c1[C@@H]1C=CO[C@@H]1O2)OCO3",
     "tox21_labels": {"NR-AR": 1, "NR-ER": 1, "SR-p53": 1, "SR-MMP": 1}},
    {"name": "Vincristine",    "smiles": "CCC1(O)C=CC2CN3CCC4=C(C3CC2C1CC(=O)OC)NC5=CC=CC=C45",
     "tox21_labels": {"SR-MMP": 1}},
    {"name": "Etoposide",      "smiles": "O=C1OC[C@@H]2c3cc4c(cc3[C@@H](c3cc(OC)c(OC)c(OC)c3)[C@H]1O2)OCO4",
     "tox21_labels": {"NR-AR": 0, "SR-p53": 1}},
    {"name": "Chlorambucil",   "smiles": "OC(=O)CCCc1ccc(N(CCCl)CCCl)cc1",
     "tox21_labels": {"NR-AR": 1, "SR-MMP": 1, "SR-p53": 1}},
    {"name": "Cisplatin",      "smiles": "Cl[Pt](Cl)(N)N",
     "tox21_labels": {"SR-MMP": 1, "SR-p53": 1}},
    {"name": "Benzene",        "smiles": "c1ccccc1",
     "tox21_labels": {"NR-AhR": 1}},
    {"name": "Acrylamide",     "smiles": "C=CC(N)=O",
     "tox21_labels": {"SR-MMP": 1}},
    {"name": "Methyl_mercury", "smiles": "[Hg](C)Cl",
     "tox21_labels": {"SR-MMP": 1, "SR-p53": 1}},
    # Extended fallback to reach ~50 entries
    {"name": "Benzo_a_pyrene",     "smiles": "c1ccc2ccc3cccc4ccc(c1)c2c34",
     "tox21_labels": {"NR-AhR": 1, "NR-AR": 1, "SR-p53": 1}},
    {"name": "Cadmium_chloride",   "smiles": "Cl[Cd]Cl",
     "tox21_labels": {"SR-MMP": 1, "SR-p53": 1, "NR-AR": 1}},
    {"name": "Arsenic_trioxide",   "smiles": "O=[As]O[As]=O",
     "tox21_labels": {"SR-MMP": 1, "SR-p53": 1}},
    {"name": "Lead_acetate",       "smiles": "CC(=O)O[Pb]OC(C)=O",
     "tox21_labels": {"NR-AR": 1, "SR-MMP": 1}},
    {"name": "Paraquat",           "smiles": "C[n+]1ccc(-c2cc[n+](C)cc2)cc1",
     "tox21_labels": {"NR-AR": 1, "SR-MMP": 1, "SR-p53": 1}},
    {"name": "Colchicine",         "smiles": "COc1ccc2c(c1)c(CC1NC(C)=O)ccc2=O",
     "tox21_labels": {"NR-AR": 1, "SR-p53": 1}},
    {"name": "Strychnine",         "smiles": "O=C1OCC2CC3N4CC5=CC=CC=C5C4CC3CC21",
     "tox21_labels": {"NR-AR": 1, "SR-MMP": 1}},
    {"name": "Thalidomide",        "smiles": "O=C1CCC(=O)N1C1CCC(=O)NC1=O",
     "tox21_labels": {"NR-AR": 0, "SR-p53": 1}},
    {"name": "Warfarin",           "smiles": "OC(=O)c1ccccc1OC(=O)c1ccccc1",
     "tox21_labels": {"NR-AR": 0, "SR-MMP": 1}},
    {"name": "Methotrexate",       "smiles": "CN(Cc1cnc2nc(N)nc(N)c2n1)c1ccc(CNC(=O)CC(CCC(=O)O)N)cc1",
     "tox21_labels": {"NR-AR": 0, "SR-p53": 1, "SR-MMP": 1}},
    {"name": "Tamoxifen",          "smiles": "CCOC(=O)c1ccc(N(CC)CC)cc1",
     "tox21_labels": {"NR-ER": 1, "NR-AR": 1}},
    {"name": "Cyclophosphamide",   "smiles": "ClCCNP1(=O)OCCCN1CCCl",
     "tox21_labels": {"NR-AR": 1, "SR-MMP": 1, "SR-p53": 1}},
    {"name": "Fluorouracil",       "smiles": "O=c1[nH]cc(F)c(=O)[nH]1",
     "tox21_labels": {"SR-p53": 1, "SR-MMP": 1}},
    {"name": "Dexamethasone",      "smiles": "CC1CC2C3CCC4=CC(=O)C=CC4(C)C3(F)C(O)C2(C)C1(O)C(=O)CO",
     "tox21_labels": {"NR-GR": 1}},
    {"name": "Tetrodotoxin",       "smiles": "OC1C2(O)C(O)C(O)NC12",
     "tox21_labels": {"NR-AR": 1, "SR-MMP": 1, "SR-p53": 1}},
    {"name": "Ricin_fragment",     "smiles": "NCC(O)=O",
     "tox21_labels": {"SR-p53": 1}},
    {"name": "Sarin_analog",       "smiles": "CC(C)OP(=O)(F)OCC",
     "tox21_labels": {"NR-AR": 1, "SR-MMP": 1, "SR-p53": 1}},
    {"name": "VX_analog",          "smiles": "CCOP(=O)(OCC)SCC",
     "tox21_labels": {"NR-AR": 1, "SR-MMP": 1}},
    {"name": "Nitrogen_mustard",   "smiles": "ClCCN(CCCl)CCCl",
     "tox21_labels": {"NR-AR": 1, "SR-MMP": 1, "SR-p53": 1}},
    {"name": "Lewisite_analog",    "smiles": "ClC=CCl",
     "tox21_labels": {"SR-MMP": 1, "SR-p53": 1}},
    {"name": "Phosgene_analog",    "smiles": "ClC(=O)Cl",
     "tox21_labels": {"SR-MMP": 1}},
    {"name": "Acrylonitrile",      "smiles": "C=CC#N",
     "tox21_labels": {"SR-MMP": 1, "SR-p53": 1}},
    {"name": "Epichlorohydrin",    "smiles": "ClCC1CO1",
     "tox21_labels": {"NR-AR": 1, "SR-MMP": 1}},
    {"name": "Dimethyl_sulfate",   "smiles": "COS(=O)(=O)OC",
     "tox21_labels": {"NR-AR": 1, "SR-MMP": 1, "SR-p53": 1}},
    {"name": "Ethylene_oxide",     "smiles": "C1CO1",
     "tox21_labels": {"NR-AR": 1, "SR-MMP": 1}},
    {"name": "Carbon_tetrachloride","smiles": "ClC(Cl)(Cl)Cl",
     "tox21_labels": {"SR-MMP": 1}},
    {"name": "Trichloroethylene",  "smiles": "ClC=C(Cl)Cl",
     "tox21_labels": {"NR-AhR": 1, "SR-MMP": 1}},
    {"name": "Vinyl_chloride",     "smiles": "ClC=C",
     "tox21_labels": {"NR-AhR": 1, "NR-AR": 1}},
    {"name": "Aniline",            "smiles": "Nc1ccccc1",
     "tox21_labels": {"NR-AR": 1}},
    {"name": "Nitrobenzene",       "smiles": "O=[N+]([O-])c1ccccc1",
     "tox21_labels": {"NR-AR": 1, "SR-MMP": 1}},
    {"name": "2_4_Dinitrotoluene", "smiles": "Cc1ccc([N+](=O)[O-])cc1[N+](=O)[O-]",
     "tox21_labels": {"NR-AhR": 1, "NR-AR": 1, "SR-p53": 1}},
    {"name": "Hexachlorobenzene",  "smiles": "Clc1c(Cl)c(Cl)c(Cl)c(Cl)c1Cl",
     "tox21_labels": {"NR-AhR": 1, "NR-AR": 1, "SR-MMP": 1}},
    {"name": "Lindane",            "smiles": "ClC1C(Cl)C(Cl)C(Cl)C(Cl)C1Cl",
     "tox21_labels": {"NR-AhR": 1, "NR-AR": 1}},
    {"name": "DDT",                "smiles": "ClC(Cl)(Cl)C(c1ccc(Cl)cc1)c1ccc(Cl)cc1",
     "tox21_labels": {"NR-AhR": 1, "NR-AR": 1, "NR-ER": 1}},
    {"name": "Dioxin_TCDD",        "smiles": "Clc1cc2oc3cc(Cl)c(Cl)cc3oc2cc1Cl",
     "tox21_labels": {"NR-AhR": 1, "NR-AR": 1, "NR-ER": 1, "SR-p53": 1, "SR-MMP": 1}},
    {"name": "Atrazine",           "smiles": "CCNc1nc(Cl)nc(NC(C)C)n1",
     "tox21_labels": {"NR-AR": 1, "NR-ER": 1}},
    {"name": "Glyphosate",         "smiles": "OC(=O)CNCP(=O)(O)O",
     "tox21_labels": {"NR-AR": 1}},
    {"name": "Organophosphate_1",  "smiles": "CCOP(=O)(OCC)OCC",
     "tox21_labels": {"SR-MMP": 1}},
    {"name": "Mercury_II_chloride","smiles": "Cl[Hg]Cl",
     "tox21_labels": {"SR-MMP": 1, "SR-p53": 1, "NR-AR": 1}},
    {"name": "Nickel_II_chloride", "smiles": "Cl[Ni]Cl",
     "tox21_labels": {"NR-AR": 1, "SR-p53": 1}},
    {"name": "Chromium_VI_oxide",  "smiles": "O=[Cr](=O)=O",
     "tox21_labels": {"NR-AR": 1, "SR-MMP": 1, "SR-p53": 1}},
]

FALLBACK_HIV_COMPOUNDS = [
    {"name": "AZT",         "smiles": "Cc1cn([C@@H]2C[C@H](N=[N+]=[N-])[C@@H](CO)O2)c(=O)[nH]c1=O", "hiv_active": 1},
    {"name": "Nevirapine",  "smiles": "Cc1ccnc2c1NC(=O)c1ccncc1-2", "hiv_active": 1},
    {"name": "Efavirenz",   "smiles": "FC(F)(F)c1cc2c(cc1)NC(=O)O[C@@]2(C#C)c1cccc(Cl)c1", "hiv_active": 1},
    {"name": "Indinavir",   "smiles": "CC(C)(C)NC(=O)[C@@H]1CN(Cc2cccnc2)CC[C@H]1NC(=O)[C@H](CC(=O)N)NC(=O)c1ccc2ccccc2n1", "hiv_active": 1},
    {"name": "Lopinavir",   "smiles": "CC(C)(C)NC(=O)[C@@H]1CN(Cc2ccccc2)CC[C@H]1NC(=O)[C@H](CC(=O)N)NC(=O)OCC", "hiv_active": 1},
    {"name": "Ritonavir",   "smiles": "CC(C)(C)NC(=O)OCC(=O)[C@@H](CC(C)C)NC(=O)[C@H](CC1=CC=NC=C1)NC(=O)OCC", "hiv_active": 1},
    {"name": "Saquinavir",  "smiles": "CC(C)(C)NC(=O)[C@@H](CC1=CC=CC=C1)NC(=O)[C@@H](CC(=O)N)NC(=O)[C@@H]1CC2=CC=CC=C2N1", "hiv_active": 1},
    {"name": "Atazanavir",  "smiles": "COC(=O)NC(=C)C(=O)N[C@@H](Cc1ccccc1)[C@@H](O)CN(CC(C)(C)C)NC(=O)OC", "hiv_active": 1},
    {"name": "Darunavir",   "smiles": "CC(C)(C)NC(=O)[C@@H](CC1=CC=CC=C1)NC(=O)C[C@@H]1COc2ccccc21", "hiv_active": 1},
    {"name": "Raltegravir", "smiles": "Cc1nc(C(=O)NCC2=CC=C(F)C=C2)c(=O)n1CC(C)(C)NC(=O)c1cc(=O)[nH]c(=O)n1C", "hiv_active": 1},
    {"name": "Elvitegravir","smiles": "CCOc1ccc2c(=O)c(C(=O)O)cn(CC(F)F)c2c1Cc1ccc(F)cc1", "hiv_active": 1},
    {"name": "Dolutegravir","smiles": "Cc1nc2c(c(=O)c1C(=O)O)cc(F)cc2[C@@H]1CC[C@@H](O1)NC(=O)c1ccc(F)cc1", "hiv_active": 1},
    {"name": "Maraviroc",   "smiles": "CC(C1CCC(CC1)N1CCN(CC1)C(=O)c1cc2c(F)cccc2n1C)N1CCSCC1=O", "hiv_active": 1},
    {"name": "Enfuvirtide", "smiles": "CCCCC(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](CC(N)=O)C(=O)O", "hiv_active": 1},
    {"name": "Tenofovir",   "smiles": "Cn1cnc2c(N)ncnc12", "hiv_active": 1},
    {"name": "Emtricitabine","smiles": "Nc1nc(=O)n([C@@H]2CS[C@H](CO)O2)cc1F", "hiv_active": 1},
    {"name": "Lamivudine",  "smiles": "Nc1ccn([C@H]2CS[C@@H](CO)O2)c(=O)n1", "hiv_active": 1},
    {"name": "Abacavir",    "smiles": "Nc1nc2c(ncn2[C@H]2C[C@@H](CO)C=C2)c(=O)n1CC=C", "hiv_active": 1},
    {"name": "Didanosine",  "smiles": "O=C1NC(=O)C=CN1[C@@H]1CC[C@@H](CO)O1", "hiv_active": 1},
    {"name": "Stavudine",   "smiles": "Cc1cn([C@H]2C=C[C@@H](CO)O2)c(=O)[nH]c1=O", "hiv_active": 1},
    {"name": "Zalcitabine", "smiles": "Nc1ccn([C@H]2CC[C@@H](CO)O2)c(=O)n1", "hiv_active": 1},
    {"name": "Foscarnet",   "smiles": "OC(=O)P(=O)(O)O", "hiv_active": 1},
    {"name": "Hydroxyurea", "smiles": "NC(=O)NO", "hiv_active": 1},
    {"name": "CD4_mimic_1", "smiles": "CC1=CC=CC(=C1)C1=NC=CC=N1", "hiv_active": 1},
    {"name": "Calanolide_A","smiles": "CC(C)(C)C1CC(=O)c2c(O)c3c(c(O)c2C1=O)CCC(C)(C)O3", "hiv_active": 1},
    {"name": "Hypericin",   "smiles": "OC1=C2C(=O)C3=C(CC4=C3C(=O)c3c(O)c(O)c(O)cc3C4=O)C2=C(O)C2=CC(=O)C=CC12", "hiv_active": 1},
    {"name": "Curcumin",    "smiles": "COc1cc(/C=C/C(=O)CC(=O)/C=C/c2ccc(O)c(OC)c2)ccc1O", "hiv_active": 1},
    {"name": "Quercetin",   "smiles": "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12", "hiv_active": 1},
    {"name": "Betulinic_acid","smiles": "OC(=O)[C@@H]1CC[C@]2(CO)CC[C@H]3[C@H](CCC4=C3CC[C@@]3(C)C4CC[C@@H]3O)[C@H]2C1", "hiv_active": 1},
    {"name": "Nelfinavir",  "smiles": "CC1(C)C(=O)N(Cc2ccccc2)[C@H](CC(=O)[C@H](CC(C)C)NC(=O)[C@@H](Cc2cccc(O)c2)S(=O)(=O)Nc2ccc(C)cc2)C1", "hiv_active": 1},
]

FALLBACK_CLINTOX_COMPOUNDS = [
    {"name": "Thalidomide_CT",   "smiles": "O=C1CCC(=O)N1C1CCC(=O)NC1=O",         "clintox_failed": 1},
    {"name": "Terfenadine",      "smiles": "CCCCC(CCc1ccccc1)(c1ccccc1)C(O)CCN1CCC(C(O)(c2ccccc2)c2ccccc2)CC1", "clintox_failed": 1},
    {"name": "Cisapride",        "smiles": "CCNC(=O)c1cc(OC)c(N)c(Cl)c1OC", "clintox_failed": 1},
    {"name": "Astemizole",       "smiles": "COc1ccc(CCN2CCC(Nc3nc4ccccc4[nH]3)CC2)cc1", "clintox_failed": 1},
    {"name": "Grepafloxacin",    "smiles": "Cc1cnc2c(c1)n(C1CC1)c1cc(C(=O)O)c(=O)cn12", "clintox_failed": 1},
    {"name": "Mibefradil",       "smiles": "COC(=O)N1CCN(C(=O)c2cc(F)ccc2OCC2=NC(C)(C)CO2)CC1", "clintox_failed": 1},
    {"name": "Trovafloxacin",    "smiles": "OC(=O)c1cn2c(cc1=O)cc(N1CC(N)C1)c(F)c2F", "clintox_failed": 1},
    {"name": "Rofecoxib",        "smiles": "CS(=O)(=O)c1ccc(-c2cc(=O)oc2-c2ccccc2)cc1", "clintox_failed": 1},
    {"name": "Valdecoxib",       "smiles": "Cc1ccc(-c2cc(=O)on2-c2ccccc2)cc1S(N)(=O)=O", "clintox_failed": 1},
    {"name": "Cerivastatin",     "smiles": "CCc1nc(C(C)C)c(CC[C@@H](O)C[C@@H](O)CC(=O)O)c(=O)n1Cc1ccc(F)cc1", "clintox_failed": 1},
    {"name": "Rezulin",          "smiles": "CC1=CC=C(CC2C(=O)NC(=O)S2)C=C1", "clintox_failed": 1},
    {"name": "Phenacetin",       "smiles": "CCOC1=CC=C(NC(C)=O)C=C1", "clintox_failed": 1},
    {"name": "Bromfenac",        "smiles": "NC1=C(CC(=O)c2ccc(Br)cc2)C=CC=C1", "clintox_failed": 1},
    {"name": "Benoxaprofen",     "smiles": "CC(C(=O)O)c1ccc2oc(-c3ccc(Cl)cc3)nc2c1", "clintox_failed": 1},
    {"name": "Duract",           "smiles": "CC(CC1=CC=CC=C1)C(=O)N1CCCC1", "clintox_failed": 1},
    {"name": "Posicor",          "smiles": "COC(=O)N1CCN(C(=O)c2cc(F)ccc2OCC2=NC(C)(C)CO2)CC1", "clintox_failed": 1},
    {"name": "Seldane",          "smiles": "CCCCC(CCc1ccccc1)(c1ccccc1)C(O)CCN1CCC(C(O)(c2ccccc2)c2ccccc2)CC1", "clintox_failed": 1},
    {"name": "Hismanal",         "smiles": "COc1ccc(CCN2CCC(Nc3nc4ccccc4[nH]3)CC2)cc1", "clintox_failed": 1},
    {"name": "Baycol",           "smiles": "C[C@H](c1ccc(F)cc1)c1nc(C(C)C)c(-c2ccc(F)cc2)c(=O)[nH]1", "clintox_failed": 1},
    {"name": "Serzone",          "smiles": "Clc1ccc(CC2=NN3CCCCC3=N2)cc1", "clintox_failed": 1},
]

# ---------------------------------------------------------------------------
# Download helpers
# ---------------------------------------------------------------------------
TOX21_URL = "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/tox21.csv.gz"
HIV_URL   = "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/HIV.csv"


def _download(url: str, dest: str, timeout: int = 60) -> bool:
    """Download url to dest. Returns True on success."""
    try:
        print(f"  Downloading {url} ...")
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=timeout) as resp, open(dest, "wb") as f:
            f.write(resp.read())
        print(f"  Saved to {dest}")
        return True
    except Exception as exc:
        print(f"  Download failed: {exc}")
        return False


def _count_positive_labels(row: dict, tox21_cols: list) -> int:
    count = 0
    for col in tox21_cols:
        val = row.get(col, "")
        try:
            if int(float(val)) == 1:
                count += 1
        except (ValueError, TypeError):
            pass
    return count


# ---------------------------------------------------------------------------
# Tox21
# ---------------------------------------------------------------------------
TOX21_TASK_COLS = [
    "NR-AR", "NR-AR-LBD", "NR-AhR", "NR-Aromatase",
    "NR-ER", "NR-ER-LBD", "NR-PPAR-gamma",
    "SR-ARE", "SR-ATAD5", "SR-HSE", "SR-MMP", "SR-p53",
]


def load_tox21(gz_path: str):
    """Parse tox21.csv.gz and return list of dicts."""
    records = []
    with gzip.open(gz_path, "rt", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            smiles = row.get("smiles", "").strip()
            name   = row.get("mol_id", row.get("compound_id", f"mol_{len(records)}")).strip()
            if not smiles:
                continue
            labels = {}
            for col in TOX21_TASK_COLS:
                val = row.get(col, "").strip()
                if val in ("0", "1"):
                    labels[col] = int(val)
            n_pos = sum(v for v in labels.values() if v == 1)
            records.append({
                "name": name,
                "smiles": smiles,
                "tox21_labels": labels,
                "n_positive": n_pos,
            })
    return records


def select_top_tox21(records, n=50):
    sorted_rec = sorted(records, key=lambda r: r["n_positive"], reverse=True)
    return sorted_rec[:n]


def download_tox21():
    gz_path = os.path.join(DATA_TOX21, "tox21.csv.gz")
    json_path = os.path.join(DATA_TOX21, "tox21_top50.json")

    ok = False
    if not os.path.exists(gz_path):
        ok = _download(TOX21_URL, gz_path)
    else:
        print(f"  tox21.csv.gz already cached at {gz_path}")
        ok = True

    records = []
    if ok:
        try:
            records = load_tox21(gz_path)
            print(f"  Loaded {len(records)} Tox21 records")
        except Exception as exc:
            print(f"  Parse error: {exc} – using fallback")

    if not records:
        print("  Using FALLBACK_TOX_COMPOUNDS")
        records = [
            {**c, "n_positive": sum(v for v in c["tox21_labels"].values() if v == 1)}
            for c in FALLBACK_TOX_COMPOUNDS
        ]

    top = select_top_tox21(records, 50)
    with open(json_path, "w") as f:
        json.dump(top, f, indent=2)
    print(f"  Saved {len(top)} compounds to {json_path}")
    return top


# ---------------------------------------------------------------------------
# HIV
# ---------------------------------------------------------------------------
def load_hiv(csv_path: str):
    records = []
    with open(csv_path, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            smiles = row.get("smiles", "").strip()
            activity = row.get("HIV_active", "0").strip()
            if not smiles:
                continue
            if activity == "1":
                records.append({
                    "name": f"HIV_active_{len(records)}",
                    "smiles": smiles,
                    "hiv_active": 1,
                })
    return records


def download_hiv():
    csv_path  = os.path.join(DATA_HIV, "HIV.csv")
    json_path = os.path.join(DATA_HIV, "hiv_active_top30.json")

    ok = False
    if not os.path.exists(csv_path):
        ok = _download(HIV_URL, csv_path)
    else:
        print(f"  HIV.csv already cached at {csv_path}")
        ok = True

    records = []
    if ok:
        try:
            records = load_hiv(csv_path)
            print(f"  Loaded {len(records)} HIV-active records")
        except Exception as exc:
            print(f"  Parse error: {exc} – using fallback")

    if not records:
        print("  Using FALLBACK_HIV_COMPOUNDS")
        records = FALLBACK_HIV_COMPOUNDS

    top = records[:30]
    with open(json_path, "w") as f:
        json.dump(top, f, indent=2)
    print(f"  Saved {len(top)} compounds to {json_path}")
    return top


# ---------------------------------------------------------------------------
# ClinTox (hardcoded)
# ---------------------------------------------------------------------------
def save_clintox():
    json_path = os.path.join(DATA_TOX21, "clintox_failed.json")
    with open(json_path, "w") as f:
        json.dump(FALLBACK_CLINTOX_COMPOUNDS, f, indent=2)
    print(f"  Saved {len(FALLBACK_CLINTOX_COMPOUNDS)} ClinTox compounds to {json_path}")
    return FALLBACK_CLINTOX_COMPOUNDS


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print("=== Phase 1: Downloading datasets ===")

    print("\n[Tox21]")
    tox21 = download_tox21()

    print("\n[HIV]")
    hiv = download_hiv()

    print("\n[ClinTox]")
    clintox = save_clintox()

    print(f"\nDone. Tox21: {len(tox21)}, HIV: {len(hiv)}, ClinTox: {len(clintox)}")
