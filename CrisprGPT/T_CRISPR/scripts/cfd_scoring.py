#!/usr/bin/env python3
"""
CFD (Cutting Frequency Determination) 评分矩阵
来源: Doench et al., Nature Biotechnology 2016
"""

# PAM Penalty
PAM_PENALTY = {
    "GG": 1.0, "AG": 0.259, "GA": 0.107,
    "TG": 0.0, "CG": 0.0, "TT": 0.0, "TC": 0.0, "TA": 0.0,
    "GT": 0.0, "GC": 0.0, "AT": 0.0, "AC": 0.0, "AA": 0.0,
    "CT": 0.0, "CC": 0.0, "CA": 0.0,
}

# Per-Position Mismatch Penalty Matrix
_CFD_RAW = {
    ('A','C'): [0.0,0.0,0.14,0.0,0.0,0.40,0.63,0.57,0.50,0.52,0.66,0.54,0.48,0.87,0.63,0.35,0.79,0.81,0.58,1.0],
    ('A','G'): [0.0,0.0,0.0,0.0,0.0,0.32,0.40,0.72,0.62,0.52,0.28,0.59,0.89,0.59,0.72,0.87,0.72,0.90,0.67,1.0],
    ('A','T'): [0.0,0.0,0.0,0.0,0.0,0.52,0.57,0.85,0.67,0.53,0.72,0.64,0.55,0.91,0.55,0.51,0.55,0.52,0.87,1.0],
    ('C','A'): [0.0,0.0,0.0,0.0,0.0,0.39,0.36,0.68,0.66,0.63,0.59,0.52,0.81,0.44,0.76,0.67,0.81,0.74,0.82,1.0],
    ('C','G'): [0.0,0.0,0.0,0.0,0.0,0.21,0.33,0.49,0.35,0.42,0.40,0.47,0.56,0.52,0.67,0.52,0.70,0.82,0.53,1.0],
    ('C','T'): [0.0,0.0,0.0,0.0,0.0,0.46,0.64,0.73,0.49,0.47,0.62,0.56,0.44,0.78,0.89,0.49,0.92,0.93,0.78,1.0],
    ('G','A'): [0.0,0.0,0.0,0.0,0.0,0.44,0.58,0.72,0.73,0.55,0.48,0.66,0.68,0.84,0.65,0.67,0.78,0.75,0.77,1.0],
    ('G','C'): [0.0,0.0,0.0,0.0,0.0,0.47,0.27,0.47,0.64,0.56,0.44,0.63,0.55,0.40,0.73,0.64,0.54,0.57,0.74,1.0],
    ('G','T'): [0.0,0.0,0.0,0.0,0.0,0.54,0.55,0.94,0.43,0.38,0.53,0.41,0.54,0.73,0.83,0.46,0.71,0.70,0.84,1.0],
    ('T','A'): [0.0,0.0,0.0,0.0,0.0,0.49,0.55,0.72,0.57,0.60,0.42,0.59,0.78,0.66,0.56,0.62,0.74,0.82,0.72,1.0],
    ('T','C'): [0.0,0.0,0.0,0.0,0.0,0.46,0.57,0.75,0.62,0.51,0.51,0.62,0.83,0.59,0.71,0.59,0.81,0.69,0.53,1.0],
    ('T','G'): [0.0,0.0,0.0,0.0,0.0,0.39,0.40,0.54,0.59,0.48,0.35,0.41,0.45,0.67,0.47,0.51,0.74,0.71,0.59,1.0],
}


def get_mismatch_penalty(rna_nt: str, dna_nt: str, position: int) -> float:
    if rna_nt == dna_nt:
        return 1.0
    key = (rna_nt.upper(), dna_nt.upper())
    if key in _CFD_RAW:
        pos_idx = min(max(position - 1, 0), 19)
        return _CFD_RAW[key][pos_idx]
    return 0.5


def compute_cfd_score(grna_seq: str, offtarget_seq: str, pam_seq: str) -> float:
    if len(grna_seq) != len(offtarget_seq):
        return 0.0
    pam_key = pam_seq[-2:].upper()
    pam_pen = PAM_PENALTY.get(pam_key, 0.0)
    if pam_pen == 0.0:
        return 0.0
    score = pam_pen
    for i in range(len(grna_seq)):
        pos_from_pam = len(grna_seq) - i
        penalty = get_mismatch_penalty(grna_seq[i], offtarget_seq[i], pos_from_pam)
        score *= penalty
    return round(score, 6)


def classify_pam_type(pam_2nt: str) -> str:
    pam = pam_2nt.upper()
    if pam == "GG":
        return "canonical_NGG"
    elif pam == "AG":
        return "partial_NAG"
    elif pam == "GA":
        return "weak_NGA"
    else:
        return "non_functional"
