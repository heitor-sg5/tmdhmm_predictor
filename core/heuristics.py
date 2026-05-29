from data.residues import POSITIVE_RESIDUES, AROMATIC_RESIDUES, KD_SCORES

# Positive-inside rule
def check_positive_inside(sequence, tm_start, tm_end, flank=30):
    """
    Apply the positive-inside rule.
    """
    n_flank = sequence[max(0, tm_start - flank): tm_start]
    c_flank = sequence[tm_end: tm_end + flank]

    n_pos = sum(1 for aa in n_flank if aa in POSITIVE_RESIDUES)
    c_pos = sum(1 for aa in c_flank if aa in POSITIVE_RESIDUES)
    delta = n_pos - c_pos

    if delta > 0:
        cytosolic_side = "N-terminus"
        orientation = "N-terminal cytosolic (type II)"
        rule_supported = True
    elif delta < 0:
        cytosolic_side = "C-terminus"
        orientation = "C-terminal cytosolic (type I)"
        rule_supported = True
    else:
        cytosolic_side = "ambiguous"
        orientation = "Cannot determine orientation"
        rule_supported = False

    return {
        "n_positive": n_pos,
        "c_positive": c_pos,
        "delta": delta,
        "cytosolic_side": cytosolic_side,
        "orientation": orientation,
        "rule_supported": rule_supported,
        "n_flank_seq": n_flank,
        "c_flank_seq": c_flank,
    }

# Aromatic anchors
def find_aromatic_anchors(sequence, tm_start, tm_end, flank=6):
    """
    Detect aromatic residues at membrane-water interface (tryptophan belt/cap).
    """
    n_zone = sequence[max(0, tm_start - flank): tm_start + 2]
    c_zone = sequence[max(0, tm_end - 2): tm_end + flank]

    n_aromatics = [(i + max(0, tm_start - flank), aa) for i, aa in enumerate(n_zone) if aa in AROMATIC_RESIDUES]
    c_aromatics = [(i + max(0, tm_end - 2), aa) for i, aa in enumerate(c_zone) if aa in AROMATIC_RESIDUES]

    all_aromatics = n_aromatics + c_aromatics
    trp_count = sum(1 for _, aa in all_aromatics if aa == "W")

    return {
        "n_aromatics": n_aromatics,
        "c_aromatics": c_aromatics,
        "total_aromatics": len(all_aromatics),
        "trp_count": trp_count,
        "belt_present": len(all_aromatics) >= 1,
        "strong_belt": trp_count >= 1,
    }

# Residue composition
def analyse_tm_composition(sequence, tm_start, tm_end):
    """
    Analyse the amino acid composition of the predicted TM region.
    """
    tm_seq = sequence[tm_start:tm_end]
    if not tm_seq:
        return {}

    n = len(tm_seq)
    hydrophobic = set("ILVFAMCWG")
    charged = set("RKDEH")

    hyd_count = sum(1 for aa in tm_seq if aa in hydrophobic)
    chg_count = sum(1 for aa in tm_seq if aa in charged)
    pro_count = tm_seq.count("P")
    mean_kd = sum(KD_SCORES.get(aa, 0.0) for aa in tm_seq) / n

    # Flag potential issues
    warnings = []
    if hyd_count / n < 0.5:
        warnings.append("Low hydrophobic content")
    if chg_count > 0:
        warnings.append(f"{chg_count} charged residue(s) in TM core")
    if pro_count > 1:
        warnings.append(f"{pro_count} proline(s) (helix breaks)")

    return {
        "tm_sequence": tm_seq,
        "length": n,
        "hydrophobic_frac": round(hyd_count / n, 2),
        "charged_count": chg_count,
        "proline_count": pro_count,
        "mean_kd": round(mean_kd, 2),
        "composition_ok": hyd_count / n >= 0.5 and chg_count == 0,
        "warnings": warnings,
    }