def _classify_score(score):
    """Map numeric score to label and confidence level."""
    if score >= 0.65:
        return "TM Helix Predicted", "High", True
    if score >= 0.40:
        return "Possible TM Helix", "Medium", True
    return "No TM Helix Detected", "Low", False


def _compute_region_score(candidate, positive_inside, aromatics, composition, motifs=None):
    """
    Score one TM region from precomputed heuristics.
    `candidate` may be None — KD-specific factors (1–2) are skipped in that case.
    """
    score = 0.0
    factors = []

    # Hydrophobic stretch (factor 1)
    if candidate:
        kd_delta = min((candidate["mean_score"] - 1.6) / 2.0, 0.45)
        kd_delta = max(kd_delta, 0.0)
        score += kd_delta
        factors.append((
            f"Hydrophobic stretch (mean KD={candidate['mean_score']:.2f}, "
            f"length={candidate['length']})",
            round(kd_delta, 2),
            True,
        ))
    else:
        factors.append(("No KD candidate metadata for this region", 0.0, False))

    # Length adequacy (factor 2)
    if candidate and candidate["length"] >= 17:
        score += 0.10
        factors.append((
            f"TM length adequate ({candidate['length']} residues >= 17)",
            0.10,
            True,
        ))
    elif candidate:
        factors.append((
            f"TM length borderline ({candidate['length']} residues)",
            0.0,
            False,
        ))

    # Composition (factor 3)
    if composition and composition.get("composition_ok"):
        score += 0.10
        factors.append((
            f"Hydrophobic composition OK "
            f"({composition.get('hydrophobic_frac', 0):.0%} hydrophobic, "
            f"0 charged residues)",
            0.10,
            True,
        ))
    elif composition and composition.get("warnings"):
        for w in composition["warnings"]:
            factors.append((f"Composition warning: {w}", 0.0, False))

    # Positive-inside rule (factor 4)
    if positive_inside.get("rule_supported"):
        pos_delta = 0.15
        score += pos_delta
        factors.append((
            f"Positive-inside rule supported "
            f"(N:{positive_inside['n_positive']} vs "
            f"C:{positive_inside['c_positive']} Arg/Lys, "
            f"delta={abs(positive_inside['delta'])})",
            pos_delta,
            True,
        ))
    else:
        factors.append((
            "Positive-inside rule ambiguous (equal Arg/Lys on both sides)",
            0.0,
            False,
        ))

    # Aromatic anchors (factor 5)
    if aromatics.get("belt_present"):
        score += 0.10
        factors.append((
            f"Aromatic anchor belt present "
            f"({aromatics['total_aromatics']} aromatic residues at interface)",
            0.10,
            True,
        ))
    else:
        factors.append(("No aromatic anchors detected at TM boundaries", 0.0, False))

    # Tryptophan bonus (factor 6)
    if aromatics.get("strong_belt"):
        score += 0.05
        factors.append((
            f"Tryptophan anchor present (Trp={aromatics['trp_count']})",
            0.05,
            True,
        ))

    # No charged residues in TM core (factor 7)
    if composition and composition.get("charged_count", 1) == 0:
        score += 0.05
        factors.append(("No charged residues in TM core", 0.05, True))

    # Motif support (factor 8)
    if motifs and motifs.get("has_motif_support"):
        motif_names = ", ".join(motifs.get("motifs_with_hits", []))
        motif_count = motifs.get("total_hits", 0)
        motif_delta = min(0.10, 0.05 + 0.02 * max(0, motif_count - 1))
        score += motif_delta
        factors.append((
            f"TM motif support ({motif_count} hit(s): {motif_names})",
            round(motif_delta, 2),
            True,
        ))
    elif motifs is not None:
        factors.append(("No TM motifs detected in this region", 0.0, False))

    score = min(score, 1.0)
    label, level, has_tm = _classify_score(score)

    return {
        "score": round(score, 3),
        "label": label,
        "confidence": level,
        "has_tm": has_tm,
        "factors": factors,
    }

def compute_region_confidence(sequence, tm_start, tm_end, candidate=None):
    """
    Confidence for a single TM region (one candidate or one display domain).
    """
    from core.heuristics import (
        check_positive_inside,
        find_aromatic_anchors,
        analyse_tm_composition,
    )
    from core.motifs import detect_tm_motifs

    pos_inside = check_positive_inside(sequence, tm_start, tm_end)
    aromatics = find_aromatic_anchors(sequence, tm_start, tm_end)
    composition = analyse_tm_composition(sequence, tm_start, tm_end)
    motifs = detect_tm_motifs(sequence, tm_start, tm_end)

    return _compute_region_score(
        candidate, pos_inside, aromatics, composition, motifs=motifs
    )

def match_candidate_for_span(candidates, span):
    """Return the KD candidate with largest overlap to a TM span, if any."""
    best = None
    best_overlap = 0
    s0, e0 = span["start"], span["end"]
    for cand in candidates:
        overlap = max(0, min(e0, cand["end"]) - max(s0, cand["start"]))
        if overlap > best_overlap:
            best_overlap = overlap
            best = cand
    return best if best_overlap > 0 else None

def compute_global_confidence(sequence, candidates):
    """
    Global confidence = mean score across all KD candidates.
    """
    if not candidates:
        label, level, has_tm = _classify_score(0.0)
        return {
            "score": 0.0,
            "label": label,
            "confidence": level,
            "has_tm": has_tm,
            "factors": [("No hydrophobic stretch >=17 residues found", 0.0, False)],
            "per_region": [],
            "candidate_count": 0,
            "domain_count": 0,
        }

    per_region = []
    for cand in candidates:
        per_region.append(
            compute_region_confidence(
                sequence, cand["start"], cand["end"], candidate=cand
            )
        )

    avg_score = sum(r["score"] for r in per_region) / len(per_region)
    avg_score = min(round(avg_score, 3), 1.0)
    label, level, has_tm = _classify_score(avg_score)

    best = max(per_region, key=lambda r: r["score"])

    return {
        "score": avg_score,
        "label": label,
        "confidence": level,
        "has_tm": has_tm,
        "factors": best["factors"],
        "per_region": per_region,
        "candidate_count": len(candidates),
        "domain_count": 0,
    }

def compute_domain_confidences(sequence, display_spans, candidates):
    """Per-domain confidence for each HMM-refined TM span."""
    domain_scores = []
    for span in display_spans:
        cand = match_candidate_for_span(candidates, span)
        domain_scores.append(
            compute_region_confidence(
                sequence, span["start"], span["end"], candidate=cand
            )
        )
    return domain_scores