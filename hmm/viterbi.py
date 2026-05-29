from hmm.model import (
    STATES,
    load_params,
    build_log_params,
    residue_class,
)

def viterbi(sequence, min_loop_gap=2):
    """
    Run Viterbi decoding on a protein sequence.
    """
    n = len(sequence)
    if n == 0:
        return {"path": [], "log_prob": 0.0, "tm_spans": [], "state_probs": []}

    # Load parameters (trained or hard-coded defaults)
    initial, transitions, emissions = load_params()
    log_init, log_trans, log_emit = build_log_params(initial, transitions, emissions)

    # Initialisation
    # vit[t][s] = best log-prob of any path ending in state s at position t
    vit = [{} for _ in range(n)]
    back = [{} for _ in range(n)]

    rc_0 = residue_class(sequence[0])
    for s in STATES:
        vit[0][s] = log_init[s] + log_emit[s][rc_0]
        back[0][s] = None

    # Recursion
    for t in range(1, n):
        rc_t = residue_class(sequence[t])
        for s in STATES:
            # Find predecessor state with highest score
            best_prev  = None
            best_score = float("-inf")
            for prev in STATES:
                score = vit[t-1][prev] + log_trans[prev][s]
                if score > best_score:
                    best_score = score
                    best_prev = prev
            vit[t][s] = best_score + log_emit[s][rc_t]
            back[t][s] = best_prev

    # Termination
    best_final = max(STATES, key=lambda s: vit[n-1][s])
    log_prob = vit[n-1][best_final]

    # Traceback
    path = [best_final]
    for t in range(n-1, 0, -1):
        path.append(back[t][path[-1]])
    path.reverse()

    # Extract TM spans and enforce minimum loop gap between domains
    tm_spans = _extract_spans(path, "M")
    tm_spans = _enforce_min_gap(tm_spans, min_gap=min_loop_gap)

    return {
        "path": path,
        "log_prob": round(log_prob, 4),
        "tm_spans": tm_spans,
        "tm_domain_count": len(tm_spans),
        "state_probs": path,
    }

def _extract_spans(path, target):
    """Extract contiguous runs of target state from a state path."""
    spans = []
    in_run = False
    start = 0
    for i, s in enumerate(path):
        if s == target and not in_run:
            in_run = True
            start = i
        elif s != target and in_run:
            spans.append({"start": start, "end": i, "length": i - start})
            in_run = False
    if in_run:
        spans.append({"start": start, "end": len(path), "length": len(path) - start})
    return spans

def _enforce_min_gap(spans, min_gap=2):
    """
    Merge neighbouring TM spans when the intervening loop is shorter than min_gap.
    """
    if not spans:
        return []

    merged = [dict(spans[0])]
    for span in spans[1:]:
        prev = merged[-1]
        gap = span["start"] - prev["end"]
        if gap < min_gap:
            prev["end"] = span["end"]
            prev["length"] = prev["end"] - prev["start"]
        else:
            merged.append(dict(span))
    return merged

def candidate_guided_tm_domains(sequence, candidates, flank=10, min_loop_gap=2):
    """
    Refine each KD candidate with local HMM decoding, then merge globally.
    """
    if not sequence or not candidates:
        return []

    refined = []
    for cand in sorted(candidates, key=lambda c: c["start"]):
        c_start = cand["start"]
        c_end = cand["end"]

        w_start = max(0, c_start - flank)
        w_end = min(len(sequence), c_end + flank)
        local_seq = sequence[w_start:w_end]
        local = viterbi(local_seq, min_loop_gap=min_loop_gap)

        best_overlap = -1
        best_span = None
        for span in local["tm_spans"]:
            g_start = w_start + span["start"]
            g_end = w_start + span["end"]
            overlap = max(0, min(g_end, c_end) - max(g_start, c_start))
            if overlap > best_overlap:
                best_overlap = overlap
                best_span = {"start": g_start, "end": g_end, "length": g_end - g_start}

        # If local HMM misses candidate, keep the KD candidate as fallback.
        if best_span is None or best_overlap <= 0:
            best_span = {"start": c_start, "end": c_end, "length": c_end - c_start}

        refined.append(best_span)

    refined = _enforce_min_gap(refined, min_gap=min_loop_gap)
    return refined

def build_alternating_path(seq_len, tm_spans, start_loop_state="C"):
    """
    Build full-sequence C/M/E path from TM domains with alternating loops.
    """
    if seq_len <= 0:
        return []

    path = [start_loop_state] * seq_len
    loop_state = start_loop_state

    for span in sorted(tm_spans, key=lambda s: s["start"]):
        s = max(0, span["start"])
        e = min(seq_len, span["end"])
        if e <= s:
            continue

        for i in range(s, e):
            path[i] = "M"

        # Flip aqueous side after each membrane crossing.
        loop_state = "E" if loop_state == "C" else "C"
        for i in range(e, seq_len):
            if path[i] != "M":
                path[i] = loop_state

    return path