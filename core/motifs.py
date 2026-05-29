MOTIF_LIBRARY = [
    "G/AXXXG/A",
    "S/TXXXS/T",
    "GXXXXXXG",
    "GXXXXXG",
    "N/QXXXG",
    "S/TXXXG",
    "PXXXP",
    "PXXP",
    "GXXP",
]

def _compile_pattern(pattern):
    """
    Compile motif string into per-position allowed residue sets.
    """
    pattern = pattern.upper().replace(" ", "")
    if not pattern:
        return []

    positions = []
    attach_to_previous = False

    for token in pattern:
        if token == "/":
            # Next residue should be added to previous position.
            attach_to_previous = True
            continue

        if token == "X":
            # Wildcard position.
            positions.append(None)
            attach_to_previous = False
            continue

        # Standard residue token.
        if attach_to_previous and positions and positions[-1] is not None:
            positions[-1].add(token)
        else:
            positions.append({token})
        attach_to_previous = False

    return positions

def _find_pattern_hits(sequence, pattern, offset=0):
    """Find all motif hits using compiled positions."""
    hits = []
    compiled = _compile_pattern(pattern)
    pattern_len = len(compiled)

    if pattern_len == 0 or len(sequence) < pattern_len:
        return hits

    for i in range(len(sequence) - pattern_len + 1):
        window = sequence[i : i + pattern_len]
        is_match = True
        for aa, allowed in zip(window, compiled):
            if allowed is not None and aa not in allowed:
                is_match = False
                break
        if is_match:
            start = offset + i
            end = start + pattern_len
            hits.append(
                {
                    "start": start,
                    "end": end,
                    "sequence": window,
                }
            )
    return hits

def detect_tm_motifs(sequence, tm_start, tm_end, motif_names=None):
    """
    Detect registered motifs inside the predicted TM segment.
    """
    tm_seq = sequence[tm_start:tm_end]
    active_patterns = motif_names or MOTIF_LIBRARY

    motifs = []
    total_hits = 0

    for pattern in active_patterns:
        if pattern not in MOTIF_LIBRARY:
            continue
        hits = _find_pattern_hits(tm_seq, pattern, offset=tm_start)
        motifs.append(
            {
                "name": pattern,
                "pattern": pattern,
                "hits": hits,
                "count": len(hits),
            }
        )
        total_hits += len(hits)

    motifs_with_hits = [m["name"] for m in motifs if m["count"] > 0]

    return {
        "tm_sequence": tm_seq,
        "tm_start": tm_start,
        "tm_end": tm_end,
        "motifs": motifs,
        "total_hits": total_hits,
        "motifs_with_hits": motifs_with_hits,
        "has_motif_support": total_hits > 0,
    }