from data.residues import KD_SCORES

TM_KD_THRESHOLD = 1.6
DEFAULT_WINDOW = 19

def score_sequence(sequence):
    """Return per-residue KD scores."""
    return [KD_SCORES.get(aa, 0.0) for aa in sequence]

def sliding_window(sequence, window=DEFAULT_WINDOW):
    """
    Compute mean KD score in a sliding window centred at each residue.
    (Sequence edges windows are truncated).
    """
    scores = score_sequence(sequence)
    half = window // 2
    result = []
    for i in range(len(scores)):
        start = max(0, i - half)
        end = min(len(scores), i + half + 1)
        window_mean = sum(scores[start:end]) / (end - start)
        result.append(round(window_mean, 4))
    return result

def find_tm_candidates(window_scores, threshold=TM_KD_THRESHOLD, min_len: int=17, max_len: int=25):
    """
    Identify contiguous runs above the KD threshold as TM candidates.
    Runs shorter than min_len are discarded (too short to span bilayer).
    Runs longer than max_len are kept but flagged (unusual).
    """
    candidates = []
    in_run = False
    run_start = 0

    for i, score in enumerate(window_scores):
        if score >= threshold and not in_run:
            in_run = True
            run_start = i
        elif score < threshold and in_run:
            run_len = i - run_start
            if run_len >= min_len:
                region = window_scores[run_start:i]
                candidates.append({
                    "start": run_start,
                    "end": i,
                    "length": run_len,
                    "mean_score": round(sum(region) / len(region), 3),
                    "peak_score": round(max(region), 3),
                    "long_flag": run_len > max_len,
                })
            in_run = False

    # Handle run extending to end of sequence
    if in_run:
        run_len = len(window_scores) - run_start
        if run_len >= min_len:
            region = window_scores[run_start:]
            candidates.append({
                "start": run_start,
                "end": len(window_scores),
                "length": run_len,
                "mean_score": round(sum(region) / len(region), 3),
                "peak_score": round(max(region), 3),
                "long_flag": run_len > max_len,
            })

    # Sort by mean score (descending)
    candidates.sort(key=lambda c: c["mean_score"], reverse=True)
    return candidates