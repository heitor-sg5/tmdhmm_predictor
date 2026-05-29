from hmm.model import STATES, EMIT_ALPHA, residue_class

PSEUDOCOUNT = 0.1 # Laplace smoothing

# Supervised estimation (when annotated data available)
def estimate_supervised(sequences, annotations):
    """
    Estimate HMM parameters by counting from labelled data.
    """
    # Initialise count tables with pseudocounts
    init_counts = {s: PSEUDOCOUNT for s in STATES}
    trans_counts = {s: {t: PSEUDOCOUNT for t in STATES} for s in STATES}
    emit_counts = {s: {c: PSEUDOCOUNT for c in EMIT_ALPHA} for s in STATES}

    for seq, ann in zip(sequences, annotations):
        if len(seq) != len(ann):
            continue  # skip mismatched pairs

        # Initial state
        init_counts[ann[0]] += 1

        for t, (aa, state) in enumerate(zip(seq, ann)):
            # Emission
            rc = residue_class(aa)
            emit_counts[state][rc] += 1

            # Transition
            if t < len(ann) - 1:
                next_state = ann[t + 1]
                trans_counts[state][next_state] += 1

    # Normalise
    init_total = sum(init_counts.values())
    initial = {s: init_counts[s] / init_total for s in STATES}

    transitions = {}
    for s in STATES:
        total = sum(trans_counts[s].values())
        transitions[s] = {t: trans_counts[s][t] / total for t in STATES}

    emissions = {}
    for s in STATES:
        total = sum(emit_counts[s].values())
        emissions[s] = {c: emit_counts[s][c] / total for c in EMIT_ALPHA}

    return initial, transitions, emissions