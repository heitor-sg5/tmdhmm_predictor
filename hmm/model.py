import json
import math
import os

from data.residues import RESIDUE_CLASS

# State definitions
STATES = ["C", "M", "E"]
EMIT_STATES = STATES # all states emit
EMIT_ALPHA  = ["H", "A", "P", "Q"] # emission alphabet (residue classes)

# Default probabilities
# Initial state probabilities: proteins usually start in cytosol
DEFAULT_INITIAL = {
    "C": 0.65,
    "M": 0.05,
    "E": 0.30,
}

# Transition probabilities P(next | current)
# M self-loop = 0.95 → expected TM length = 1/(1-0.95) = 20 residues
# C/E self-loops = 0.96 → expected loop length = 25 residues
DEFAULT_TRANSITIONS = {
    "C": {"C": 0.96, "M": 0.04, "E": 0.00},
    "M": {"C": 0.025,"M": 0.95, "E": 0.025},
    "E": {"C": 0.00, "M": 0.04, "E": 0.96},
}

# Emission probabilities P(residue_class | state)
# M heavily favours H; C/E are enriched in Q and P
DEFAULT_EMISSIONS: dict[str, dict[str, float]] = {
    # HAPQ
    "C": {"H": 0.25, "A": 0.05, "P": 0.40, "Q": 0.30},
    "M": {"H": 0.75, "A": 0.10, "P": 0.13, "Q": 0.02},
    "E": {"H": 0.22, "A": 0.06, "P": 0.42, "Q": 0.30},
}

# Load trained parameters
TRAINED_PARAMS_PATH = os.path.join(os.path.dirname(__file__), "trained_params.json")

def load_params():
    """
    If hmm/trained_params.json exists (written by the CLI trainer),
    use those parameters. Otherwise fall back to hand-crafted defaults.
    """
    if os.path.exists(TRAINED_PARAMS_PATH):
        with open(TRAINED_PARAMS_PATH) as f:
            params = json.load(f)
        return (
            params["initial"],
            params["transitions"],
            params["emissions"],
        )
    return DEFAULT_INITIAL, DEFAULT_TRANSITIONS, DEFAULT_EMISSIONS

def save_params(initial, transitions, emissions):
    """Save trained parameters (called by the CLI trainer)."""
    os.makedirs(os.path.dirname(TRAINED_PARAMS_PATH), exist_ok=True)
    with open(TRAINED_PARAMS_PATH, "w") as f:
        json.dump(
            {"initial": initial, "transitions": transitions, "emissions": emissions},
            f, indent=2,
        )
    print(f"Saved trained parameters to {TRAINED_PARAMS_PATH}")

# Log probability helpers
def _log(p):
    """Pseudocount; avoids log(0)."""
    return math.log(p + 1e-12)

def build_log_params(initial, transitions, emissions):
    """Convert all probabilities to log-space (used by Viterbi)."""
    log_init = {s: _log(initial[s]) for s in STATES}
    log_trans = {s: {t: _log(transitions[s][t]) for t in STATES} for s in STATES}
    log_emit = {s: {c: _log(emissions[s][c]) for c in EMIT_ALPHA} for s in STATES}
    return log_init, log_trans, log_emit

def residue_class(aa):
    """Map a single amino acid to its emission class."""
    return RESIDUE_CLASS.get(aa.upper(), "P") # default to P (polar) if unknown