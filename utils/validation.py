from data.residues import VALID_AA

def validate_sequence(raw):
    """
    Clean and validate a protein sequence string.
    """
    if not raw or not raw.strip():
        return None

    # Strip FASTA header lines (">")
    lines = raw.strip().splitlines()
    seq_lines = [ln for ln in lines if not ln.startswith(">")]
    seq = "".join(seq_lines).upper().replace(" ", "").replace("\t", "")

    # Remove digits
    seq = "".join(c for c in seq if not c.isdigit())

    if len(seq) < 10:
        return None
    if len(seq) > 2000:
        return None

    invalid = set(seq) - VALID_AA
    if invalid:
        return None

    return seq

def validation_error_message(raw):
    """Return validation error message."""
    if not raw or not raw.strip():
        return "Please enter a protein sequence."

    lines = raw.strip().splitlines()
    seq_lines = [ln for ln in lines if not ln.startswith(">")]
    seq = "".join(seq_lines).upper().replace(" ", "").replace("\t", "")
    seq = "".join(c for c in seq if not c.isdigit())

    if len(seq) < 10:
        return f"Sequence too short ({len(seq)} residues). Minimum is 10."
    if len(seq) > 2000:
        return f"Sequence too long ({len(seq)} residues). Maximum is 2000."

    invalid = sorted(set(seq) - VALID_AA)
    if invalid:
        return f"Unrecognised amino acid codes: {', '.join(invalid)}. Use standard single-letter codes."

    return "Unknown error."