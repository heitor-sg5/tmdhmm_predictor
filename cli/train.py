import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from hmm.model import DEFAULT_INITIAL, DEFAULT_TRANSITIONS, DEFAULT_EMISSIONS, save_params, load_params
from hmm.trainer import estimate_supervised

def build_state_sequence(seq_len, tm_ranges):
    """
    Convert one or more TM annotations to a per-residue state list.
    Assumes the N-terminus is cytosolic and alternates loop labels
    between C and E across successive membrane spans.
    """
    states = [None] * seq_len
    current_side = "C"
    last_end = 0

    for tm_start, tm_end in sorted(tm_ranges):
        # Annotate loop region before the TM segment.
        for i in range(last_end, tm_start - 1):
            states[i] = current_side

        # Annotate TM segment.
        for i in range(tm_start - 1, min(tm_end, seq_len)):
            states[i] = "M"

        # Flip aqueous side for the next loop region.
        current_side = "E" if current_side == "C" else "C"
        last_end = min(tm_end, seq_len)

    # Annotate remaining tail after the last TM segment.
    for i in range(last_end, seq_len):
        states[i] = current_side

    # Fill any unannotated positions conservatively as C.
    return [state if state is not None else "C" for state in states]

def parse_uniprot_txt(filepath):
    """
    Parse a UniProt flatfile (.txt) to extract sequences and TM annotations.
    """
    sequences   = []
    annotations = []

    current_seq   = []
    current_tm    = []
    in_sequence   = False

    with open(filepath) as f:
        for line in f:
            line = line.rstrip()

            if line.startswith("FT   TRANSMEM"):
                parts = line.split()
                if len(parts) >= 3 and ".." in parts[2]:
                    try:
                        start, end = parts[2].split("..")
                        current_tm.append((int(start), int(end)))
                    except ValueError:
                        pass

            elif line.startswith("SQ "):
                in_sequence = True
                current_seq = []

            elif line.startswith("//"):
                # End of entry
                in_sequence = False
                if current_seq and current_tm:
                    seq = "".join(current_seq).upper()
                    if tm_ranges := sorted(current_tm):
                        if tm_ranges[-1][1] <= len(seq) and tm_ranges[0][0] >= 1:
                            ann = build_state_sequence(len(seq), tm_ranges)
                            sequences.append(seq)
                            annotations.append(ann)
                current_seq = []
                current_tm  = []

            elif in_sequence and line.startswith(" "):
                # Sequence lines: remove spaces and digits
                seq_part = line.replace(" ", "").replace("\t", "")
                seq_part = "".join(c for c in seq_part if c.isalpha())
                current_seq.append(seq_part)

    print(f"Parsed {len(sequences)} sequences with TM annotations from {filepath}")
    return sequences, annotations

# Main CLI logic
def main():
    parser = argparse.ArgumentParser(
        description="Train the TM Helix Predictor HMM and save parameters to JSON.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--uniprot", help="UniProt TXT file (supervised training)")
    parser.add_argument("--show", action="store_true", help="Print current parameters and exit")
    args = parser.parse_args()

    if args.show:
        init, trans, emis = load_params()
        print("\n-- Initial probabilities --")
        for s, p in init.items():
            print(f"  {s}: {p:.3f}")
        print("\n-- Transition probabilities --")
        for s in ["C", "M", "E"]:
            row = "  ".join(f"{t}={trans[s][t]:.3f}" for t in ["C", "M", "E"])
            print(f"  {s} -> {row}")
        print("\n-- Emission probabilities --")
        for s in ["C", "M", "E"]:
            row = "  ".join(f"{c}={emis[s][c]:.3f}" for c in ["H", "A", "P", "Q"])
            print(f"  {s}: {row}")
        return

    elif args.uniprot:
        print(f"Loading UniProt TXT: {args.uniprot}")
        seqs, anns = parse_uniprot_txt(args.uniprot)
        if not seqs:
            print("No annotated sequences found. Check file format.")
            sys.exit(1)
        else:
            print(f"Running supervised estimation on {len(seqs)} annotated sequences...")
            initial, transitions, emissions = estimate_supervised(seqs, anns)

    else:
        parser.print_help()
        print("\nNo input specified.")
        sys.exit(0)

    save_params(initial, transitions, emissions)
    print("\nTrained parameters saved to `hmm/trained_params.json`.")
    print("To revert to defaults, delete: hmm/trained_params.json")

if __name__ == "__main__":
    main()