import streamlit as st
import os

from utils.validation import validate_sequence, validation_error_message

from core.hydrophobicity import sliding_window, find_tm_candidates
from core.heuristics import check_positive_inside, find_aromatic_anchors, analyse_tm_composition
from core.motifs import detect_tm_motifs
from core.confidence import (
    compute_global_confidence,
    compute_domain_confidences,
    compute_region_confidence,
)

from hmm.viterbi import viterbi, candidate_guided_tm_domains, build_alternating_path, merge_hmm_and_guided_spans
from hmm.model import load_params, TRAINED_PARAMS_PATH

from ui.plots import (
    plot_hydrophobicity,
    plot_hmm_path,
    plot_sequence_coloured,
    plot_topology_diagram,
    plot_emission_table,
    plot_transition_table,
)

# Page configs
st.set_page_config(page_title="TM Helix Predictor", layout="wide")

st.title("Transmembrane Helix Predictor")
st.caption("Kyte-Doolittle heuristics + 3-state HMM Viterbi decoding")

with st.sidebar:
    st.header("Input")

    seq_input = st.text_area(
        "Protein sequence:",
        value="MYGKIIFVLLLLSLEVLGSASTITTVVIPAVIGILVSLGVIAGTITVWRVRSSKPKSNGRHD"
              "PFKEEPVHETAPTESKSTTSKVNEISDSAIVAGVVIGLLLIIVYLFFRCLKSIAQIEESLTTQ",
        height=240,
        help="Standard single-letter amino acid codes.",
    )

    window_size = st.slider(
        "Window size:",
        min_value=5, max_value=25, value=19, step=1,
        help="Window size for KD averaging",
    )

    kd_threshold = st.slider(
        "KD threshold:",
        min_value=1.0, max_value=2.5, value=1.6, step=0.1,
        help="Mean KD score above which a window is considered hydrophobic.",
    )

    st.divider()
    predict_btn = st.button("Run", type="primary", use_container_width=True)

# Main logic
if predict_btn or "last_result" in st.session_state:

    # Validate
    seq = validate_sequence(seq_input)
    if seq is None:
        st.error(f"{validation_error_message(seq_input)}")
        st.stop()

    # Run pipeline
    with st.spinner("Running prediction..."):
        # Heuristics
        window_scores = sliding_window(seq, window_size)
        candidates = find_tm_candidates(window_scores, threshold=kd_threshold)

        # HMM
        hmm_result = viterbi(seq, min_loop_gap=2)
        hmm_path = hmm_result["path"]
        hmm_spans = hmm_result["tm_spans"]

        # Candidate-guided multi-domain refinement (local HMM per candidate).
        guided_spans = candidate_guided_tm_domains(seq, candidates, flank=10, min_loop_gap=2)
        # Merge global HMM spans with candidate-guided refined spans to avoid missing
        # domains that fall just below the KD threshold but are detected by the HMM
        display_spans = merge_hmm_and_guided_spans(hmm_spans, guided_spans, min_gap=2)
        display_path = build_alternating_path(
            len(seq),
            display_spans,
            start_loop_state=hmm_path[0] if hmm_path else "C",
        )

        confidence = compute_global_confidence(seq, candidates)
        confidence["domain_count"] = len(display_spans)
        domain_confidences = compute_domain_confidences(seq, display_spans, candidates)

        # Precompute per-domain statistics once for UI display
        display_domains = []
        for span, d_conf in zip(display_spans, domain_confidences):
            d_start = span["start"]
            d_end = span["end"]
            display_domains.append({
                "span": span,
                "confidence": d_conf,
                "composition": analyse_tm_composition(seq, d_start, d_end),
                "aromatics": find_aromatic_anchors(seq, d_start, d_end),
                "motifs": detect_tm_motifs(seq, d_start, d_end),
                "positive_inside": check_positive_inside(seq, d_start, d_end),
            })

    # Store result in session state
    st.session_state["last_result"] = {
        "seq": seq, "window_scores": window_scores, "candidates": candidates,
        "confidence": confidence,
        "domain_confidences": domain_confidences,
        "display_domains": display_domains,
        "hmm_path": hmm_path, "hmm_spans": hmm_spans,
        "display_path": display_path, "display_spans": display_spans,
        "hmm_log_prob": hmm_result["log_prob"],
        "tm_domain_count": len(display_spans),
    }

r = st.session_state.get("last_result")
if r is None:
    st.info("Enter a protein sequence to run.")
    st.stop()

seq = r["seq"]
confidence = r["confidence"]
candidates = r["candidates"]
hmm_path = r["hmm_path"]
hmm_spans = r["hmm_spans"]
display_path = r.get("display_path", hmm_path)
display_spans = r.get("display_spans", hmm_spans)
domain_confidences = r.get("domain_confidences", [])
display_domains = r.get("display_domains", [])
tm_domain_count = r.get("tm_domain_count", len(display_spans))
candidate_count = confidence.get("candidate_count", len(candidates))

# Results
col1, col2, col3, col4 = st.columns(4)
with col1:
    st.metric("Prediction", f"{confidence['label']}")
with col2:
    st.metric("Confidence score", f"{confidence['score']:.0%}")
with col3:
    st.metric("Confidence level", confidence["confidence"])
with col4:
    tm_count = len([s for s in display_path if s == "M"])
    st.metric("HMM TM domains", tm_domain_count)

st.caption(
    f"Confidence score = mean over {candidate_count} KD candidate(s). "
    f"TM domains merged when separated by <2 loop residues."
)

st.divider()

# Plot tabs
tab_hydro, tab_hmm, tab_topo, tab_evidence, tab_model = st.tabs([
    "Hydrophobicity",
    "HMM State Path",
    "Topology",
    "Statistics",
    "HMM Model",
])

with tab_hydro:
    st.subheader("Sliding-Window Hydrophobicity (Kyte-Doolittle)")
    st.caption(
        "Orange = heuristic TM candidate (mean KD ≥ threshold). "
        "Blue = HMM-predicted TM domain(s). "
        "Red = KD threshold."
    )
    plot_hydrophobicity(seq, r["window_scores"], candidates, display_spans, kd_threshold)
    st.markdown(f"**HMM domain count:** {tm_domain_count}")

    if candidates:
        best = candidates[0]
        st.markdown(
            f"**Best heuristic candidate:** residues {best['start']+1}–{best['end']} "
            f"(length {best['length']}), mean KD = **{best['mean_score']:.2f}**, "
            f"peak KD = {best['peak_score']:.2f}"
        )
        if best.get("long_flag"):
            st.warning("Candidate is longer than typical TM helices (>25 residues).")
    else:
        st.info("No hydrophobic stretch above threshold detected.")

with tab_hmm:
    st.subheader("Viterbi Decoding (Most Probable State Path)")
    st.caption(
        "The HMM assigns each residue to one of three states: "
        "Cytosolic (C), TM Helix (M), or Extracellular (E), "
        "using the most probable path through the model."
    )
    plot_hmm_path(display_path)
    st.markdown(f"**Log-probability of best path:** {r['hmm_log_prob']:.1f}")
    st.markdown(f"**TM domains found:** {tm_domain_count}")

    if display_spans:
        st.markdown("**TM spans (candidate-guided + local HMM refinement):**")
        for span in display_spans:
            seq_slice = seq[span["start"]:span["end"]]
            st.code(
                f"Residues {span['start']+1}–{span['end']}  "
                f"(length {span['length']})  {seq_slice}"
            )
    else:
        st.info("HMM did not predict a membrane-spanning segment.")

    st.subheader("Colour-coded sequence")
    plot_sequence_coloured(seq, display_path)

with tab_topo:
    st.subheader("Predicted Membrane Topology")
    st.caption(
        "Topology is determined directly from the HMM state path, "
        "allowing multi-domain C/M/E alternation."
    )
    plot_topology_diagram(seq, display_path)

    spans_for_rule = display_spans if display_spans else (candidates[:1] if candidates else [])

    for i, domain in enumerate(display_domains, start=1):
        span = domain["span"]
        d_start, d_end = span["start"], span["end"]
        d_pos_inside = domain["positive_inside"]
        label = f"TM domain {i} (residues {d_start+1}–{d_end})"
        with st.expander(label):
            st.markdown(f"""
            **Positive-inside rule**:

            | Region | Arg+Lys count |
            |--------|--------------|
            | N-terminal | {d_pos_inside['n_positive']} |
            | C-terminal | {d_pos_inside['c_positive']} |
            | Δ (N − C) | {d_pos_inside['delta']:+d} |

            **Predicted cytosolic side:** {d_pos_inside['cytosolic_side']}  
            **Orientation:** {d_pos_inside['orientation']}
            """)

    if not spans_for_rule:
        st.info("No TM domains detected. Positive-inside rule cannot be applied.")

with tab_evidence:
    st.subheader("Per-Domain Statistics")
    if not display_spans:
        if candidates:
            st.info(
                "No HMM-refined TM domains identified, but KD candidates exist."
            )
        else:
            st.info("No TM domains or KD candidates detected.")

    for i, domain in enumerate(display_domains, start=1):
        span = domain["span"]
        d_start = span["start"]
        d_end = span["end"]
        d_conf = domain["confidence"]
        d_comp = domain["composition"]
        d_aro = domain["aromatics"]
        d_motifs = domain["motifs"]

        with st.expander(
            f"TM domain {i}: residues {d_start+1}–{d_end} "
        ):
            st.metric("Domain confidence", f"{d_conf['score']:.0%}")
            st.caption(f"{d_conf['label']} — {d_conf['confidence']} confidence")

            col_a, col_b, col_c, col_d = st.columns(4)
            col_a.metric("Length", d_comp.get("length", "—"))
            col_b.metric("% Hydrophobic", f"{d_comp.get('hydrophobic_frac', 0):.0%}")
            col_c.metric("Charged residues", d_comp.get("charged_count", "—"))
            col_d.metric("Mean KD", d_comp.get("mean_kd", "—"))

            if d_comp.get("warnings"):
                for w in d_comp["warnings"]:
                    st.warning(w)
            if d_comp.get("tm_sequence"):
                st.code(f"TM sequence: {d_comp['tm_sequence']}")

            n_aro = d_aro.get("n_aromatics", [])
            c_aro = d_aro.get("c_aromatics", [])
            st.markdown("**Aromatic anchors:**")
            if n_aro or c_aro:
                st.markdown(
                    f"N-terminal interface: {', '.join(f'{aa}@{pos+1}' for pos, aa in n_aro) or 'none'} \n"
                    f"C-terminal interface: {', '.join(f'{aa}@{pos+1}' for pos, aa in c_aro) or 'none'}"
                )
                if d_aro.get("strong_belt"):
                    st.success("Tryptophan anchor present; strong membrane localisation signal.")
            else:
                st.info("No aromatic anchors detected for this domain.")

            st.markdown("**TM motifs:**")
            motif_rows = d_motifs.get("motifs", [])
            if motif_rows:
                any_hits = False
                for motif in motif_rows:
                    hits = motif.get("hits", [])
                    if hits:
                        any_hits = True
                        hit_text = ", ".join(
                            f"{h['sequence']}@{h['start']+1}-{h['end']}"
                            for h in hits
                        )
                        st.markdown(f"- **{motif['name']}**: {hit_text}")
                if not any_hits:
                    st.info("No motifs detected for this domain.")
            else:
                st.info("Motif detector unavailable for this domain.")

with tab_model:
    st.subheader("HMM Parameters")
    st.caption(
        "These probabilities drive the Viterbi decoder trained on TM protein data. "
        "To retrain: `python -m cli.train --uniprot data/uniprot_tm.txt`"
    )

    initial, transitions, emissions = load_params()

    if os.path.exists(TRAINED_PARAMS_PATH):
        st.success("Using trained parameters from `hmm/trained_params.json`")
    else:
        st.info("Using default parameters (no trained_params.json found)")

    col_e, col_t = st.columns(2)
    with col_e:
        st.markdown("**Emission probabilities** P(residue class | state)")
        plot_emission_table(emissions)
    with col_t:
        st.markdown("**Transition probabilities** P(next state | current state)")
        plot_transition_table(transitions)

    with st.expander("Initial state probabilities"):
        for s, p in initial.items():
            st.markdown(f"  **{s}**: {p:.3f}")
            