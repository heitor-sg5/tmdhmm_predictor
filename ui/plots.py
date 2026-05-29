import plotly.graph_objects as go
import streamlit as st

from data.residues import KD_SCORES

# Colour pallete
STATE_COLOURS = {
    "C":"#5b9bd5", # blue — cytosolic
    "M":"#e8872a", # amber — membrane
    "E":"#55a868", # green — extracellular
}

STATE_LABELS = {
    "C": "Cytosolic",
    "M": "TM Helix",
    "E": "Extracellular",
}

def plot_hydrophobicity(sequence, window_scores, candidates, hmm_spans, threshold=1.6):
    """
    Sliding-window KD hydrophobicity plot.
    """
    x = list(range(len(window_scores)))

    fig = go.Figure()

    # Per-residue KD bar chart
    per_res = [KD_SCORES.get(aa, 0.0) for aa in sequence]
    fig.add_trace(go.Bar(
        x=x, y=per_res,
        name="Per-residue KD",
        marker_color=["#e8872a" if v >= 1.8 else "#aecde8" for v in per_res],
        opacity=0.35,
        showlegend=True,
    ))

    # Sliding window line
    fig.add_trace(go.Scatter(
        x=x, y=window_scores,
        mode="lines",
        name="Sliding window (w=19)",
        line=dict(color="#1a4f8a", width=2.5),
    ))

    # Threshold line
    fig.add_hline(
        y=threshold,
        line_dash="dash",
        line_color="red",
        line_width=1.5,
        annotation_text=f"KD threshold ({threshold})",
        annotation_position="top right",
    )

    # Shade heuristic TM candidates
    for cand in candidates:
        fig.add_vrect(
            x0=cand["start"], x1=cand["end"],
            fillcolor="#e8872a", opacity=0.18,
            layer="below", line_width=0,
            annotation_text=f"KD candidate (μ={cand['mean_score']:.2f})",
            annotation_position="top left",
        )

    # Outline HMM-predicted TM spans
    for idx, span in enumerate(hmm_spans, start=1):
        fig.add_vrect(
            x0=span["start"], x1=span["end"],
            fillcolor="rgba(0,0,0,0)", opacity=1.0,
            layer="above", line_width=2, line_color="#5b9bd5",
            annotation_text=f"HMM TM{idx}",
            annotation_position="bottom right",
        )

    fig.update_layout(
        title="Sliding-Window Hydrophobicity (Kyte-Doolittle)",
        xaxis_title="Residue position",
        yaxis_title="Mean KD score",
        height=360,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        margin=dict(t=60),
    )
    st.plotly_chart(fig, use_container_width=True)

def plot_hmm_path(path):
    """
    Visualise the Viterbi state path as a colour-coded step plot.
    """
    state_y = {"C": 0, "M": 1, "E": 2}
    y_vals  = [state_y[s] for s in path]
    colours = [STATE_COLOURS[s] for s in path]

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=list(range(len(path))),
        y=y_vals,
        mode="lines+markers",
        line=dict(color="#888", width=1, shape="hv"),
        marker=dict(color=colours, size=5, line=dict(width=0)),
        name="Viterbi path",
    ))

    # State legend
    for state, y in state_y.items():
        fig.add_trace(go.Scatter(
            x=[None], y=[None],
            mode="markers",
            marker=dict(size=10, color=STATE_COLOURS[state]),
            name=STATE_LABELS[state],
            showlegend=True,
        ))

    fig.update_layout(
        title="HMM Viterbi State Path",
        xaxis_title="Residue position",
        yaxis=dict(
            tickvals=[0, 1, 2],
            ticktext=["Cytosolic (C)", "TM Helix (M)", "Extracellular (E)"],
        ),
        height=320,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        margin=dict(t=60),
    )
    st.plotly_chart(fig, use_container_width=True)

def plot_sequence_coloured(sequence, path):
    """
    Colour-coded sequence display with each residue highlighted by its predicted state.
    """
    bg_colours = {
        "C": "#aed6f1",
        "M": "#f0b27a",
        "E": "#a9dfbf",
    }
    html = (
        "<div style='font-family: monospace; font-size: 15px; "
        "line-height: 2.2; word-break: break-all; padding: 8px;'>"
    )
    for i, (aa, state) in enumerate(zip(sequence, path)):
        bg = bg_colours[state]
        tip = f"Pos {i+1}: {aa} ({STATE_LABELS[state]})"
        html += (
            f"<span title='{tip}' style='background:{bg}; color:#222; "
            f"padding: 2px 5px; margin: 1px; border-radius: 4px; "
            f"display: inline-block;'>{aa}</span>"
        )
    html += "</div>"

    # Legend
    legend = " &nbsp; ".join(
        f"<span style='background:{bg_colours[s]};padding:2px 8px;"
        f"border-radius:4px'>{STATE_LABELS[s]}</span>"
        for s in ["C", "M", "E"]
    )
    st.markdown(
        f"<div style='font-size:13px;margin-bottom:6px'>{legend}</div>",
        unsafe_allow_html=True,
    )
    st.markdown(html, unsafe_allow_html=True)

def plot_topology_diagram(sequence, path):
    """
    Simple topology sketch of cytosolic/TM/extracellular regions as labelled segments.
    """
    if not path:
        st.info("No HMM state path available for topology plot.")
        return

    # Build segment list: (state, start, end)
    segments = []
    current_state = path[0]
    seg_start = 0
    for i, s in enumerate(path[1:], 1):
        if s != current_state:
            segments.append((current_state, seg_start, i))
            current_state = s
            seg_start = i
    segments.append((current_state, seg_start, len(path)))

    # Assign y positions: C=bottom, M=middle (gradient), E=top
    y_map = {"C": 0.1, "M": 0.5, "E": 0.9}

    fig = go.Figure()

    # Draw membrane band
    fig.add_shape(type="rect", x0=0, x1=len(sequence), y0=0.35, y1=0.65, fillcolor="#ffe8c0", opacity=0.5, line_width=0, layer="below")
    fig.add_annotation(x=len(sequence) * 0.02, y=0.5, text="Lipid bilayer", showarrow=False, font=dict(size=11, color="#a0522d"))

    # Draw segments as coloured lines
    for state, start, end in segments:
        mid_x = (start + end) / 2
        y = y_map[state]
        fig.add_trace(go.Scatter(
            x=[start, end],
            y=[y, y],
            mode="lines",
            line=dict(color=STATE_COLOURS[state], width=8),
            name=STATE_LABELS[state],
            showlegend=False,
        ))
        label = f"{STATE_LABELS[state]}\n{start+1}–{end}"
        fig.add_annotation(x=mid_x, y=y + 0.07, text=label, showarrow=False, font=dict(size=10))

    # Side labels and orientation from decoded sequence endpoints.
    start_state = path[0] if path else "C"
    end_state = path[-1] if path else "E"
    fig.add_annotation(
        x=0, y=-0.05,
        text=f"N-terminus: {STATE_LABELS.get(start_state, start_state)}",
        showarrow=False, font=dict(size=10, color="#444"),
    )
    fig.add_annotation(
        x=len(sequence), y=-0.05,
        text=f"C-terminus: {STATE_LABELS.get(end_state, end_state)}",
        showarrow=False, font=dict(size=10, color="#444"),
    )

    fig.update_layout(
        title="Predicted Membrane Topology",
        xaxis=dict(title="Residue position", showgrid=False),
        yaxis=dict(visible=False, range=[-0.15, 1.15]),
        height=300,
        showlegend=False,
        margin=dict(t=50, b=40),
    )
    st.plotly_chart(fig, use_container_width=True)

def plot_emission_table(emissions: dict) -> None:
    """Show HMM emission probabilities as a colour-scaled table."""
    import pandas as pd
    rows = []
    for state in ["C", "M", "E"]:
        rows.append({
            "State": STATE_LABELS[state],
            "H (hydrophobic)": round(emissions[state]["H"], 3),
            "A (aromatic)": round(emissions[state]["A"], 3),
            "P (polar)": round(emissions[state]["P"], 3),
            "Q (charged)": round(emissions[state]["Q"], 3),
        })
    df = pd.DataFrame(rows).set_index("State")
    st.dataframe(df.style.background_gradient(cmap="YlOrRd", axis=None), use_container_width=True)

def plot_transition_table(transitions: dict) -> None:
    """Show HMM transition probabilities as a table."""
    import pandas as pd
    rows = []
    for s in ["C", "M", "E"]:
        rows.append({
            "From \\ To": STATE_LABELS[s],
            "→ Cytosolic (C)": round(transitions[s]["C"], 3),
            "→ TM Helix (M)": round(transitions[s]["M"], 3),
            "→ Extracellular (E)": round(transitions[s]["E"], 3),
        })
    df = pd.DataFrame(rows).set_index("From \\ To")
    st.dataframe(df.style.background_gradient(cmap="Blues", axis=None), use_container_width=True)