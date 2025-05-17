from graphviz import Digraph
import re


def parse_hmm_transitions(hmm_file_path):
    with open(hmm_file_path, 'r') as file:
        lines = file.readlines()

    transitions = []
    capture = False
    state_index = 0

    def safe_float(val):
        # Convert to float, return 0.0 for '*' or any non-numeric value
        return float(val) if val != '*' else 0.0  

    for i, line in enumerate(lines):
        if line.strip().startswith("HMM"):
            capture = True
            continue
        if capture:
            if line.strip() == "//":
                break
            if re.match(r'^\s*\d+', line):  # Match lines with a numeric index
                parts = line.strip().split()
                if len(parts) >= 23:
                    state_index += 1
                    transition_probs_line = lines[i + 2].strip().split()
                    if len(transition_probs_line) >= 7:
                        transitions.append({
                            "state": f"M{state_index}",
                            "m->m": safe_float(transition_probs_line[0]),
                            "m->i": safe_float(transition_probs_line[1]),
                            "m->d": safe_float(transition_probs_line[2]),
                            "i->m": safe_float(transition_probs_line[3]),
                            "i->i": safe_float(transition_probs_line[4]),
                            "d->m": safe_float(transition_probs_line[5]),
                            "d->d": safe_float(transition_probs_line[6]),
                        })
    return transitions
def build_graph(transitions):
    dot = Digraph(format='png')
    dot.attr(rankdir='LR', size='10,5')
    dot.node("Start", shape="circle", style="filled", fillcolor="lightgray")
    dot.node("End", shape="circle", style="filled", fillcolor="lightgray")

    for i, trans in enumerate(transitions[:4]):  # Use first 4 states
        m = f"M{i+1}"
        i_state = f"I{i+1}"
        d = f"D{i+1}"

        dot.node(m, shape="square", style="filled", fillcolor="lightblue")
        dot.node(i_state, shape="diamond", style="filled", fillcolor="lightgreen")
        dot.node(d, shape="circle", style="filled", fillcolor="orange")

        if i > 0:
            prev_m = f"M{i}"
            dot.edge(prev_m, m, label=f"m->m: {trans['m->m']:.2f}")
            dot.edge(prev_m, i_state, label=f"m->i: {trans['m->i']:.2f}")
            dot.edge(prev_m, d, label=f"m->d: {trans['m->d']:.2f}")
        else:
            dot.edge("Start", m, label="Start")

        dot.edge(i_state, i_state, label=f"i->i: {trans['i->i']:.2f}")
        dot.edge(i_state, m, label=f"i->m: {trans['i->m']:.2f}")
        dot.edge(d, m, label=f"d->m: {trans['d->m']:.2f}")
        dot.edge(d, d, label=f"d->d: {trans['d->d']:.2f}")

    dot.edge("M10", "End", label="End")
    dot.render("hmm_transition_diagram2", view=True)


transitions = parse_hmm_transitions("structural_model.hmm")
build_graph(transitions)