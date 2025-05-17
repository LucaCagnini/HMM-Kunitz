import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logomaker
import math

AA_ORDER = ['A','C','D','E','F','G','H','I','K','L',
            'M','N','P','Q','R','S','T','V','W','Y']

def parse_hmm_emissions(hmm_path):
    emissions = []  # emission list
    with open(hmm_path) as f:
        lines = f.readlines()
        #print(lines)

    in_matrix = False
    for line in lines:
        if line.strip().startswith("HMM"):
            in_matrix = True
            continue
        if in_matrix:
            if line.strip() == "//":
                break
            if line.strip() == "" or line.startswith("  COMPO") or line.strip().startswith("m->"):
                continue
            if line.strip().split()[0].isdigit():
                parts = line.strip().split()
                emission_scores = parts[1:21]  # values
                probs = []
                for score in emission_scores:
                    if score == '*':
                        probs.append(0.0)
                    else:
                        val = float(score)
                        prob = math.exp(-val)  # Conversion -ln(p) in p
                        probs.append(prob)

                total = sum(probs)
                if total > 0:
                    probs = [p / total for p in probs]
                else:
                    probs = [0.0] * 20

                emissions.append(probs)

    if len(emissions) == 0:
        raise ValueError("⚠️ Nessuna emissione trovata. Controlla se il file HMM è corretto.")

    df = pd.DataFrame(emissions, columns=AA_ORDER).astype(float)
    return df

def compute_information_content(df):
    """
    Converte la matrice di frequenza in bits (contenuto informativo).
    """
    ic_matrix = []
    for i, row in df.iterrows():
        row_probs = row.values.astype(float)
        entropy = -sum([p * math.log2(p) for p in row_probs if p > 0])
        max_entropy = math.log2(20)  # 4.32 bit for 20 aa
        info = max_entropy - entropy
        ic_row = [p * info for p in row_probs]
        ic_matrix.append(ic_row)

    ic_df = pd.DataFrame(ic_matrix, columns=AA_ORDER).fillna(0.0).astype(float)
    return ic_df

def plot_logo(df, output_file='hmm_logo_bits.png'):
    plt.figure(figsize=(len(df)/2.5, 3.5))
    custom_colors = {
    'A': 'green',
    'C': 'blue',
    'D': 'red',
    'E': 'red',
    'F': 'purple',
    'G': 'orange',
    'H': 'pink',
    'I': 'brown',
    'K': 'cyan',
    'L': 'brown',
    'M': 'olive',
    'N': 'gray',
    'P': 'yellow',
    'Q': 'gray',
    'R': 'cyan',
    'S': 'gold',
    'T': 'gold',
    'V': 'brown',
    'W': 'purple',
    'Y': 'purple',
    '-': 'white',
}
    logomaker.Logo(df, shade_below=.5, fade_below=.5, color_scheme = custom_colors)
    plt.xticks(rotation=90)
    plt.ylabel("Bits")
    plt.ylim(0, 4.5)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()
    print(f"Sequence logo in bits salvato in {output_file}")

def main():
    hmm_file = 'structural_model.hmm'  
    freqs = parse_hmm_emissions(hmm_file)
    bits = compute_information_content(freqs)
    plot_logo(bits)

if __name__ == '__main__':
    main()
