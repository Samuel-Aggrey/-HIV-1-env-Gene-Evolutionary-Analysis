# ============================================
# HIV-1 env Evolution + Functional Mapping Pipeline
# ============================================

from Bio import AlignIO
from Bio.Seq import Seq
from collections import Counter
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math

# ============================================
# STEP 1 — Load Alignment
# ============================================

def load_alignment(file_path):
    alignment = AlignIO.read(file_path, "clustal")
    print(f"Loaded {len(alignment)} sequences, length {alignment.get_alignment_length()}")
    return alignment

# ============================================
# STEP 2 — Metrics
# ============================================

def compute_entropy(column):
    counts = Counter(column)
    total = sum(counts.values())

    entropy = 0
    for count in counts.values():
        p = count / total
        entropy -= p * math.log2(p)

    return entropy

def compute_metrics(alignment):
    num_seqs = len(alignment)
    aln_length = alignment.get_alignment_length()

    conservation_scores = []
    entropy_scores = []
    gap_fractions = []

    valid_bases = {'A', 'T', 'G', 'C'}

    for i in range(aln_length):
        raw_column = list(alignment[:, i])

        # GAP FRACTION
        gap_fraction = raw_column.count('-') / num_seqs
        gap_fractions.append(gap_fraction)

        # FILTER VALID BASES
        column = [b for b in raw_column if b in valid_bases]

        if len(column) == 0:
            conservation_scores.append(0)
            entropy_scores.append(0)
            continue

        # CONSERVATION
        counts = Counter(column)
        most_common = counts.most_common(1)[0][1]
        conservation = most_common / len(column)
        conservation_scores.append(conservation)

        # ENTROPY
        entropy_scores.append(compute_entropy(column))

    return conservation_scores, entropy_scores, gap_fractions

# ============================================
# STEP 3 — Quantile-Based Categorization
# ============================================

def categorize_sites(conservation_scores):
    high_thresh = np.percentile(conservation_scores, 75)
    low_thresh = np.percentile(conservation_scores, 25)

    categories = []

    for score in conservation_scores:
        if score >= high_thresh:
            categories.append("Highly Conserved")
        elif score <= low_thresh:
            categories.append("Variable")
        else:
            categories.append("Moderately Conserved")

    print("\nThresholds:")
    print(f"High ≥ {high_thresh:.3f}")
    print(f"Low ≤ {low_thresh:.3f}")

    return categories

# ============================================
# STEP 4 — Translate Alignment (DNA → Protein)
# ============================================

def translate_alignment(alignment):
    protein_seqs = []

    for record in alignment:
        seq = str(record.seq).replace('-', '')  # remove gaps
        trimmed_len = len(seq) - (len(seq) % 3)
        seq = seq[:trimmed_len]

        protein = str(Seq(seq).translate())
        protein_seqs.append(protein)

    return protein_seqs

# ============================================
# STEP 5 — Map Positions to gp120 Regions
# ============================================

def map_to_gp120(position_nt):
    aa_pos = position_nt // 3  # convert nucleotide → amino acid

    if 131 <= aa_pos <= 157:
        return "V1 Loop"
    elif 158 <= aa_pos <= 196:
        return "V2 Loop"
    elif 296 <= aa_pos <= 331:
        return "V3 Loop"
    elif 385 <= aa_pos <= 418:
        return "V4 Loop"
    elif 460 <= aa_pos <= 469:
        return "V5 Loop"
    else:
        return "Conserved/Core"

# ============================================
# STEP 6 — Save Results
# ============================================

def save_results(conservation, entropy, gaps, categories, output="results.csv"):
    regions = [map_to_gp120(i) for i in range(len(conservation))]

    df = pd.DataFrame({
        "Position": range(1, len(conservation)+1),
        "Conservation": conservation,
        "Entropy": entropy,
        "GapFraction": gaps,
        "Category": categories,
        "gp120_region": regions
    })

    df.to_csv(output, index=False)
    print(f"Saved to {output}")

# ============================================
# STEP 7 — Plot
# ============================================

def plot_metrics(conservation, entropy, gaps):
    import matplotlib.pyplot as plt

    positions = range(1, len(conservation) + 1)

    fig, axes = plt.subplots(
        3, 1,
        figsize=(12, 8),
        sharex=True
    )

    # Define colors (colorblind-friendly palette)
    colors = {
        "conservation": "#1f77b4",  # blue
        "entropy": "#ff7f0e",  # orange
        "gaps": "#2ca02c"  # green
    }

    # =============================
    # 1. Conservation
    # =============================
    axes[0].plot(positions, conservation, linewidth=1.5, color=colors["conservation"])
    axes[0].set_ylabel("Conservation")
    axes[0].set_ylim(0, 1.05)

    axes[0].text(
        0.98, 0.85,
        "Conservation",
        transform=axes[0].transAxes,
        ha='right',
        fontsize=10,
        color=colors["conservation"]
    )

    # =============================
    # 2. Entropy
    # =============================
    axes[1].plot(positions, entropy, linewidth=1.5, color=colors["entropy"])
    axes[1].set_ylabel("Entropy")

    axes[1].text(
        0.98, 0.85,
        "Entropy",
        transform=axes[1].transAxes,
        ha='right',
        fontsize=10,
        color=colors["entropy"]
    )

    # =============================
    # 3. Gap Fraction
    # =============================
    axes[2].plot(positions, gaps, linewidth=1.5, color=colors["gaps"])
    axes[2].set_ylabel("Gap Fraction")
    axes[2].set_xlabel("Position")
    axes[2].set_ylim(0, 1)

    axes[2].text(
        0.98, 0.85,
        "Gap Fraction",
        transform=axes[2].transAxes,
        ha='right',
        fontsize=10,
        color=colors["gaps"]
    )

    # =============================
    # Styling
    # =============================
    for ax in axes:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    plt.suptitle("HIV-1 env Gene Evolutionary Profile", fontsize=14)

    plt.tight_layout()
    plt.savefig("env_evolution_profile.png", dpi=300, bbox_inches='tight')
    plt.show()

# ============================================
# MAIN
# ============================================

if __name__ == "__main__":

    alignment_file = "env_alignment.aln"

    alignment = load_alignment(alignment_file)

    conservation, entropy, gaps = compute_metrics(alignment)

    categories = categorize_sites(conservation)

    save_results(conservation, entropy, gaps, categories)

    plot_metrics(conservation, entropy, gaps)