"""Generate presentation figures (PNGs) from the real engine output."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import engine as e

OUT = os.path.join("..", "figures")
os.makedirs(OUT, exist_ok=True)
D = os.path.join("..", "Alignment", "Data", "Coronaviruses")
H = os.path.join("..", "Alignment", "Data", "Hemoglobin")

GREEN, ORANGE, GRAY = "#bfe3c0", "#ffd08a", "#e0e0e0"


def color_for(a, b):
    if a == "-" or b == "-":
        return GRAY
    return GREEN if a == b else ORANGE


def draw_alignment(ax, row1, row2, title, maxcols=None):
    if maxcols:
        row1, row2 = row1[:maxcols], row2[:maxcols]
    n = len(row1)
    ax.set_xlim(0, n)
    ax.set_ylim(0, 2)
    ax.axis("off")
    ax.set_title(title, fontsize=11, loc="left")
    for i, (x, y) in enumerate(zip(row1, row2)):
        c = color_for(x, y)
        ax.add_patch(Rectangle((i, 1), 1, 1, facecolor=c, edgecolor="white", linewidth=0.4))
        ax.add_patch(Rectangle((i, 0), 1, 1, facecolor=c, edgecolor="white", linewidth=0.4))
        if n <= 70:
            ax.text(i + 0.5, 1.5, x, ha="center", va="center", fontsize=8, family="monospace")
            ax.text(i + 0.5, 0.5, y, ha="center", va="center", fontsize=8, family="monospace")


# ---------- Figure 1: toy good vs bad ----------
seq1, seq2 = "GTTTTGCTTGGAT", "GTTTTGATATTAGTTGGAT"
g1, g2, gs = e.align_global(seq1, seq2, e.dna_scorer(2, 2), 3)
b1, b2, bs = e.align_global(seq1, seq2, e.dna_scorer(2, 2), 0.1)
fig, axes = plt.subplots(2, 1, figsize=(9, 3.4))
draw_alignment(axes[0], g1, g2,
               f"good parameters (gap = 3):  score {gs:g},  identity {e.percent_identity(g1,g2):.1f}%,  "
               f"gaps {e.gap_stats(g1,g2)[0]}  ->  one clean insertion")
draw_alignment(axes[1], b1, b2,
               f"bad parameters (gap = 0.1):  score {bs:g},  identity {e.percent_identity(b1,b2):.1f}%,  "
               f"gaps {e.gap_stats(b1,b2)[0]}  ->  scattered, not real")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig1_toy_good_bad.png"), dpi=160, bbox_inches="tight")
plt.close(fig)

# ---------- Figure 2: HBB good vs bad (first 60 cols) ----------
human = e.read_fasta(os.path.join(H, "Homo_sapiens_hemoglobin.fasta"))
zeb = e.read_fasta(os.path.join(H, "Danio_rerio_hemoglobin.fasta"))
hg1, hg2, hgs = e.align_global(human, zeb, e.matrix_scorer(e.BLOSUM62), 10)
hb1, hb2, hbs = e.align_global(human, zeb, e.matrix_scorer(e.BLOSUM62), 0)
fig, axes = plt.subplots(2, 1, figsize=(9.5, 3.6))
draw_alignment(axes[0], hg1, hg2,
               f"good (BLOSUM62, gap = 10):  score {hgs:g},  identity {e.percent_identity(hg1,hg2):.1f}%,  "
               f"gap runs {e.gap_stats(hg1,hg2)[0]}   (first 45 columns)", maxcols=45)
draw_alignment(axes[1], hb1, hb2,
               f"bad (BLOSUM62, gap = 0):  score {hbs:g},  identity {e.percent_identity(hb1,hb2):.1f}%,  "
               f"gap runs {e.gap_stats(hb1,hb2)[0]}   shredded into many gaps", maxcols=45)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig2_hbb_good_bad.png"), dpi=170, bbox_inches="tight")
plt.close(fig)

# ---------- Figure 3: SARS codon breakdown ----------
g1n = e.read_fasta(os.path.join(D, "SARS-CoV_genome_spike_protein.fasta"))
g2n = e.read_fasta(os.path.join(D, "SARS-CoV-2_genome_spike_protein.fasta"))
p1, p2 = e.translate(g1n), e.translate(g2n)
pr1, pr2, _ = e.align_global(p1, p2, e.matrix_scorer(e.BLOSUM62), 10)
i1 = i2 = 0
identical = synonymous = nonsyn = 0
for a, b in zip(pr1, pr2):
    if a == "-" or b == "-":
        if a != "-": i1 += 1
        if b != "-": i2 += 1
        continue
    c1, c2 = g1n[i1*3:i1*3+3], g2n[i2*3:i2*3+3]
    if c1 == c2: identical += 1
    elif a == b: synonymous += 1
    else: nonsyn += 1
    i1 += 1; i2 += 1
fig, ax = plt.subplots(figsize=(7, 4))
labels = ["identical\ncodon", "synonymous\n(silent in protein)", "nonsynonymous\n(changes protein)"]
vals = [identical, synonymous, nonsyn]
colors = ["#9ec7e8", "#bfe3c0", "#ef9a9a"]
bars = ax.bar(labels, vals, color=colors, edgecolor="white")
for bar, v in zip(bars, vals):
    ax.text(bar.get_x()+bar.get_width()/2, v+6, str(v), ha="center", fontsize=12)
ax.set_ylabel("codon columns")
ax.set_title("SARS-CoV vs SARS-CoV-2 spike: where DNA changes hide\n"
             f"DNA identity 77.3%   vs   protein identity 78.0%   "
             f"({synonymous} of {synonymous+nonsyn} DNA changes are silent)", fontsize=10)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig3_sars_codons.png"), dpi=160, bbox_inches="tight")
plt.close(fig)

# ---------- Figure 4: spike domain map with changed residues ----------
i2 = 0
changed = []
for a, b in zip(pr1, pr2):
    if b == "-":
        continue
    i2 += 1
    if a == "-" or a != b:
        changed.append(i2)
L = len(p2)
domains = [("NTD", 14, 305, "#c5cae9"), ("RBD", 319, 541, "#ffcc80"),
           ("S1/S2", 542, 685, "#e1bee7"), ("S2 fusion", 686, 1273, "#b2dfdb")]
fig, ax = plt.subplots(figsize=(13, 3.0))
ax.set_xlim(-10, L + 10); ax.set_ylim(0, 3.2); ax.axis("off")
for r in changed:
    ax.plot([r, r], [2.15, 2.78], color="#d32f2f", linewidth=0.7)
ax.text(0, 2.98, f"Red ticks: the {len(changed)} amino acids that changed "
        f"({100*len(changed)/L:.0f}% of the spike)", fontsize=15, weight="bold", color="#222")
for name, s, en, col in domains:
    ax.add_patch(Rectangle((s, 0.7), en - s, 1.15, facecolor=col, edgecolor="white"))
    ax.text((s + en) / 2, 1.27, name, ha="center", va="center", fontsize=16, weight="bold", color="#222")
ax.text(0, 0.35, "Functional regions along the 1273-amino-acid spike chain", fontsize=13, color="#444")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig4_spike_domains.png"), dpi=170, bbox_inches="tight")
plt.close(fig)

# ---------- Figure 5: gap problem linear vs affine ----------
anc = "ACGTACGTAACCGGTTACGTACGTACGTTTGGCCAATTACGTACGTGGCCATTAACG"
desc = "ACGTCCGTAACCGTTTACGTACTTACGTCCCCCAAAATTGGCGAATTACGTAGGTGGCCATAAACG"
l1, l2, ls = e.align_global(anc, desc, e.dna_scorer(2, 3), 1)
a1, a2, asc = e.align_global_affine(anc, desc, e.dna_scorer(2, 3), 8, 1)
fig, axes = plt.subplots(2, 1, figsize=(11, 3.4))
draw_alignment(axes[0], l1, l2,
               f"linear gap = 1:  {e.gap_stats(l1,l2)[0]} separate gaps (the one real insertion gets shattered)")
draw_alignment(axes[1], a1, a2,
               f"affine (open 8, extend 1):  {e.gap_stats(a1,a2)[0]} gap (one clean insertion kept together)")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig5_gap_problem.png"), dpi=160, bbox_inches="tight")
plt.close(fig)

# ---------- Figure 7: DP scoring grid with the best path ----------
def dp_grid(s1, s2, match=1, mismatch=1, gap=1):
    n, m = len(s1), len(s2)
    S = [[0.0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        S[i][0] = -i * gap
    for j in range(1, m + 1):
        S[0][j] = -j * gap
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = S[i - 1][j - 1] + (match if s1[i - 1] == s2[j - 1] else -mismatch)
            S[i][j] = max(diag, S[i - 1][j] - gap, S[i][j - 1] - gap)
    i, j = n, m
    path = [(i, j)]
    while i > 0 or j > 0:
        if i > 0 and j > 0 and S[i][j] == S[i - 1][j - 1] + (match if s1[i - 1] == s2[j - 1] else -mismatch):
            i, j = i - 1, j - 1
        elif i > 0 and S[i][j] == S[i - 1][j] - gap:
            i -= 1
        else:
            j -= 1
        path.append((i, j))
    return np.array(S), path

s1d, s2d = "GCATGCU", "GATTACA"
Sg, pathg = dp_grid(s1d, s2d)
fig, ax = plt.subplots(figsize=(6.2, 5.8))
ax.imshow(Sg, cmap="Blues")
for i in range(Sg.shape[0]):
    for j in range(Sg.shape[1]):
        ax.text(j, i, int(Sg[i, j]), ha="center", va="center", fontsize=11, color="#333")
ax.plot([p[1] for p in pathg], [p[0] for p in pathg], color="#d32f2f", linewidth=2.4, marker="o", markersize=6)
ax.set_xticks(range(len(s2d) + 1)); ax.set_xticklabels(["-"] + list(s2d), fontsize=12)
ax.set_yticks(range(len(s1d) + 1)); ax.set_yticklabels(["-"] + list(s1d), fontsize=12)
ax.set_title("the scoring grid, and the best path back through it", fontsize=12)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig7_dp_grid.png"), dpi=170, bbox_inches="tight")
plt.close(fig)

# ---------- Figure 8: codon bar + conservation along the spike ----------
i1 = i2 = 0
iden = syn = non = 0
for a, b in zip(pr1, pr2):
    if a == "-" or b == "-":
        if a != "-": i1 += 1
        if b != "-": i2 += 1
        continue
    if g1n[i1*3:i1*3+3] == g2n[i2*3:i2*3+3]: iden += 1
    elif a == b: syn += 1
    else: non += 1
    i1 += 1; i2 += 1
cols = list(zip(pr1, pr2))
resi_at = []
r = 0
for a, b in cols:
    if b != "-": r += 1
    resi_at.append(r)
pos, ident = [], []
for start in range(0, len(cols) - 30 + 1, 4):
    both = [(a, b) for a, b in cols[start:start+30] if a != "-" and b != "-"]
    if both:
        pos.append(resi_at[start]); ident.append(100*sum(1 for a, b in both if a == b)/len(both))
fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
axes[0].bar(["identical", "synonymous", "nonsynonymous"], [iden, syn, non],
            color=["#9ec7e8", "#bfe3c0", "#ef9a9a"], edgecolor="white")
for i, v in enumerate([iden, syn, non]):
    axes[0].text(i, v + 5, str(v), ha="center", fontsize=11)
axes[0].set_ylabel("codon columns"); axes[0].set_title("where the DNA changes hide")
axes[0].spines["top"].set_visible(False); axes[0].spines["right"].set_visible(False)
axes[1].axvspan(319, 541, color="#ffe0b2", alpha=0.6, label="receptor-binding domain")
axes[1].plot(pos, ident, color="#2E7D46", linewidth=1.3)
axes[1].set_xlabel("position along the spike (amino acid number)")
axes[1].set_ylabel("percent identity (window of 30)")
axes[1].set_title("which parts of the spike stayed the same")
axes[1].legend(loc="lower right", fontsize=9)
axes[1].spines["top"].set_visible(False); axes[1].spines["right"].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig8_conservation.png"), dpi=160, bbox_inches="tight")
plt.close(fig)

# ---------- Figure 10: conservation line on its own ----------
fig, ax = plt.subplots(figsize=(8, 3.7))
ax.axvspan(319, 541, color="#ffe0b2", alpha=0.6, label="receptor-binding domain")
ax.plot(pos, ident, color="#2E7D46", linewidth=1.4)
ax.set_xlabel("position along the spike (amino acid number)")
ax.set_ylabel("percent identity (window of 30)")
ax.set_title("which parts of the spike stayed the same")
ax.legend(loc="lower right", fontsize=9)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig10_conservation_line.png"), dpi=160, bbox_inches="tight")
plt.close(fig)

# ---------- Figure 9: dot plot ----------
def dot_points(a, b, k=2):
    seen = {}
    for j in range(len(b) - k + 1):
        seen.setdefault(b[j:j + k], []).append(j)
    xs, ys = [], []
    for i in range(len(a) - k + 1):
        for j in seen.get(a[i:i + k], []):
            xs.append(j); ys.append(i)
    return xs, ys

dx, dy = dot_points(human, zeb, 2)
fig, ax = plt.subplots(figsize=(5.6, 5.6))
ax.scatter(dx, dy, s=5, color="#2E7D46")
ax.set_xlabel("zebrafish HBB position")
ax.set_ylabel("human HBB position")
ax.set_title("dot plot: matching 2-letter windows")
ax.invert_yaxis()
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig9_dotplot.png"), dpi=160, bbox_inches="tight")
plt.close(fig)

print("wrote figures:")
for f in sorted(os.listdir(OUT)):
    print("  figures/" + f, os.path.getsize(os.path.join(OUT, f)), "bytes")
print(f"\ncodon counts: identical={identical} synonymous={synonymous} nonsyn={nonsyn}; changed residues={len(changed)} ({100*len(changed)/L:.1f}%)")
