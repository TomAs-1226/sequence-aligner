"""Render a static 3D image of the SARS-CoV-2 spike trimer, colored by which
residues changed vs SARS-CoV. This renders as a normal image, so it shows in
Colab, in the exported HTML, and locally (unlike the live py3Dmol viewer)."""
import os
import urllib.request
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
import engine as e

OUT = os.path.join("..", "figures")
os.makedirs(OUT, exist_ok=True)

# 1. which residues changed (from our own aligner)
p1 = e.translate(e.read_fasta("../Alignment/Data/Coronaviruses/SARS-CoV_genome_spike_protein.fasta"))
p2 = e.translate(e.read_fasta("../Alignment/Data/Coronaviruses/SARS-CoV-2_genome_spike_protein.fasta"))
pr1, pr2, _ = e.align_global(p1, p2, e.matrix_scorer(e.BLOSUM62), 10)
changed = set()
resi2 = 0
for a, b in zip(pr1, pr2):
    if b == "-":
        continue
    resi2 += 1
    if a == "-" or a != b:
        changed.add(resi2)

# 2. get the structure (cache locally)
pdb_path = "6vxx.pdb"
if not os.path.exists(pdb_path):
    req = urllib.request.Request("https://files.rcsb.org/download/6VXX.pdb",
                                 headers={"User-Agent": "Mozilla/5.0"})
    open(pdb_path, "w").write(urllib.request.urlopen(req, timeout=30).read().decode("utf-8", "ignore"))

# 3. parse CA atoms per chain
chains = {}
for line in open(pdb_path):
    if line.startswith("ATOM") and line[12:16].strip() == "CA":
        ch = line[21]
        resi = int(line[22:26])
        xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        chains.setdefault(ch, []).append((resi, xyz))

# 4. draw two views of the trimer
fig = plt.figure(figsize=(11, 6))
for k, (elev, azim, label) in enumerate([(15, 0, "side view"), (88, 0, "top view (looking down the spike)")]):
    ax = fig.add_subplot(1, 2, k + 1, projection="3d")
    for atoms in chains.values():
        xs = [a[1][0] for a in atoms]
        ys = [a[1][1] for a in atoms]
        zs = [a[1][2] for a in atoms]
        ax.plot(xs, ys, zs, color="#d3dad7", linewidth=1.0)
        rx = [a[1][0] for a in atoms if a[0] in changed]
        ry = [a[1][1] for a in atoms if a[0] in changed]
        rz = [a[1][2] for a in atoms if a[0] in changed]
        ax.scatter(rx, ry, rz, color="#d32f2f", s=9, depthshade=False)
    ax.view_init(elev=elev, azim=azim)
    ax.set_axis_off()
    ax.set_title(label, fontsize=12)
fig.suptitle("SARS-CoV-2 spike trimer  |  red = the 297 amino acids that changed from SARS-CoV,  gray = conserved",
             fontsize=13, y=0.97)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(os.path.join(OUT, "fig6_spike_3d.png"), dpi=160, bbox_inches="tight")
plt.close(fig)
print("wrote figures/fig6_spike_3d.png; changed residues:", len(changed),
      "; chains:", list(chains.keys()), "; CA per chain:", [len(v) for v in chains.values()])
