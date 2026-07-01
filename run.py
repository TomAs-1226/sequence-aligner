"""
One-click runner. Runs every analysis and writes ALL results into the Output
folder: alignment text files, a summary, and every chart as a PNG.

Usage:  python run.py      (or just double-click run.bat on Windows)
"""
import os
import sys
import shutil
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
APP = os.path.join(HERE, "app")
sys.path.insert(0, APP)
import engine as e  # noqa: E402

OUT = os.path.join(HERE, "Output")
FIG_OUT = os.path.join(OUT, "figures")
os.makedirs(FIG_OUT, exist_ok=True)

CORONA = os.path.join(HERE, "Alignment", "Data", "Coronaviruses")
HEMO = os.path.join(HERE, "Alignment", "Data", "Hemoglobin")

summary = []


def log(line=""):
    print(line)
    summary.append(line)


def write_alignment(path, title, row1, row2, score, scorer, gap):
    ng, gl = e.gap_stats(row1, row2)
    with open(path, "w", encoding="utf-8") as f:
        f.write(title + "\n" + "=" * len(title) + "\n\n")
        f.write(e.format_alignment(row1, row2) + "\n")
        f.write(f"score            = {score:g}\n")
        f.write(f"percent identity = {e.percent_identity(row1, row2):.1f}%\n")
        f.write(f"gaps             = {ng}  lengths {gl}\n")
    return ng, gl


log("Sequence aligner: running everything and writing to the Output folder.")
log("")

# ---------------- Toy good vs bad ----------------
log("[1] Toy DNA example (good vs bad parameters)")
s1, s2 = "GTTTTGCTTGGAT", "GTTTTGATATTAGTTGGAT"
g = e.align_global(s1, s2, e.dna_scorer(2, 2), 3)
b = e.align_global(s1, s2, e.dna_scorer(2, 2), 0.1)
ng_g, _ = write_alignment(os.path.join(OUT, "toy_good.txt"), "Toy, GOOD params (gap=3)", g[0], g[1], g[2], e.dna_scorer(2, 2), 3)
ng_b, _ = write_alignment(os.path.join(OUT, "toy_bad.txt"), "Toy, BAD params (gap=0.1)", b[0], b[1], b[2], e.dna_scorer(2, 2), 0.1)
log(f"    good: score {g[2]:g}, identity {e.percent_identity(g[0],g[1]):.1f}%, gaps {ng_g}")
log(f"    bad : score {b[2]:g}, identity {e.percent_identity(b[0],b[1]):.1f}%, gaps {ng_b}")
log("")

# ---------------- HBB good vs bad ----------------
log("[2] Hemoglobin beta protein, human vs zebrafish (good vs bad)")
human = e.read_fasta(os.path.join(HEMO, "Homo_sapiens_hemoglobin.fasta"))
zeb = e.read_fasta(os.path.join(HEMO, "Danio_rerio_hemoglobin.fasta"))
hg = e.align_global(human, zeb, e.matrix_scorer(e.BLOSUM62), 10)
hb = e.align_global(human, zeb, e.matrix_scorer(e.BLOSUM62), 0)
ng_hg, _ = write_alignment(os.path.join(OUT, "hbb_good.txt"), "HBB, GOOD (BLOSUM62, gap=10)", hg[0], hg[1], hg[2], e.matrix_scorer(e.BLOSUM62), 10)
ng_hb, _ = write_alignment(os.path.join(OUT, "hbb_bad.txt"), "HBB, BAD (BLOSUM62, gap=0)", hb[0], hb[1], hb[2], e.matrix_scorer(e.BLOSUM62), 0)
log(f"    good: score {hg[2]:g}, identity {e.percent_identity(hg[0],hg[1]):.1f}%, gap runs {ng_hg}")
log(f"    bad : score {hb[2]:g}, identity {e.percent_identity(hb[0],hb[1]):.1f}%, gap runs {ng_hb}")
log("")

# ---------------- SARS DNA vs protein ----------------
log("[3] SARS-CoV vs SARS-CoV-2 spike (DNA and protein)")
g1n = e.read_fasta(os.path.join(CORONA, "SARS-CoV_genome_spike_protein.fasta"))
g2n = e.read_fasta(os.path.join(CORONA, "SARS-CoV-2_genome_spike_protein.fasta"))
p1, p2 = e.translate(g1n), e.translate(g2n)
d = e.align_global(g1n, g2n, e.dna_scorer(1, 1), 2)
pr = e.align_global(p1, p2, e.matrix_scorer(e.BLOSUM62), 10)
write_alignment(os.path.join(OUT, "sars_dna.txt"), "SARS spike DNA alignment", d[0], d[1], d[2], e.dna_scorer(1, 1), 2)
write_alignment(os.path.join(OUT, "sars_protein.txt"), "SARS spike protein alignment", pr[0], pr[1], pr[2], e.matrix_scorer(e.BLOSUM62), 10)
i1 = i2 = 0
iden = syn = non = 0
for a, bb in zip(pr[0], pr[1]):
    if a == "-" or bb == "-":
        if a != "-": i1 += 1
        if bb != "-": i2 += 1
        continue
    if g1n[i1*3:i1*3+3] == g2n[i2*3:i2*3+3]: iden += 1
    elif a == bb: syn += 1
    else: non += 1
    i1 += 1; i2 += 1
log(f"    DNA     identity {e.percent_identity(d[0],d[1]):.1f}%, gaps {e.gap_stats(d[0],d[1])[0]}")
log(f"    protein identity {e.percent_identity(pr[0],pr[1]):.1f}%, gaps {e.gap_stats(pr[0],pr[1])[0]}")
log(f"    codons: identical {iden}, synonymous {syn}, nonsynonymous {non}")
log("")

# ---------------- Gap problem ----------------
log("[4] Gap problem: linear vs affine")
anc = "ACGTACGTAACCGGTTACGTACGTACGTTTGGCCAATTACGTACGTGGCCATTAACG"
desc = "ACGTCCGTAACCGTTTACGTACTTACGTCCCCCAAAATTGGCGAATTACGTAGGTGGCCATAAACG"
lin = e.align_global(anc, desc, e.dna_scorer(2, 3), 1)
aff = e.align_global_affine(anc, desc, e.dna_scorer(2, 3), 8, 1)
write_alignment(os.path.join(OUT, "gap_linear.txt"), "Gap problem, linear gap=1 (match 2, mismatch 3)", lin[0], lin[1], lin[2], e.dna_scorer(2, 3), 1)
with open(os.path.join(OUT, "gap_affine.txt"), "w", encoding="utf-8") as f:
    f.write("Gap problem, affine (open 8, extend 1)\n\n")
    f.write(e.format_alignment(aff[0], aff[1]) + "\n")
    f.write(f"gaps = {e.gap_stats(aff[0],aff[1])}\n")
log(f"    linear gap=1: {e.gap_stats(lin[0],lin[1])[0]} separate gaps")
log(f"    affine      : {e.gap_stats(aff[0],aff[1])[0]} gap")
log("")

# ---------------- Figures ----------------
log("[5] Generating charts...")
for script in ("make_figures.py", "make_3d.py"):
    try:
        subprocess.run([sys.executable, script], cwd=APP, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    except Exception as ex:
        log(f"    (skipped {script}: {ex})")
src_fig = os.path.join(HERE, "figures")
if os.path.isdir(src_fig):
    for f in os.listdir(src_fig):
        if f.endswith(".png"):
            shutil.copy(os.path.join(src_fig, f), os.path.join(FIG_OUT, f))
    log(f"    copied {len([f for f in os.listdir(FIG_OUT) if f.endswith('.png')])} charts to Output/figures")
log("")

with open(os.path.join(OUT, "summary.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(summary))

log("Done. Open the Output folder to see everything.")
print("\nOutput folder:", OUT)
