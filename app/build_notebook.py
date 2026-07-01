"""
Assembles the final Colab notebook (SequenceAligner.ipynb) from:
  - the validated engine (engine.py), embedded as a code cell,
  - the real pack sequences (HBB proteins + SARS spike genes), baked into a data cell,
  - the agent-built/verified sections in sections.json (keyed by section id).

Run:  python build_notebook.py
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
PACK = os.path.join(HERE, "..", "Alignment", "Data")


def read_fasta(path):
    seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith(">"):
                seq.append(line)
    return "".join(seq)


def md(source):
    return {"cell_type": "markdown", "metadata": {}, "source": source.splitlines(keepends=True)}


def code(source):
    return {
        "cell_type": "code",
        "metadata": {},
        "execution_count": None,
        "outputs": [],
        "source": source.splitlines(keepends=True),
    }


# ---- read the engine source ----
with open(os.path.join(HERE, "engine.py")) as f:
    engine_src = f.read()
# the notebook cell does not need the module docstring's file reference; keep as-is.

# ---- read the real sequences ----
H = os.path.join(PACK, "Hemoglobin")
C = os.path.join(PACK, "Coronaviruses")
HBB_HUMAN = read_fasta(os.path.join(H, "Homo_sapiens_hemoglobin.fasta"))
HBB_GORILLA = read_fasta(os.path.join(H, "Gorilla_gorilla_hemoglobin.fasta"))
HBB_COW = read_fasta(os.path.join(H, "Bos_taurus_hemoglobin.fasta"))
HBB_ZEBRAFISH = read_fasta(os.path.join(H, "Danio_rerio_hemoglobin.fasta"))
SARS_COV_SPIKE_DNA = read_fasta(os.path.join(C, "SARS-CoV_genome_spike_protein.fasta"))
SARS_COV2_SPIKE_DNA = read_fasta(os.path.join(C, "SARS-CoV-2_genome_spike_protein.fasta"))


def wrap(seq, width=70):
    return "\n".join('    "%s"' % seq[i:i + width] for i in range(0, len(seq), width))


data_cell = f'''# ============================================================
# The real sequences (baked in, so this notebook runs anywhere
# with no file uploads). HBB = hemoglobin-beta protein per species;
# SARS spike = the spike GENE nucleotides for each virus.
# ============================================================

HBB_HUMAN = "{HBB_HUMAN}"
HBB_GORILLA = "{HBB_GORILLA}"
HBB_COW = "{HBB_COW}"
HBB_ZEBRAFISH = "{HBB_ZEBRAFISH}"

SARS_COV_SPIKE_DNA = (
{wrap(SARS_COV_SPIKE_DNA)}
)
SARS_COV2_SPIKE_DNA = (
{wrap(SARS_COV2_SPIKE_DNA)}
)

# Translate the spike GENES into the spike PROTEINS (standard genetic code).
SARS_COV_SPIKE_PROT = translate(SARS_COV_SPIKE_DNA)
SARS_COV2_SPIKE_PROT = translate(SARS_COV2_SPIKE_DNA)

print("HBB lengths (aa):", len(HBB_HUMAN), len(HBB_GORILLA), len(HBB_COW), len(HBB_ZEBRAFISH))
print("SARS-CoV  spike: ", len(SARS_COV_SPIKE_DNA), "nt ->", len(SARS_COV_SPIKE_PROT), "aa")
print("SARS-CoV-2 spike:", len(SARS_COV2_SPIKE_DNA), "nt ->", len(SARS_COV2_SPIKE_PROT), "aa")
'''

# ---- load agent sections ----
sections_path = os.path.join(HERE, "sections.json")
if os.path.exists(sections_path):
    with open(sections_path, encoding="utf-8") as f:
        sections = json.load(f)
else:
    sections = {}


def section_cells(section_id, fallback_title):
    data = sections.get(section_id)
    if not data or not data.get("cells"):
        return [md(f"> _Section `{section_id}` not yet assembled._")]
    out = []
    for c in data["cells"]:
        if c["kind"] == "markdown":
            out.append(md(c["source"]))
        else:
            out.append(code(c["source"]))
    return out


cells = []

cells.append(md(
    "# Build and Tune Your Own Sequence Aligner\n"
    "**Module 2 Project - CMU Pre-College Computational Biology**\n\n"
    "This notebook builds a small, *tunable* sequence-alignment app from my own global and "
    "local alignment code, then uses it to make an argument about what separates a good "
    "alignment from a bad one.\n\n"
    "**How it's organized**\n\n"
    "1. **The engine** - my alignment code (the algorithms I wrote this module, generalized).\n"
    "2. **The sequences** - the real DNA/protein I'll align.\n"
    "3. **The interactive app** - sliders for match / mismatch / gap, global vs local, DNA vs protein.\n"
    "4. **How the alignment is built** - pictures of the scoring grid and a dot plot.\n"
    "5. **Good vs bad on a toy** - the same two sequences, two parameter sets, two very different answers.\n"
    "6. **Good vs bad on a real gene (HBB)** - the same story with real biology.\n"
    "7. **DNA vs protein (SARS spike)** - why synonymous mutations show up in DNA but vanish in protein.\n"
    "8. **What makes an alignment good?** - the central question, answered with biology.\n"
    "9. **Stretch: the gap problem** - why one gap penalty isn't enough (and affine gaps).\n"
    "10. **Stretch: sanity check** - does my aligner agree with a professional tool?\n"
    "11. **Stretch: the mutations in 3D** - the spike changes drawn on the real structure.\n"
))

cells.append(md(
    "## Part 1 - The engine (my alignment code)\n\n"
    "This is the dynamic-programming alignment I wrote in the module. The DNA path "
    "(match reward / mismatch penalty / one gap penalty) is exactly the global and local "
    "alignment from class; I generalized it so it can *also* score proteins with a substitution "
    "matrix (BLOSUM62 / PAM250) and use affine gaps. It has no external dependencies, so "
    "**Run All** works in a fresh Colab.\n\n"
    "Run this cell first."
))
cells.append(code(engine_src))

cells.append(md("## Part 2 - The sequences"))
cells.append(code(data_cell))

cells.append(md("## Part 3 - The interactive app"))
cells.extend(section_cells("ui_app", "Interactive app"))

cells.append(md("## Part 4 - How the alignment is built (pictures)"))
cells.extend(section_cells("viz_algorithm", "Algorithm pictures"))

cells.append(md("## Part 5 - Good vs bad parameters: a toy example"))
cells.extend(section_cells("toy_goodbad", "Toy good vs bad"))

cells.append(md("## Part 6 - Good vs bad parameters: a real gene (hemoglobin beta)"))
cells.extend(section_cells("hbb_goodbad", "HBB good vs bad"))

cells.append(md("## Part 7 - DNA vs protein: the SARS spike"))
cells.extend(section_cells("sars_dna_protein", "SARS DNA vs protein"))

cells.append(md("## Part 8 - What makes an alignment good? (the central question)"))
cells.extend(section_cells("quality_essay", "Quality essay"))

cells.append(md("## Part 9 (stretch) - The gap problem"))
cells.extend(section_cells("gap_problem", "Gap problem"))

cells.append(md("## Part 10 (stretch) - Sanity check vs a professional tool"))
cells.extend(section_cells("sanity_check", "Sanity check"))

cells.append(md("## Part 11 (stretch) - Seeing the mutations in 3D"))
cells.extend(section_cells("spike_3d", "3D spike viewer"))

cells.append(md(
    "## Wrap-up\n\n"
    "Every number in this notebook came from my own engine in Part 1. The big lesson: "
    "**changing the scoring changes the alignment**, so the raw score alone can't tell you "
    "which alignment is right - biology has to (percent identity, sensible gaps, conserved "
    "regions, and agreement with independent evidence)."
))

notebook = {
    "cells": cells,
    "metadata": {
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python", "version": "3.12"},
        "colab": {"provenance": []},
    },
    "nbformat": 4,
    "nbformat_minor": 5,
}

out_path = os.path.join(HERE, "..", "SequenceAligner.ipynb")
with open(out_path, "w", encoding="utf-8") as f:
    json.dump(notebook, f, indent=1, ensure_ascii=False)

print("Wrote", os.path.abspath(out_path), "with", len(cells), "cells.")
print("Sections present:", sorted(sections.keys()) if sections else "(none yet)")
