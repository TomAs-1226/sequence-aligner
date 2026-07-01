# Sequence Aligner: plain-English documentation

This explains, with no jargon, what the app does, how it does it, and what it gives you back.
If you can read this page, you can explain the whole project.

---

## 1. The big idea in one paragraph

You have two strings of letters. They are either DNA (`A C G T`) or protein (the 20 amino-acid
letters). They are similar but not identical: one has extra letters, missing letters, or swapped
letters compared to the other. Aligning them means sliding them against each other and inserting
dashes (`-`) until the matching letters line up as well as possible. The app finds the best
possible lineup for a scoring rule you choose, and then tells you how good that lineup is.

Example:

```
GCA-TGCU
G-ATTACA
```

The dash means the other sequence has a letter here that this one does not.

---

## 2. The three knobs you control

Every alignment is a trade-off, and you control it with three numbers:

| Knob | What it means | Effect when you raise it |
|------|---------------|--------------------------|
| **Match reward** | Points you gain each time two letters are the same | Encourages stacking up matches |
| **Mismatch penalty** | Points you lose when two different letters are paired | Discourages pairing unlike letters |
| **Gap penalty** | Points you lose every time you insert a dash | Discourages dashes (shifting things around) |

The whole point of the project is that changing these three numbers changes the answer. There is
no single correct alignment. There is only the best alignment for the scoring you picked. Picking
good numbers is a biology judgment, not just a math one.

---

## 3. How it actually finds the best alignment (the method)

The app uses a technique called dynamic programming. You do not need the textbook name. Here is
the intuition:

1. Make a grid. One sequence runs down the left side, the other across the top.
2. Each box in the grid asks: what is the best score for lining up everything up to here?
3. You fill the grid one box at a time. Each box only looks at three neighbors (the box above,
   the box to the left, and the box diagonally up-left) and picks the best of three choices:
   match or mismatch (go diagonal), gap in one sequence (go down), or gap in the other (go right).
4. The bottom-right box holds the score of the best possible alignment.
5. Then you walk backwards from that box, retracing the choices you made, and that path is the
   alignment.

This always finds the real best-scoring alignment. It never just guesses.

### Two modes

- **Global alignment** (Needleman-Wunsch): line up the sequences end to end, all of both. Use
  this when the two sequences are roughly the same thing, like the same gene in two species.
- **Local alignment** (Smith-Waterman): find the single best-matching stretch and ignore the
  ragged ends. Same grid, with two changes: no box is ever allowed to go below zero, and you start
  the traceback from the highest box anywhere in the grid instead of the corner. Use this when you
  are hunting for a shared region buried inside two otherwise different sequences.

---

## 4. What the app gives you back (the outputs)

For any alignment, you get:

- **The alignment itself.** The two sequences with dashes inserted, stacked so you can read down
  each column. A middle line marks each column: `|` is a match, `.` is a mismatch, a space is a gap.
- **The score.** The single number the scoring rule produced. Higher is better, but it is only
  comparable between runs that used the same three knobs.
- **Percent identity.** Of the lined-up columns, the fraction that are exact matches. This is the
  number biologists trust, because it does not depend on your knob settings.
- **Gap count and gap lengths.** How many separate runs of dashes there are, and how long each is.
  One long gap usually means a chunk was inserted or deleted, which is normal in biology. Many tiny
  scattered gaps usually mean the scoring is bad.

---

## 4b. The pictures the notebook draws

The notebook does not just print letters. It draws several pictures so the results are easy to see:

- **Colored alignment.** Each column is shaded: green for a match, orange for a mismatch, gray for a gap.
- **Scoring matrix heatmap.** The protein scoring table (BLOSUM62) drawn as a grid of colored numbers.
- **The scoring grid with its path.** The dynamic programming grid for a short example, with the best path drawn back through it in red. This is the algorithm, shown as a picture.
- **Dot plot.** A dot wherever the two sequences share a short piece. Related sequences make a diagonal line.
- **Codon breakdown and conservation.** A bar chart of identical, synonymous, and nonsynonymous codons, plus a line showing percent identity along the spike so you can see which parts changed.
- **The spike in 3D.** The 297 changed amino acids drawn on the real 3D spike structure, both as a plain picture and as a viewer you can rotate in Colab.

## 5. Why score and percent identity are different (this is the key insight)

The score is whatever your knobs say it is. Turn the match reward up high enough and any alignment
scores great, so the number is flattering itself. A high score does not prove a good alignment.

Percent identity, gap structure, and whether the aligned regions are biologically meaningful are
what tell you if the alignment is real. That is the project's central question: how do you judge
quality when the score depends on your own settings? The answer is that you judge it with biology
(identity, sensible gaps, known conserved regions), not with the score alone.

---

## 6. DNA vs protein (the SARS comparison)

The same gene can be aligned two ways:

- as DNA (the raw `A C G T` letters), or
- as the protein it codes for (translate every 3 DNA letters into 1 amino-acid letter).

The interesting result: many DNA differences disappear when you compare the proteins. That is
because the genetic code is redundant. Several different 3-letter DNA codons spell the same amino
acid. A DNA change that swaps one codon for another that means the same amino acid is called a
synonymous (silent) mutation. It shows up in the DNA alignment but vanishes in the protein
alignment. Where the two alignments disagree, you are seeing exactly which mutations changed the
protein and which were silent. That tells you which changes biology cares about.

---

## 6c. Extra features (what changed and how to explain it)

The engine and notebook grew some new abilities. Here is each one in plain words, so you can explain them when presenting.

- **Semi-global alignment (free end gaps).** Normal global alignment charges you for gaps at the very start and end. If one sequence just has extra letters hanging off the front or back, that is unfair and gives a bad (even negative) score. Semi-global alignment makes those end gaps free, so only the middle overlap is scored. Use it when one sequence overlaps or sits inside the other. In code: `align_semiglobal(seq1, seq2, scorer, gap)`.
- **Reverse complement.** DNA has two strands. `reverse_complement(dna)` flips a sequence to its partner strand (A pairs with T, C with G, and the order reverses). Handy when a match is hiding on the opposite strand.
- **Banded alignment (for very long sequences).** The normal method needs a grid as big as one sequence times the other. For two 30,000-letter genomes that is almost a billion boxes, which runs out of memory. Banded alignment only fills a narrow strip near the diagonal of the grid, because for similar sequences the best path stays close to the diagonal. It gives the same answer while using a tiny fraction of the memory, and it aligns two whole genomes in about 20 to 30 seconds. In code: `align_banded(seq1, seq2, scorer, gap, band)`.

New notebook parts, all built on the same engine:

- **Parameter sweep.** Slides the gap penalty and plots the score, the percent identity, and the number of gaps together, so you can see the score and the biology disagree.
- **Is it better than chance?** Shuffles a sequence 200 times and compares the real score to the shuffled ones. A z-score far above 3 means the alignment is real, not luck.
- **A family tree.** Turns the pairwise alignments of the four hemoglobins into distances and joins the closest ones (UPGMA). The tree matches biology: human and gorilla together, then cow, with zebrafish as the outsider.
- **The big one.** Aligns two whole 30,000-letter coronavirus genomes with the banded method and plots conservation across the whole genome.

---

## 7. What each code file does

| File | Plain-English job |
|------|-------------------|
| `globalAlignmentScore.py` | Builds the scoring grid for end-to-end alignment |
| `globalAlignment.py` | Walks the grid backwards to produce the end-to-end alignment |
| `localAlignmentScore.py` | Builds the scoring grid for best-stretch alignment (floored at 0) |
| `localAlignment.py` | Finds the best box, walks back to produce the best-stretch alignment |
| `editDistance*.py` | Counts the fewest edits (swap, insert, delete) to turn one string into another |
| `lcs*.py` | Finds the longest run of letters two strings share in order |
| `sharedKMers.py` | Counts short chunks two strings have in common |
| `io_python.py` | Reads sequence files (FASTA format) and writes alignments out |
| `helperFunctions.py` | Small shared utilities |
| `manhattanTourist.py` | A grid-path warm-up exercise (not used by the aligner) |

The `Data/` folder holds the real sequences (hemoglobin in several animals, the SARS-CoV and
SARS-CoV-2 genomes, and their spike proteins). The `Output/` folder is where results get written.

---

## 8. How to run the core right now

From inside the `Alignment` folder:

```python
from globalAlignment import global_alignment
from localAlignment import local_alignment

# match reward = 1, mismatch penalty = 1, gap penalty = 1
print(global_alignment("GCATGCU", "GATTACA", 1, 1, 1))
print(local_alignment("GGTTGACTA", "TGTTACGG", 1, 1, 1))
```

That is the engine. Everything else we build on top (the tunable app, the readouts, the visuals,
the protein scoring matrices) sits on top of these functions.
