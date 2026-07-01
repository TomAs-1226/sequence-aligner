"""
Sequence alignment engine for the tunable aligner app.

This is the generalized layer that sits on top of the assignment's global/local
alignment functions. It adds:
  - DNA scoring (match reward / mismatch penalty) AND protein scoring matrices
    (BLOSUM62, PAM250),
  - linear gaps AND affine gaps (open + extend) for the gap-problem stretch goal,
  - readouts: percent identity, gap count and gap lengths,
  - DNA -> protein translation for the SARS DNA-vs-protein comparison.

Everything is plain Python with no third-party dependencies, so it runs anywhere
(including a fresh Google Colab) with no `pip install`.
"""

from __future__ import annotations

NEG_INF = float("-inf")

# ---------------------------------------------------------------------------
# Protein scoring matrices (exact values, identical to NCBI / Biopython).
# Stored as text blocks and parsed once at import time.
# ---------------------------------------------------------------------------

_BLOSUM62_TEXT = """
   A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   *
A   4  -1  -2  -2   0  -1  -1   0  -2  -1  -1  -1  -1  -2  -1   1   0  -3  -2   0  -2  -1   0  -4
R  -1   5   0  -2  -3   1   0  -2   0  -3  -2   2  -1  -3  -2  -1  -1  -3  -2  -3  -1   0  -1  -4
N  -2   0   6   1  -3   0   0   0   1  -3  -3   0  -2  -3  -2   1   0  -4  -2  -3   3   0  -1  -4
D  -2  -2   1   6  -3   0   2  -1  -1  -3  -4  -1  -3  -3  -1   0  -1  -4  -3  -3   4   1  -1  -4
C   0  -3  -3  -3   9  -3  -4  -3  -3  -1  -1  -3  -1  -2  -3  -1  -1  -2  -2  -1  -3  -3  -2  -4
Q  -1   1   0   0  -3   5   2  -2   0  -3  -2   1   0  -3  -1   0  -1  -2  -1  -2   0   3  -1  -4
E  -1   0   0   2  -4   2   5  -2   0  -3  -3   1  -2  -3  -1   0  -1  -3  -2  -2   1   4  -1  -4
G   0  -2   0  -1  -3  -2  -2   6  -2  -4  -4  -2  -3  -3  -2   0  -2  -2  -3  -3  -1  -2  -1  -4
H  -2   0   1  -1  -3   0   0  -2   8  -3  -3  -1  -2  -1  -2  -1  -2  -2   2  -3   0   0  -1  -4
I  -1  -3  -3  -3  -1  -3  -3  -4  -3   4   2  -3   1   0  -3  -2  -1  -3  -1   3  -3  -3  -1  -4
L  -1  -2  -3  -4  -1  -2  -3  -4  -3   2   4  -2   2   0  -3  -2  -1  -2  -1   1  -4  -3  -1  -4
K  -1   2   0  -1  -3   1   1  -2  -1  -3  -2   5  -1  -3  -1   0  -1  -3  -2  -2   0   1  -1  -4
M  -1  -1  -2  -3  -1   0  -2  -3  -2   1   2  -1   5   0  -2  -1  -1  -1  -1   1  -3  -1  -1  -4
F  -2  -3  -3  -3  -2  -3  -3  -3  -1   0   0  -3   0   6  -4  -2  -2   1   3  -1  -3  -3  -1  -4
P  -1  -2  -2  -1  -3  -1  -1  -2  -2  -3  -3  -1  -2  -4   7  -1  -1  -4  -3  -2  -2  -1  -2  -4
S   1  -1   1   0  -1   0   0   0  -1  -2  -2   0  -1  -2  -1   4   1  -3  -2  -2   0   0   0  -4
T   0  -1   0  -1  -1  -1  -1  -2  -2  -1  -1  -1  -1  -2  -1   1   5  -2  -2   0  -1  -1   0  -4
W  -3  -3  -4  -4  -2  -2  -3  -2  -2  -3  -2  -3  -1   1  -4  -3  -2  11   2  -3  -4  -3  -2  -4
Y  -2  -2  -2  -3  -2  -1  -2  -3   2  -1  -1  -2  -1   3  -3  -2  -2   2   7  -1  -3  -2  -1  -4
V   0  -3  -3  -3  -1  -2  -2  -3  -3   3   1  -2   1  -1  -2  -2   0  -3  -1   4  -3  -2  -1  -4
B  -2  -1   3   4  -3   0   1  -1   0  -3  -4   0  -3  -3  -2   0  -1  -4  -3  -3   4   1  -1  -4
Z  -1   0   0   1  -3   3   4  -2   0  -3  -3   1  -1  -3  -1   0  -1  -3  -2  -2   1   4  -1  -4
X   0  -1  -1  -1  -2  -1  -1  -1  -1  -1  -1  -1  -1  -1  -2   0   0  -2  -1  -1  -1  -1  -1  -4
*  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4   1
"""

_PAM250_TEXT = """
   A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   *
A   2  -2   0   0  -2   0   0   1  -1  -1  -2  -1  -1  -3   1   1   1  -6  -3   0   0   0   0  -8
R  -2   6   0  -1  -4   1  -1  -3   2  -2  -3   3   0  -4   0   0  -1   2  -4  -2  -1   0  -1  -8
N   0   0   2   2  -4   1   1   0   2  -2  -3   1  -2  -3   0   1   0  -4  -2  -2   2   1   0  -8
D   0  -1   2   4  -5   2   3   1   1  -2  -4   0  -3  -6  -1   0   0  -7  -4  -2   3   3  -1  -8
C  -2  -4  -4  -5  12  -5  -5  -3  -3  -2  -6  -5  -5  -4  -3   0  -2  -8   0  -2  -4  -5  -3  -8
Q   0   1   1   2  -5   4   2  -1   3  -2  -2   1  -1  -5   0  -1  -1  -5  -4  -2   1   3  -1  -8
E   0  -1   1   3  -5   2   4   0   1  -2  -3   0  -2  -5  -1   0   0  -7  -4  -2   3   3  -1  -8
G   1  -3   0   1  -3  -1   0   5  -2  -3  -4  -2  -3  -5   0   1   0  -7  -5  -1   0   0  -1  -8
H  -1   2   2   1  -3   3   1  -2   6  -2  -2   0  -2  -2   0  -1  -1  -3   0  -2   1   2  -1  -8
I  -1  -2  -2  -2  -2  -2  -2  -3  -2   5   2  -2   2   1  -2  -1   0  -5  -1   4  -2  -2  -1  -8
L  -2  -3  -3  -4  -6  -2  -3  -4  -2   2   6  -3   4   2  -3  -3  -2  -2  -1   2  -3  -3  -1  -8
K  -1   3   1   0  -5   1   0  -2   0  -2  -3   5   0  -5  -1   0   0  -3  -4  -2   1   0  -1  -8
M  -1   0  -2  -3  -5  -1  -2  -3  -2   2   4   0   6   0  -2  -2  -1  -4  -2   2  -2  -2  -1  -8
F  -3  -4  -3  -6  -4  -5  -5  -5  -2   1   2  -5   0   9  -5  -3  -3   0   7  -1  -4  -5  -2  -8
P   1   0   0  -1  -3   0  -1   0   0  -2  -3  -1  -2  -5   6   1   0  -6  -5  -1  -1   0  -1  -8
S   1   0   1   0   0  -1   0   1  -1  -1  -3   0  -2  -3   1   2   1  -2  -3  -1   0   0   0  -8
T   1  -1   0   0  -2  -1   0   0  -1   0  -2   0  -1  -3   0   1   3  -5  -3   0   0  -1   0  -8
W  -6   2  -4  -7  -8  -5  -7  -7  -3  -5  -2  -3  -4   0  -6  -2  -5  17   0  -6  -5  -6  -4  -8
Y  -3  -4  -2  -4   0  -4  -4  -5   0  -1  -1  -4  -2   7  -5  -3  -3   0  10  -2  -3  -4  -2  -8
V   0  -2  -2  -2  -2  -2  -2  -1  -2   4   2  -2   2  -1  -1  -1   0  -6  -2   4  -2  -2  -1  -8
B   0  -1   2   3  -4   1   3   0   1  -2  -3   1  -2  -4  -1   0   0  -5  -3  -2   3   2  -1  -8
Z   0   0   1   3  -5   3   3   0   2  -2  -3   0  -2  -5   0   0  -1  -6  -4  -2   2   3  -1  -8
X   0  -1   0  -1  -3  -1  -1  -1  -1  -1  -1  -1  -1  -2  -1   0   0  -4  -2  -1  -1  -1  -1  -8
*  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8  -8   1
"""


def _parse_matrix(text: str) -> dict[tuple[str, str], int]:
    """Turn an NCBI-style matrix text block into a {(a, b): score} dictionary."""
    lines = [ln for ln in text.strip().splitlines() if ln.strip()]
    columns = lines[0].split()
    scores: dict[tuple[str, str], int] = {}
    for ln in lines[1:]:
        parts = ln.split()
        row_letter = parts[0]
        for col_letter, value in zip(columns, parts[1:]):
            scores[(row_letter, col_letter)] = int(value)
    return scores


BLOSUM62 = _parse_matrix(_BLOSUM62_TEXT)
PAM250 = _parse_matrix(_PAM250_TEXT)


# ---------------------------------------------------------------------------
# Scoring: a "scorer" is just a function score(a, b) -> number.
# ---------------------------------------------------------------------------

def dna_scorer(match: float, mismatch: float):
    """Returns a scoring function for DNA: +match if equal, -mismatch if not."""
    def score(a: str, b: str) -> float:
        return match if a == b else -mismatch
    return score


def matrix_scorer(matrix: dict[tuple[str, str], int]):
    """Returns a scoring function that looks pairs up in a substitution matrix.

    Unknown letters fall back to the matrix's 'X' (any) row/column."""
    def score(a: str, b: str) -> float:
        if (a, b) in matrix:
            return matrix[(a, b)]
        a = a if (a, "A") in matrix else "X"
        b = b if ("A", b) in matrix else "X"
        return matrix[(a, b)]
    return score


# ---------------------------------------------------------------------------
# Linear-gap alignment (one flat penalty per gap symbol).
# ---------------------------------------------------------------------------

def align_global(seq1: str, seq2: str, scorer, gap: float):
    """Needleman-Wunsch global alignment with a linear gap penalty.

    Returns (row1, row2, score)."""
    n, m = len(seq1), len(seq2)
    s = [[0.0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        s[i][0] = s[i - 1][0] - gap
    for j in range(1, m + 1):
        s[0][j] = s[0][j - 1] - gap
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = s[i - 1][j - 1] + scorer(seq1[i - 1], seq2[j - 1])
            up = s[i - 1][j] - gap
            left = s[i][j - 1] - gap
            s[i][j] = max(diag, up, left)

    row1, row2, _, _ = _traceback_linear(seq1, seq2, scorer, gap, s, n, m, local=False)
    return row1, row2, s[n][m]


def align_local(seq1: str, seq2: str, scorer, gap: float):
    """Smith-Waterman local alignment with a linear gap penalty.

    Returns (row1, row2, score, start1, end1, start2, end2)."""
    n, m = len(seq1), len(seq2)
    s = [[0.0] * (m + 1) for _ in range(n + 1)]
    best, bi, bj = 0.0, 0, 0
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = s[i - 1][j - 1] + scorer(seq1[i - 1], seq2[j - 1])
            up = s[i - 1][j] - gap
            left = s[i][j - 1] - gap
            s[i][j] = max(0.0, diag, up, left)
            if s[i][j] > best:
                best, bi, bj = s[i][j], i, j

    row1, row2, si, sj = _traceback_linear(seq1, seq2, scorer, gap, s, bi, bj, local=True)
    return row1, row2, best, si, bi, sj, bj


def _traceback_linear(seq1, seq2, scorer, gap, s, i, j, local):
    eps = 1e-9
    row1, row2 = [], []
    while (i > 0 or j > 0) and not (local and s[i][j] <= eps):
        if i > 0 and j > 0 and abs(s[i][j] - (s[i - 1][j - 1] + scorer(seq1[i - 1], seq2[j - 1]))) < eps:
            row1.append(seq1[i - 1]); row2.append(seq2[j - 1]); i -= 1; j -= 1
        elif i > 0 and abs(s[i][j] - (s[i - 1][j] - gap)) < eps:
            row1.append(seq1[i - 1]); row2.append("-"); i -= 1
        else:
            row1.append("-"); row2.append(seq2[j - 1]); j -= 1
    return "".join(reversed(row1)), "".join(reversed(row2)), i, j


# ---------------------------------------------------------------------------
# Affine-gap alignment (Gotoh): opening a gap costs `gap_open`, each extra
# symbol in the same gap costs `gap_extend`. This is the gap-problem stretch.
# ---------------------------------------------------------------------------

def align_global_affine(seq1: str, seq2: str, scorer, gap_open: float, gap_extend: float):
    """Global alignment with affine gaps. Returns (row1, row2, score)."""
    n, m = len(seq1), len(seq2)
    M = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
    Ix = [[NEG_INF] * (m + 1) for _ in range(n + 1)]  # gap in seq2 (consume seq1)
    Iy = [[NEG_INF] * (m + 1) for _ in range(n + 1)]  # gap in seq1 (consume seq2)
    M[0][0] = 0.0
    for i in range(1, n + 1):
        Ix[i][0] = -gap_open - (i - 1) * gap_extend
    for j in range(1, m + 1):
        Iy[0][j] = -gap_open - (j - 1) * gap_extend
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            sub = scorer(seq1[i - 1], seq2[j - 1])
            M[i][j] = max(M[i - 1][j - 1], Ix[i - 1][j - 1], Iy[i - 1][j - 1]) + sub
            Ix[i][j] = max(M[i - 1][j] - gap_open, Ix[i - 1][j] - gap_extend)
            Iy[i][j] = max(M[i][j - 1] - gap_open, Iy[i][j - 1] - gap_extend)

    row1, row2 = _traceback_affine(seq1, seq2, scorer, gap_open, gap_extend, M, Ix, Iy, n, m)
    return row1, row2, max(M[n][m], Ix[n][m], Iy[n][m])


def _traceback_affine(seq1, seq2, scorer, gap_open, gap_extend, M, Ix, Iy, i, j):
    eps = 1e-9
    row1, row2 = [], []
    # choose the starting matrix (which one holds the best end score)
    state = max(("M", M[i][j]), ("Ix", Ix[i][j]), ("Iy", Iy[i][j]), key=lambda t: t[1])[0]
    while i > 0 or j > 0:
        if state == "M":
            row1.append(seq1[i - 1]); row2.append(seq2[j - 1])
            sub = scorer(seq1[i - 1], seq2[j - 1])
            prev = M[i][j] - sub
            i -= 1; j -= 1
            if abs(prev - M[i][j]) < eps:
                state = "M"
            elif abs(prev - Ix[i][j]) < eps:
                state = "Ix"
            else:
                state = "Iy"
        elif state == "Ix":
            row1.append(seq1[i - 1]); row2.append("-")
            if i >= 1 and abs(Ix[i][j] - (M[i - 1][j] - gap_open)) < eps:
                state = "M"
            else:
                state = "Ix"
            i -= 1
        else:  # Iy
            row1.append("-"); row2.append(seq2[j - 1])
            if j >= 1 and abs(Iy[i][j] - (M[i][j - 1] - gap_open)) < eps:
                state = "M"
            else:
                state = "Iy"
            j -= 1
    return "".join(reversed(row1)), "".join(reversed(row2))


# ---------------------------------------------------------------------------
# Semi-global alignment (free end gaps): overhangs at the start and end of
# either sequence are not penalized. Good when one sequence sits inside or
# overlaps the other, so global alignment's end-gap penalty is unfair.
# ---------------------------------------------------------------------------

def align_semiglobal(seq1: str, seq2: str, scorer, gap: float):
    """Global-style alignment where leading and trailing gaps are free.

    Returns (row1, row2, score) with the full sequences shown, but the score
    only counts the aligned core (end overhangs are not penalized)."""
    n, m = len(seq1), len(seq2)
    s = [[0.0] * (m + 1) for _ in range(n + 1)]  # first row/col stay 0 -> free leading gaps
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = s[i - 1][j - 1] + scorer(seq1[i - 1], seq2[j - 1])
            up = s[i - 1][j] - gap
            left = s[i][j - 1] - gap
            s[i][j] = max(diag, up, left)

    # best end is anywhere in the last row or last column (free trailing gaps)
    best, bi, bj = s[n][0], n, 0
    for j in range(m + 1):
        if s[n][j] >= best:
            best, bi, bj = s[n][j], n, j
    for i in range(n + 1):
        if s[i][m] > best:
            best, bi, bj = s[i][m], i, m

    eps = 1e-9
    core1, core2 = [], []
    i, j = bi, bj
    while i > 0 and j > 0:
        if abs(s[i][j] - (s[i - 1][j - 1] + scorer(seq1[i - 1], seq2[j - 1]))) < eps:
            core1.append(seq1[i - 1]); core2.append(seq2[j - 1]); i -= 1; j -= 1
        elif abs(s[i][j] - (s[i - 1][j] - gap)) < eps:
            core1.append(seq1[i - 1]); core2.append("-"); i -= 1
        else:
            core1.append("-"); core2.append(seq2[j - 1]); j -= 1
    si, sj = i, j
    core1.reverse(); core2.reverse()

    # add the unaligned overhangs as free end gaps
    lead1 = seq1[:si] + "-" * sj
    lead2 = "-" * si + seq2[:sj]
    tail1 = seq1[bi:] + "-" * (m - bj)
    tail2 = "-" * (n - bi) + seq2[bj:]
    return lead1 + "".join(core1) + tail1, lead2 + "".join(core2) + tail2, best


def reverse_complement(dna: str) -> str:
    """Reverse complement of a DNA string (A<->T, C<->G)."""
    comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N",
            "a": "t", "t": "a", "c": "g", "g": "c", "n": "n", "U": "A", "u": "a"}
    return "".join(comp.get(c, c) for c in reversed(dna))


# ---------------------------------------------------------------------------
# Readouts
# ---------------------------------------------------------------------------

def percent_identity(row1: str, row2: str) -> float:
    """Fraction of aligned columns (excluding gap columns) that are exact matches."""
    matches = aligned = 0
    for a, b in zip(row1, row2):
        if a == "-" or b == "-":
            continue
        aligned += 1
        if a == b:
            matches += 1
    return 100.0 * matches / aligned if aligned else 0.0


def gap_stats(row1: str, row2: str) -> tuple[int, list[int]]:
    """Returns (number of separate gaps, list of each gap's length)."""
    lengths: list[int] = []
    run = 0
    for a, b in zip(row1, row2):
        if a == "-" or b == "-":
            run += 1
        else:
            if run:
                lengths.append(run)
            run = 0
    if run:
        lengths.append(run)
    return len(lengths), lengths


def alignment_score(row1: str, row2: str, scorer, gap: float) -> float:
    """Recompute an alignment's score under a linear gap model (a sanity check)."""
    total = 0.0
    for a, b in zip(row1, row2):
        if a == "-" or b == "-":
            total -= gap
        else:
            total += scorer(a, b)
    return total


def format_alignment(row1: str, row2: str, width: int = 60) -> str:
    """Plain-text alignment block with a match line (| match, . mismatch, space gap)."""
    out = []
    for start in range(0, len(row1), width):
        a = row1[start:start + width]
        b = row2[start:start + width]
        mid = "".join(
            "|" if x == y and x != "-" else (" " if x == "-" or y == "-" else ".")
            for x, y in zip(a, b)
        )
        out.append(a)
        out.append(mid)
        out.append(b)
        out.append("")
    return "\n".join(out)


# ---------------------------------------------------------------------------
# DNA -> protein translation (standard genetic code)
# ---------------------------------------------------------------------------

_CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def translate(dna: str, stop_at_stop: bool = True) -> str:
    """Translate a DNA coding sequence into one-letter amino acids."""
    dna = dna.upper().replace("U", "T")
    protein = []
    for i in range(0, len(dna) - 2, 3):
        aa = _CODON_TABLE.get(dna[i:i + 3], "X")
        if aa == "*":
            if stop_at_stop:
                break
            aa = "*"
        protein.append(aa)
    return "".join(protein)


# ---------------------------------------------------------------------------
# FASTA reading
# ---------------------------------------------------------------------------

def read_fasta(path: str) -> str:
    """Read a single-record FASTA file and return its sequence (header stripped)."""
    seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith(">"):
                seq.append(line)
    return "".join(seq)


def read_multi_fasta(path: str) -> tuple[list[str], list[str]]:
    """Read a multi-record FASTA file. Returns (names, sequences)."""
    names: list[str] = []
    seqs: list[str] = []
    cur: list[str] = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if names:
                    seqs.append("".join(cur))
                names.append(line[1:])
                cur = []
            else:
                cur.append(line)
    if names:
        seqs.append("".join(cur))
    return names, seqs
