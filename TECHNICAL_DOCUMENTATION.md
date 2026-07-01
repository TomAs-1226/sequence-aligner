# Sequence Aligner: full technical documentation

This document explains every algorithm and model in the project, in detail, for both
the Python engine and the Windows app. It is written so a student can read the code and
explain exactly what each part does, why it works, and what it costs. A shorter plain
English overview lives in `Alignment/DOCUMENTATION.md`; this is the deep reference.

## Contents
1. The alignment algorithms (global, local, semi-global, affine, banded)
2. Supporting algorithms and data (edit distance, LCS, k-mers, translation, scoring matrices, readouts)
3. GPU accelerated alignment
4. The on-device AI models (secondary structure, protein language model)
5. The Windows app (architecture and features)

---

# Part 1. The alignment algorithms


All of the aligners in `app/engine.py` are built on one idea: **dynamic programming (DP)**. Aligning two sequences means choosing, at every position, whether to line up a letter from `seq1` with a letter from `seq2` (a diagonal move), skip a letter of `seq2` by putting a gap in `seq1` (a left move), or skip a letter of `seq1` by putting a gap in `seq2` (an up move). There are exponentially many ways to interleave these choices, so instead of trying them all we fill a table `s` where `s[i][j]` is the best score achievable for the first `i` letters of `seq1` against the first `j` letters of `seq2`. Each cell is computed once from a few already-solved smaller cells (its neighbors up, left, and diagonal), so the whole table costs `O(n*m)` work and stores the answer to every subproblem. Once the table is full, we walk backwards from the finish cell to the start (the **traceback**), and every step we take reconstructs one column of the final alignment. All five aligners below share this skeleton; they differ only in how a cell is scored, where the alignment is allowed to start and end, and how many tables are needed.

Throughout the code a **scorer** is just a function `score(a, b)` that returns a number for lining up letter `a` with letter `b`. For DNA, `dna_scorer(match, mismatch)` returns `match` when `a == b` and `-mismatch` otherwise. For proteins, `matrix_scorer(matrix)` looks the pair up in `BLOSUM62` or `PAM250`. The `gap` parameter is a **positive** number: gaps and mismatches are penalties that get *subtracted*, and matches are rewards that get *added*. Keep that sign convention in mind because it is the reason you see `+ scorer(...)` but `- gap` everywhere.

---

## 1. Global alignment (Needleman-Wunsch): `align_global`

**What it does.** Aligns the two sequences *end to end*. Every letter of `seq1` and every letter of `seq2` must appear in the result, so gaps at the very start or end are still charged. This is the right choice when the two inputs are believed to be the same thing from one end to the other (for example two homologous genes of similar length).

**The table.** `s` has shape `(n+1) x (m+1)` where `n = len(seq1)` and `m = len(seq2)`. `s[i][j]` is the best score aligning `seq1[:i]` with `seq2[:j]`.

**Initialization.** The first row and column are pure gaps, so they accumulate the linear penalty. Because `gap` is positive and subtracted:

```
s[0][0] = 0
s[i][0] = s[i-1][0] - gap      # i letters of seq1 all aligned to gaps
s[0][j] = s[0][j-1] - gap      # j letters of seq2 all aligned to gaps
```

**Recurrence.** For every interior cell, take the best of the three moves:

```
diag = s[i-1][j-1] + scorer(seq1[i-1], seq2[j-1])   # match/mismatch, reward added or penalty subtracted inside scorer
up   = s[i-1][j]   - gap                             # gap in seq2 (consume a seq1 letter)
left = s[i][j-1]   - gap                             # gap in seq1 (consume a seq2 letter)
s[i][j] = max(diag, up, left)
```

Note the sign convention the prompt asks about: the match reward is *added* through `scorer` (a positive number for equal DNA letters), while a mismatch and a gap are positive penalties that are *subtracted* (`scorer` returns `-mismatch`, and `gap` is subtracted directly).

**Traceback.** The final score is `s[n][m]`, the bottom-right corner, because a global alignment must use all of both sequences. The reconstruction is done by the shared helper `_traceback_linear(..., n, m, local=False)`. Starting at `(n, m)` it repeatedly asks "which of the three moves produced this cell's value?" by recomputing each candidate and comparing with a tolerance `eps = 1e-9`:

```
while i > 0 or j > 0:
    if s[i][j] == s[i-1][j-1] + scorer(seq1[i-1], seq2[j-1]):   # diagonal
        emit (seq1[i-1], seq2[j-1]); i -= 1; j -= 1
    elif s[i][j] == s[i-1][j] - gap:                            # up
        emit (seq1[i-1], "-"); i -= 1
    else:                                                       # left
        emit ("-", seq2[j-1]); j -= 1
```

The emitted columns are collected in reverse and then reversed at the end, giving `row1, row2`. `align_global` returns `(row1, row2, s[n][m])`.

**Complexity.** Time `O(n*m)` to fill the table. Space `O(n*m)` because the full table is kept for the traceback.

**Where it is used.** It is the default aligner in the app and the base case that local, semi-global, and banded alignment specialize.

---

## 2. Local alignment (Smith-Waterman): `align_local`

**What it does.** Finds the single best-scoring *substring-to-substring* match. Instead of forcing the whole sequences to line up, it hunts for the one region where they agree most. This is what you use to find a conserved domain buried inside two otherwise dissimilar sequences.

**The one change that makes it local: floor at 0.** The recurrence is identical to global, except a fourth option, `0.0`, is added to the `max`:

```
diag = s[i-1][j-1] + scorer(seq1[i-1], seq2[j-1])
up   = s[i-1][j]   - gap
left = s[i][j-1]   - gap
s[i][j] = max(0.0, diag, up, left)
```

The `0.0` means "an alignment can start fresh here." As soon as the running score would go negative, the cell resets to 0, which is exactly saying "throw away everything before this point and begin a new local alignment." Because of this floor, the entire first row and column are already 0, so no explicit initialization loop is needed (the table starts filled with `0.0`).

**Finding the end.** The best local alignment does not have to end at the corner, so the code tracks the maximum cell seen anywhere while filling the table:

```
if s[i][j] > best:
    best, bi, bj = s[i][j], i, j
```

**Traceback.** It starts at that best cell `(bi, bj)` and walks back using the same `_traceback_linear`, but with `local=True`. The extra condition `not (local and s[i][j] <= eps)` stops the walk the moment it reaches a cell whose value is 0. So the traceback runs "from the max cell back to a 0 cell," which precisely carves out the high-scoring region. The helper also returns the coordinates `(si, sj)` where it stopped. `align_local` returns `(row1, row2, best, si, bi, sj, bj)`, i.e. the aligned strings, the score, and the start/end indices in each sequence.

**Complexity.** Time `O(n*m)`, space `O(n*m)`, same as global.

**Where it is used.** The local-alignment mode of the app, for finding the best matching sub-region rather than an end-to-end alignment.

---

## 3. Semi-global alignment (free end gaps): `align_semiglobal`

**What it does.** A middle ground between global and local. The *interior* still has to align continuously, but gaps hanging off either end (overhangs) are **free**. This is the correct model when one sequence is expected to sit inside the other, or when two reads overlap: penalizing the overhang the way global alignment does would be unfair.

**Why end gaps are free (two separate tricks).**

*Free leading gaps* come from the initialization. Unlike global (where the first row and column ramp down by `-gap`), here the table is left as all zeros and never seeded:

```
s = [[0.0]*(m+1) for _ in range(n+1)]   # first row and column stay 0
```

A path can therefore reach any cell along the top or left edge for free, meaning a run of gaps at the *start* of either sequence costs nothing. The interior recurrence is the ordinary global one (`max(diag, up, left)`, no floor at 0).

*Free trailing gaps* come from where we start the traceback. In global we always begin at the corner `s[n][m]`; here we begin at the **best cell in the entire last row or last column**:

```
best, bi, bj = s[n][0], n, 0
for j in range(m+1):
    if s[n][j] >= best: best, bi, bj = s[n][j], n, j     # best in last row
for i in range(n+1):
    if s[i][m] > best:  best, bi, bj = s[i][m], i, m     # best in last column
```

Ending in the last row means all remaining letters of `seq2` become a trailing gap; ending in the last column means all remaining letters of `seq1` do. Either way that trailing overhang is never charged, because the score we report is `best`, the value at that chosen end cell.

**Traceback and stitching the overhangs.** From `(bi, bj)` it walks back with the standard linear rule until it hits the top or left edge (`while i > 0 and j > 0`), producing the aligned **core** `core1/core2` and the start coordinates `(si, sj)`. The unaligned overhangs are then pasted on explicitly as free end gaps so the returned strings still show the full sequences:

```
lead1 = seq1[:si] + "-" * sj          # seq1's leading letters, then gaps for seq2's lead
lead2 = "-" * si  + seq2[:sj]
tail1 = seq1[bi:] + "-" * (m - bj)     # seq1's trailing letters, then gaps for seq2's tail
tail2 = "-" * (n - bi) + seq2[bj:]
return lead1 + core + tail1,  lead2 + core + tail2,  best
```

So the two returned rows display the complete `seq1` and `seq2`, but the score counts only the aligned core; the leading and trailing overhangs contribute nothing.

**Complexity.** Time `O(n*m)` to fill the table plus `O(n + m)` to scan the last row and column; space `O(n*m)`.

**Where it is used.** The semi-global / "overlap" mode of the app, for cases where one sequence overlaps or is contained in the other.

---

## 4. Affine gaps (Gotoh): `align_global_affine`

**What it does.** A more biologically realistic gap model. Under the linear model every gap symbol costs the same `gap`, so a single long deletion is charged the same as many scattered short ones. Affine gaps instead charge `gap_open` to *start* a gap and `gap_extend` for each *additional* symbol in that same gap, so one long indel is cheaper than several short ones. This is the "gap-problem" stretch goal.

**Three matrices (states).** Because the cost of adding a gap symbol now depends on whether we are already inside a gap, one table is not enough. Gotoh's method uses three tables of shape `(n+1) x (m+1)`, one per "state" of the last move:

```
M[i][j]  = best score for seq1[:i] vs seq2[:j] ending in a match/mismatch (diagonal)
Ix[i][j] = best score ending in a gap in seq2  (a seq1 letter aligned to '-')
Iy[i][j] = best score ending in a gap in seq1  (a seq2 letter aligned to '-')
```

**Initialization.** `M[0][0] = 0`, everything else starts at `NEG_INF` (an unreachable score). The first column of `Ix` and first row of `Iy` are seeded with one open plus extends:

```
Ix[i][0] = -gap_open - (i-1)*gap_extend      # a leading run of i gaps in seq2
Iy[0][j] = -gap_open - (j-1)*gap_extend       # a leading run of j gaps in seq1
```

**Recurrences (including the Ix<->Iy transitions).**

```
sub = scorer(seq1[i-1], seq2[j-1])

M[i][j]  = max(M[i-1][j-1], Ix[i-1][j-1], Iy[i-1][j-1]) + sub

Ix[i][j] = max( M[i-1][j] - gap_open,      # open a new gap in seq2 from a match
                Ix[i-1][j] - gap_extend,   # extend an existing gap in seq2
                Iy[i-1][j] - gap_open )    # open a gap in seq2 right after a gap in seq1

Iy[i][j] = max( M[i][j-1] - gap_open,      # open a new gap in seq1 from a match
                Iy[i][j-1] - gap_extend,   # extend an existing gap in seq1
                Ix[i][j-1] - gap_open )    # open a gap in seq1 right after a gap in seq2
```

Two subtle points the prompt calls out:

- **The `Ix <-> Iy` transition.** A gap in one sequence can open *immediately after* a gap in the other. That is the third term in each recurrence (`Iy[i-1][j] - gap_open` inside `Ix`, and `Ix[i][j-1] - gap_open` inside `Iy`). Switching from a gap in `seq1` to a gap in `seq2` is a brand-new gap, so it pays `gap_open`, not `gap_extend`. Without these terms the aligner could not place two adjacent indels of opposite kind.
- **`M` can come from any of the three states.** Right after any kind of gap you can resume matching, so `M[i][j]` takes the best of `M`, `Ix`, `Iy` at the diagonal predecessor and adds `sub`.
- **Why open == extend reduces to linear.** If you set `gap_open == gap_extend == gap`, then opening and extending cost exactly the same, the distinction between "inside a gap" and "just starting one" disappears, and the three-state recurrence collapses to the single-table linear Needleman-Wunsch of section 1.

**Traceback across states.** The final score is the best of the three corner cells, `max(M[n][m], Ix[n][m], Iy[n][m])`. `_traceback_affine` first picks the starting **state** (the matrix holding that best end value), then walks back cell by cell, and at each step decides which state it *came from* by checking which predecessor equation is satisfied (again within `eps = 1e-9`):

- In state `M`: emit a real pair `(seq1[i-1], seq2[j-1])`, subtract `sub` to get `prev`, step diagonally (`i -= 1; j -= 1`), and set the new state to whichever of `M/Ix/Iy` at the new cell equals `prev`.
- In state `Ix`: emit `(seq1[i-1], "-")`, decide whether this gap was opened from `M`, opened from `Iy`, or extended from `Ix`, then step up (`i -= 1`).
- In state `Iy`: emit `("-", seq2[j-1])`, decide open-from-`M` / open-from-`Ix` / extend-from-`Iy`, then step left (`j -= 1`).

Columns are collected reversed and flipped at the end. `align_global_affine` returns `(row1, row2, max(M[n][m], Ix[n][m], Iy[n][m]))`.

**Complexity.** Time `O(n*m)` (a constant amount of work per cell across three tables). Space `O(n*m)` because all three full tables are retained for the traceback.

**Where it is used.** The affine-gap mode of the app, the gap-problem stretch goal, where a long single indel should be preferred over many short ones.

---

## 5. Banded alignment: `align_banded`

**What it does.** A memory- and time-saving version of global alignment for **long, similar** sequences. If the two sequences are nearly identical, the optimal path hugs the main diagonal of the DP table, so there is no point computing cells far off that diagonal. Banded alignment only fills a diagonal stripe of width `band` on either side.

**The band idea and its width.** A full table cell `(i, j)` lies on the main diagonal when `i == j`. The band keeps only cells with `|i - j| <= band`. If the caller does not pass a `band`, it is auto-sized to cover the length difference plus a cushion, then capped so it never exceeds the longer sequence:

```
if band is None:
    band = abs(n - m) + max(64, max(n, m) // 20)
band = min(band, max(n, m))
W = 2*band + 1          # number of columns actually stored per row
```

**Band-relative coordinates.** Instead of storing a full row of length `m+1`, each row stores only `W = 2*band + 1` cells. A real column index `j` maps to a band-relative index `c` by:

```
c = j - (i - band)      # equivalently  j = i - band + c
```

So within row `i`, band position `c = 0` corresponds to `j = i - band` (the far-left edge of the stripe) and `c = W-1` to `j = i + band`. Cells whose real `j` falls outside `0..m` are marked invalid.

**Only two rows in memory.** The score arrays are just `prev` and `cur`, each of length `W`, swapped after every row (`prev, cur = cur, prev`). That is why the memory for scores is `O(band)` per row rather than `O(m)`.

**Initialization (row 0 and column 0).** Row 0 is the pure-gap top edge: for each valid band column, `prev[c] = -j*gap` (direction "left", code `2`), and out-of-range columns get `NEG` with direction `3` (invalid). Inside the main loop, whenever a cell's real column is `j == 0` it is the pure-gap left edge, so `cur[c] = -i*gap` with direction "up", code `1`.

**Recurrence (in band-relative form).** For a valid interior cell (`si = seq1[i-1]`):

```
diag = prev[c] + scorer(si, seq2[j-1])          # same c: (i-1, j-1)
up   = prev[c+1] - gap   if c+1 < W  else NEG    # (i-1, j) shifts +1 in band coords
left = cur[c-1]  - gap   if c-1 >= 0 else NEG    # (i, j-1) shifts -1 in band coords
best = diag; d = 0
if up   > best: best = up;   d = 1
if left > best: best = left; d = 2
cur[c] = best; drow[c] = d
```

The reason `up` reads `prev[c+1]` (not `prev[c]`) is the coordinate shift: moving from row `i-1` to row `i` slides the band by one, so cell `(i-1, j)` sits at band position `c+1` in the previous row. Likewise `(i, j-1)` is at `c-1` in the current row. Off-band neighbors are treated as `NEG` so the path is forced to stay inside the stripe.

**Direction bytes for traceback.** Full score history is not kept, but the *direction* taken in each cell is. `dirs` is a list of `bytearray`s (one per row), each of length `W`, storing `0 = diag`, `1 = up`, `2 = left`, `3 = invalid`. This is cheap: one byte per band cell.

**Reading off the score.** The finish cell is the real corner `(n, m)`. Its band-relative column is `cend = m - (n - band)`, and the score is `prev[cend]` (after the loop, `prev` holds the last row) when that index is in range, else `NEG`.

**Traceback.** Starting at `(i, j) = (n, m)`, each step recomputes `c = j - (i - band)`, reads `d = dirs[i][c]`, and moves accordingly: `d == 0` diagonal, `d == 1` up (gap in seq2), `d == 2` left (gap in seq1). There are fallback branches for the case `d == 3` (a cell that was never validly computed, which can happen if the true path left the band): it then just takes whatever move is still possible (diagonal if both indices remain, else up, else left) so the traceback always terminates at `(0, 0)`. The columns are collected reversed and flipped. `align_banded` returns `(row1, row2, score)`.

**Why memory is `O(n*band)`, not `O(n*m)`.** The scores use only two rows of width `W = 2*band + 1`, i.e. `O(band)`. The direction history `dirs` keeps one `bytearray` of width `W` per row, `n+1` rows total, so `O(n*band)`. For similar sequences `band` is small and roughly constant, so this is far less than the `O(n*m)` of a full table. Time is `O(n*band)` as well, since each row only touches `W` cells.

**When it is exact.** Banded alignment gives the *same* answer as full Needleman-Wunsch **exactly when the optimal path stays inside the band** (`|i - j| <= band` for every cell on that path). Because off-band neighbors are `NEG`, any optimal path that would need to wander farther than `band` off the diagonal is forbidden, and the result may then be suboptimal. The auto-sizing (`abs(n - m) + max(64, max(n, m) // 20)`) is a heuristic that makes the band wide enough to be exact for typical similar sequences, but it is not a guarantee for arbitrary inputs.

**Complexity.** Time `O(n*band)`, score memory `O(band)`, direction memory `O(n*band)`.

**Where it is used.** The banded mode of the app, for aligning long and highly similar sequences where the full `O(n*m)` table would be wasteful.

---

**Source file:** `C:\Users\thomas\Desktop\Alignment\app\engine.py`

---

# Part 2. Supporting algorithms and data


This section documents the helper algorithms and the scoring data that sit around the core alignment routines. Two groups of code are covered:

1. The **original assignment files** in `Alignment/Alignment/` (`editDistanceTable.py`, `lcsLength.py`, `lcs.py`, `sharedKMers.py`, `manhattanTourist.py`). These are the classic dynamic programming (DP) building blocks from the course.
2. The **general engine** in `app/engine.py`, which adds DNA and protein scoring, the substitution matrices, DNA to protein translation, reverse complement, and the readout functions (`percent_identity`, `gap_stats`, `alignment_score`).

Throughout, `n` and `m` are the lengths of the two input sequences unless stated otherwise.

---

## 1. Edit distance table (`edit_distance_table`)

**File:** `Alignment/Alignment/editDistanceTable.py`

### What it does
Builds the full DP table for the **Levenshtein edit distance** between two strings `str1` and `str2`. Edit distance is the minimum number of single-character insertions, deletions, and substitutions needed to turn one string into the other. Entry `(i, j)` of the returned table is the edit distance between the first `i` characters of `str1` and the first `j` characters of `str2`, so the bottom-right corner is the distance between the whole strings.

### Setup and rules
The function raises `ValueError` if either string is empty. It then allocates a table `scoring_matrix` with `num_rows = len(str1) + 1` rows and `num_cols = len(str2) + 1` columns, filled with zeros.

The first row and first column are the base cases. Turning an empty prefix into a length-`j` prefix costs `j` insertions, and turning a length-`i` prefix into empty costs `i` deletions:

```
scoring_matrix[0][j] = j     for every column j
scoring_matrix[i][0] = i     for every row i
```

Every interior cell is the minimum over the three edit moves:

```
up   = scoring_matrix[row-1][col] + 1        # delete str1[row-1]
left = scoring_matrix[row][col-1] + 1        # insert str2[col-1]
diag = scoring_matrix[row-1][col-1]          # match/substitute
if str1[row-1] != str2[col-1]:
    diag += 1                                # substitution costs 1
scoring_matrix[row][col] = min(up, left, diag)
```

The key detail is that `diag` only adds 1 when the two characters differ. Equal characters are copied through for free (a free match), so the diagonal move is either a cost-0 match or a cost-1 substitution.

### Output
The function returns the entire table (a `list[list[int]]`), not just a single number. Callers read `scoring_matrix[len(str1)][len(str2)]` for the final distance. There is no traceback in this file; it only produces the score table.

### Complexity
- **Time:** `O(n * m)`, one constant-time `min` per cell.
- **Space:** `O(n * m)` because the full table is built and returned.

### Where it is used in this project
This is the course version of edit distance. In the app, `main.py` computes SARS-CoV vs SARS-CoV-2 edit distance through the related `edit_distance` module, and builds a hemoglobin pairwise distance matrix. The `edit_distance_table` function is the table-returning form of the same DP. It is conceptually the same recurrence the engine's global aligner uses, except edit distance minimizes a cost while the engine maximizes a score.

---

## 2. Longest common subsequence

**Files:** `Alignment/Alignment/lcsLength.py` (score matrix and length) and `Alignment/Alignment/lcs.py` (backtrack to recover the actual subsequence).

A **subsequence** keeps characters in order but allows skipping. The longest common subsequence (LCS) of two strings is the longest string that is a subsequence of both. Unlike a substring, the characters need not be contiguous.

### 2a. Score matrix (`lcs_score_matrix`) and length (`lcs_length`)

`lcs_score_matrix(str1, str2)` raises `ValueError` on an empty input, then builds a zero-filled table of `len(str1) + 1` rows and `len(str2) + 1` columns. The first row and column stay 0 (an empty prefix shares nothing). Each interior cell is:

```
up   = scoring_matrix[row-1][col]
left = scoring_matrix[row][col-1]
diag = scoring_matrix[row-1][col-1]
if str1[row-1] == str2[col-1]:
    diag += 1                              # a matching character extends the LCS
scoring_matrix[row][col] = max(up, left, diag)
```

Here the diagonal move rewards a **match** with `+1`; mismatches gain nothing. Because we take the `max`, the value never decreases as we move down or right, and `scoring_matrix[i][j]` equals the LCS length of the two prefixes.

`lcs_length(str1, str2)` simply calls `lcs_score_matrix` and returns the corner value `scoring_matrix[len(str1)][len(str2)]`.

### 2b. Backtrack to recover the subsequence (`longest_common_subsequence`)

**File:** `Alignment/Alignment/lcs.py`. It imports `lcs_score_matrix` from `lcsLength`.

For empty input it returns `""`. Otherwise it rebuilds the score matrix, sets `r = len(str1)`, `c = len(str2)`, and walks from the bottom-right corner back toward the origin, prepending matched characters as it goes:

```
while r != 0 and c != 0:
    up   = matrix[r-1][c]
    left = matrix[r][c-1]
    diag = matrix[r-1][c-1]
    if str1[r-1] == str2[c-1]:
        diag += 1

    if matrix[r][c] == up:          # came from above -> move up
        r -= 1
    elif matrix[r][c] == left:      # came from the left -> move left
        c -= 1
    elif matrix[r][c] == diag:      # came from the diagonal
        if str1[r-1] == str2[c-1]:
            lcs = str1[r-1] + lcs   # only prepend on a real match
        r -= 1; c -= 1
    else:
        raise ValueError("Error: bad scoring matrix state.")
return lcs
```

How the output is built: the character is prepended (`lcs = str1[r-1] + lcs`) rather than appended, so the subsequence comes out in forward order without a final reversal. A character is only added on a diagonal step **and** only when `str1[r-1] == str2[c-1]`. The `if/elif` order checks `up` first, so when several moves tie in value the backtrack prefers up, then left, then diagonal, which is a valid (though implementation-specific) choice of one optimal LCS.

### Complexity
- **Time:** `O(n * m)` to fill the matrix; the backtrack adds `O(n + m)` steps, so total is `O(n * m)`.
- **Space:** `O(n * m)` for the score matrix.

### Where it is used in this project
`lcs_length` is called from `main.py` to report the LCS length between the two SARS coronavirus genomes as a quick similarity measure. The `longest_common_subsequence` backtrack (in `lcs.py`) is the course exercise that shows how to turn a score matrix into an actual answer, and it is the same idea the engine's alignment tracebacks use, except the engine's tracebacks also emit gap symbols and use scores instead of pure match counts.

---

## 3. Count shared k-mers (`count_shared_kmers`)

**File:** `Alignment/Alignment/sharedKMers.py`

### What it does
A **k-mer** is a substring of length `k`. This counts how many k-mers the two strings share, counting multiplicity. If a k-mer appears 3 times in one string and 2 times in the other, it contributes `min(3, 2) = 2` to the total. This is an alignment-free similarity measure: it does not build a DP table, it just compares the multisets of length-`k` windows.

### Rule
```
if k <= 0 or k > len(str1) or k > len(str2):
    return 0

count every length-k window of str1 into counts1   # dict k-mer -> count
count every length-k window of str2 into counts2

shared = 0
for each kmer in counts1:
    if kmer in counts2:
        shared += min(counts1[kmer], counts2[kmer])
return shared
```

A string of length `L` has `L - k + 1` overlapping windows, generated by the loop `for i in range(len(str) - k + 1)`. The guard at the top returns `0` for a nonpositive `k` or a `k` larger than either string, avoiding empty or negative-range loops.

### Output
A single integer: the total number of shared k-mer occurrences (the summed minima).

### Complexity
Let `L1 = len(str1)`, `L2 = len(str2)`.
- **Time:** `O(L1 * k + L2 * k)` because building each k-mer slice costs `O(k)`, and there are about `L1` and `L2` of them. The final comparison loop is `O(number of distinct k-mers in str1)`, each a `O(k)` dictionary lookup, which is dominated by the counting cost.
- **Space:** `O(L1 + L2)` k-mer entries across the two dictionaries (fewer if there are repeats).

### Where it is used in this project
It is the standalone course "shared k-mers" exercise. In the app it appears in `app/test_all.py`, which imports `count_shared_kmers` and checks that `count_shared_kmers("GATTACA", "GATTACA", 3) == 5` (the 5 length-3 windows of a 7-character string, each shared with itself). It is a fast, coarse similarity check that is independent of the alignment machinery.

---

## 4. Manhattan tourist (`manhattan_tourist`)

**File:** `Alignment/Alignment/manhattanTourist.py`

### What it does
Solves the **Manhattan Tourist Problem**: on an `n` by `m` grid of intersections where every edge carries a weight, find the maximum total weight of any path from the source `(0, 0)` to the sink `(n, m)` when you may only move **down** or **right**. This is the grid-graph longest-path DP that underlies alignment. Alignment is exactly this problem on a specially weighted grid where diagonal edges are also allowed.

### Inputs
- `n`, `m`: grid dimensions (number of rows and columns of edges).
- `down`: a matrix where `down[i-1][j]` is the weight of the edge going downward into cell `(i, j)`.
- `right`: a matrix where `right[i][j-1]` is the weight of the edge going rightward into cell `(i, j)`.

### Rule
Allocate an `(n+1)` by `(m+1)` table `s` of zeros. The first column can only be reached by going straight down, and the first row only by going straight right, so those are running sums:

```
s[i][0] = s[i-1][0] + down[i-1][0]      # first column: only downward edges
s[0][j] = s[0][j-1] + right[0][j-1]     # first row: only rightward edges
```

Every interior cell takes the better of arriving from above or from the left:

```
from_above = s[i-1][j] + down[i-1][j]
from_left  = s[i][j-1] + right[i][j-1]
s[i][j] = max(from_above, from_left)
```

### Output
Returns the single integer `s[n][m]`, the weight of a longest path to the sink. This version returns only the best weight; it does not reconstruct the path.

### Complexity
- **Time:** `O(n * m)`, one `max` per cell.
- **Space:** `O(n * m)` for the table.

### Where it is used in this project
It is the course grid-DP exercise. Its role in the project is conceptual: it is the same "best value into each cell from a small set of predecessors" pattern that the engine's aligners use. In `engine.py`, `align_global` computes each cell as `max(diag, up, left)`, which is the Manhattan tourist recurrence plus a diagonal (match/mismatch) edge and gap penalties on the down and right edges.

---

## 5. Scoring functions (`app/engine.py`)

A **scorer** in the engine is just a function `score(a, b) -> number` that returns how good it is to line up character `a` against character `b`. Passing the scorer as an argument is what lets the same alignment code handle both DNA and protein.

### 5a. DNA scorer (`dna_scorer`)
```python
def dna_scorer(match, mismatch):
    def score(a, b):
        return match if a == b else -mismatch
    return score
```
`dna_scorer(match, mismatch)` is a factory: it returns a closure that gives `+match` when the two bases are equal and `-mismatch` when they differ. Note that `mismatch` is passed as a positive number and negated inside, so `dna_scorer(1, 1)` rewards a match with `+1` and penalizes a mismatch with `-1`. Each call is `O(1)`.

### 5b. Matrix scorer (`matrix_scorer`)
```python
def matrix_scorer(matrix):
    def score(a, b):
        if (a, b) in matrix:
            return matrix[(a, b)]
        a = a if (a, "A") in matrix else "X"
        b = b if ("A", b) in matrix else "X"
        return matrix[(a, b)]
    return score
```
`matrix_scorer(matrix)` returns a closure that looks the pair `(a, b)` up in a substitution matrix dictionary (see BLOSUM62 and PAM250 below). If the exact pair is present it is returned directly. Otherwise the fallback logic replaces any letter not found in the matrix with `"X"`, the "any residue" symbol, so unknown or ambiguous characters still score against the matrix's `X` row and column instead of crashing. The membership tests `(a, "A")` and `("A", b)` check whether the individual letters exist as row and column keys. Each lookup is `O(1)`.

Both scorers are the same interface, so `align_global`, `align_local`, and the affine and banded variants can take either one unchanged.

---

## 6. Substitution matrices: BLOSUM62 and PAM250

**File:** `app/engine.py`, module-level text blocks `_BLOSUM62_TEXT` and `_PAM250_TEXT`, parsed by `_parse_matrix`.

### What they mean
A **substitution matrix** scores how likely one amino acid is to be replaced by another over evolutionary time. Both matrices are `20 x 20` amino acids plus the ambiguity/special symbols `B` (Asx), `Z` (Glx), `X` (any), and `*` (stop or gap). The values are the exact standard NCBI/Biopython numbers.

- **BLOSUM62** (BLOcks SUbstitution Matrix, 62% identity clustering) is the common default for protein alignment. Diagonal (identical residue) entries are positive, for example `W/W = 11`, `C/C = 9`, `H/H = 8`. Conservative swaps that preserve chemistry are near zero or mildly positive, and unlikely swaps are negative, for example `W/D = -4`.
- **PAM250** (Point Accepted Mutation, 250 mutations per 100 residues) reflects a larger evolutionary distance, so it is more tolerant of substitutions. Its diagonal values are larger and more spread out, for example `W/W = 17`, `C/C = 12`, and off-diagonal penalties are steeper, for example `*` rows are `-8`.

**Positive scores** mean the pairing is favorable (identical or chemically similar residues, so the alignment should keep them together). **Negative scores** mean the pairing is unfavorable (chemically dissimilar residues that rarely substitute), which pushes the aligner toward a gap or a different alignment instead.

### How they are parsed (`_parse_matrix`)
```python
def _parse_matrix(text):
    lines = [ln for ln in text.strip().splitlines() if ln.strip()]
    columns = lines[0].split()                      # header letters, in order
    scores = {}
    for ln in lines[1:]:
        parts = ln.split()
        row_letter = parts[0]                       # this row's amino acid
        for col_letter, value in zip(columns, parts[1:]):
            scores[(row_letter, col_letter)] = int(value)
    return scores
```
The parser drops blank lines, reads the first non-blank line as the ordered list of column letters, then for each remaining line takes `parts[0]` as the row letter and pairs each remaining number with its column letter via `zip`. Every cell becomes a dictionary entry `{(row_letter, col_letter): int(value)}`. The result is a flat `dict[tuple[str, str], int]` that `matrix_scorer` can index in `O(1)`.

The module runs the parser once at import time:
```python
BLOSUM62 = _parse_matrix(_BLOSUM62_TEXT)
PAM250   = _parse_matrix(_PAM250_TEXT)
```

### Complexity
Parsing is `O(R * C)` where `R` and `C` are the numbers of rows and columns (24 each here, so effectively constant). It runs exactly once per matrix at import. Each subsequent scorer lookup is `O(1)`. Storing the matrix is `O(R * C)` dictionary entries.

### Where they are used in this project
They are the protein scoring option. The docstring at the top of `engine.py` names them as the protein path alongside DNA match/mismatch scoring. `make_figures.py` uses BLOSUM62 for the hemoglobin protein alignment figures (for example "good (BLOSUM62, gap = 10)" versus "bad (BLOSUM62, gap = 0)"), passing the matrix through `matrix_scorer` into the aligners.

---

## 7. DNA to protein translation (`translate`) and the codon table

**File:** `app/engine.py`, dictionary `_CODON_TABLE` and function `translate`.

### The codon table
`_CODON_TABLE` is the **standard genetic code**: a dictionary mapping each of the 64 three-letter DNA codons to its one-letter amino acid, using `T` (not `U`). The three stop codons `TAA`, `TAG`, and `TGA` map to `"*"`. For example `ATG -> M` (Methionine, the start codon), `TTT -> F`, `GGG -> G`.

### Reading frame and rule (`translate`)
```python
def translate(dna, stop_at_stop=True):
    dna = dna.upper().replace("U", "T")     # normalize: uppercase, RNA U -> DNA T
    protein = []
    for i in range(0, len(dna) - 2, 3):     # step by 3: reading frame 0
        aa = _CODON_TABLE.get(dna[i:i+3], "X")
        if aa == "*":
            if stop_at_stop:
                break                       # stop translation at a stop codon
            aa = "*"
        protein.append(aa)
    return "".join(protein)
```

Key details:
- **Reading frame:** the loop starts at index 0 and steps by 3, so it reads frame 0 only (codons `dna[0:3]`, `dna[3:6]`, and so on). It does not try all three frames or the reverse strand.
- **Range bound:** `range(0, len(dna) - 2, 3)` guarantees every slice `dna[i:i+3]` has all three bases, so a trailing 1 or 2 leftover bases are ignored rather than producing a short codon.
- **Normalization:** input is uppercased and `U` is replaced by `T`, so RNA sequences translate correctly.
- **Unknown codons:** `.get(codon, "X")` returns `"X"` for any codon not in the table (for example one containing `N`), so ambiguous DNA yields `X` rather than an error.
- **Stop handling:** by default (`stop_at_stop=True`) translation halts at the first stop codon and the `*` is not appended. With `stop_at_stop=False` the `*` is kept and translation continues.

### Output
A protein string of one-letter amino acid codes.

### Complexity
- **Time:** `O(L)` for a DNA string of length `L` (one dictionary lookup and one 3-character slice per codon, `L/3` codons).
- **Space:** `O(L)` for the output protein (about `L/3` characters).

### Where it is used in this project
`translate` powers the "DNA vs protein" comparison mentioned in the engine docstring. `make_3d.py` and `make_figures.py` translate the SARS-CoV and SARS-CoV-2 spike **DNA** coding sequences into protein before aligning them, for example `p1 = e.translate(e.read_fasta(".../SARS-CoV_genome_spike_protein.fasta"))`. The test suite verifies it on small cases (`translate("ATGGCCTGA") == "MA"`, where `TGA` is a stop that halts output) and confirms it reproduces the provided spike protein FASTA files exactly.

---

## 8. Reverse complement (`reverse_complement`)

**File:** `app/engine.py`

### What it does
Returns the **reverse complement** of a DNA string: the sequence you would read on the opposite strand in the 5' to 3' direction. Each base is swapped for its Watson-Crick partner (`A<->T`, `C<->G`) and the order is reversed.

### Rule
```python
comp = {"A":"T","T":"A","C":"G","G":"C","N":"N",
        "a":"t","t":"a","c":"g","g":"c","n":"n","U":"A","u":"a"}
return "".join(comp.get(c, c) for c in reversed(dna))
```
The dictionary handles both cases, treats `N` as its own complement (unknown stays unknown), and maps `U`/`u` (RNA uracil) to `A`/`a`. Any character not in the map is passed through unchanged via `comp.get(c, c)`. Iterating over `reversed(dna)` while complementing does both operations in a single pass.

### Complexity
- **Time:** `O(L)` for length `L`, one dictionary lookup per base.
- **Space:** `O(L)` for the output string.

### Where it is used in this project
It is available for comparing a sequence against the opposite strand (useful because a gene can lie on either strand). `test_all.py` checks a concrete case (`reverse_complement("ATTTAC") == "GTAAAT"`) and verifies the **involution** property: applying it twice returns the original sequence for random DNA.

---

## 9. Readouts: percent identity, gap statistics, and score recomputation

These three functions in `app/engine.py` summarize an already-computed alignment. Each takes two aligned rows (`row1`, `row2`) that include gap characters `"-"`.

### 9a. Percent identity (`percent_identity`)
```python
def percent_identity(row1, row2):
    matches = aligned = 0
    for a, b in zip(row1, row2):
        if a == "-" or b == "-":
            continue                 # skip gap columns entirely
        aligned += 1
        if a == b:
            matches += 1
    return 100.0 * matches / aligned if aligned else 0.0
```
It walks the two rows column by column, **ignoring any column that contains a gap**, counts how many of the remaining aligned columns are exact matches, and returns the percentage. If there are no non-gap columns it returns `0.0` (avoiding division by zero). Test example: `percent_identity("AC-GT", "ACTGT")` is `100.0` because the one gap column is excluded, leaving 4 of 4 columns matching.
- **Time:** `O(len(row))`. **Space:** `O(1)`.

### 9b. Gap statistics (`gap_stats`)
```python
def gap_stats(row1, row2):
    lengths = []
    run = 0
    for a, b in zip(row1, row2):
        if a == "-" or b == "-":
            run += 1                 # extend the current gap run
        else:
            if run:
                lengths.append(run) # a gap run just ended
            run = 0
    if run:
        lengths.append(run)         # flush a trailing gap
    return len(lengths), lengths
```
It counts **contiguous runs** of gap columns. A column is "gap" if either row has `"-"` there. It tracks the current run length, records it when a non-gap column ends the run, and flushes a final run at the end. It returns a tuple `(number of separate gaps, list of each gap's length)`. Test example: `gap_stats("A--GT-", "ABCGTX") == (2, [2, 1])`, meaning two separate gap runs of lengths 2 and 1.
- **Time:** `O(len(row))`. **Space:** `O(g)` where `g` is the number of gap runs.

This distinction matters for the affine-gap part of the project: one long insertion should count as **one** gap, not many. `make_figures.py` uses `gap_stats(...)[0]` (the gap count) to contrast a linear-gap alignment that "shatters" an insertion into many gaps against an affine-gap alignment that keeps it as one.

### 9c. Alignment score recomputation (`alignment_score`)
```python
def alignment_score(row1, row2, scorer, gap):
    total = 0.0
    for a, b in zip(row1, row2):
        if a == "-" or b == "-":
            total -= gap             # linear gap penalty per gap symbol
        else:
            total += scorer(a, b)    # match/mismatch or matrix score
    return total
```
Given a finished alignment, a scorer, and a **linear** gap penalty, it recomputes the score from scratch: each gap column subtracts `gap`, each aligned column adds `scorer(a, b)`. This models a flat penalty per gap symbol (it does **not** model affine open/extend costs), so it is used as a sanity check that the DP aligners returned a score consistent with their output. In `test_all.py` it verifies, for example, that `align_global`'s returned score matches `alignment_score(r[0], r[1], dna_scorer(m, mm), gp)` to within `1e-6`.
- **Time:** `O(len(row))` (each `scorer` call is `O(1)`). **Space:** `O(1)`.

---

## Summary table

| Function | File | Time | Space | Returns |
| --- | --- | --- | --- | --- |
| `edit_distance_table` | `editDistanceTable.py` | `O(n*m)` | `O(n*m)` | full DP cost table |
| `lcs_score_matrix` / `lcs_length` | `lcsLength.py` | `O(n*m)` | `O(n*m)` | LCS table / LCS length |
| `longest_common_subsequence` | `lcs.py` | `O(n*m)` | `O(n*m)` | an LCS string |
| `count_shared_kmers` | `sharedKMers.py` | `O((L1+L2)*k)` | `O(L1+L2)` | count of shared k-mers |
| `manhattan_tourist` | `manhattanTourist.py` | `O(n*m)` | `O(n*m)` | max path weight |
| `dna_scorer` / `matrix_scorer` | `engine.py` | `O(1)` per call | `O(1)` | a `score(a,b)` closure |
| `_parse_matrix` (BLOSUM62, PAM250) | `engine.py` | `O(R*C)` once | `O(R*C)` | `{(a,b): score}` dict |
| `translate` | `engine.py` | `O(L)` | `O(L)` | protein string |
| `reverse_complement` | `engine.py` | `O(L)` | `O(L)` | reverse-complement DNA |
| `percent_identity` | `engine.py` | `O(len)` | `O(1)` | identity percent |
| `gap_stats` | `engine.py` | `O(len)` | `O(g)` | `(count, [lengths])` |
| `alignment_score` | `engine.py` | `O(len)` | `O(1)` | recomputed linear score |

Absolute file paths referenced:
- `C:/Users/thomas/Desktop/Alignment/app/engine.py`
- `C:/Users/thomas/Desktop/Alignment/Alignment/editDistanceTable.py`
- `C:/Users/thomas/Desktop/Alignment/Alignment/lcsLength.py`
- `C:/Users/thomas/Desktop/Alignment/Alignment/lcs.py`
- `C:/Users/thomas/Desktop/Alignment/Alignment/sharedKMers.py`
- `C:/Users/thomas/Desktop/Alignment/Alignment/manhattanTourist.py`
- `C:/Users/thomas/Desktop/Alignment/Alignment/main.py` (usage of `lcs_length` / edit distance)
- `C:/Users/thomas/Desktop/Alignment/app/make_figures.py`, `make_3d.py`, `test_all.py` (usage of engine readouts, translate, matrices)

---

# Part 3. GPU accelerated alignment


Source file: `C:/Users/thomas/Desktop/Alignment/app/gpu_align.py`

This module computes global (Needleman-Wunsch) alignment **scores** for one query sequence against a large batch of equal-length target sequences, all at once, on whatever GPU is available. It is a scores-only accelerator: it never builds a traceback and never returns an alignment, only the final numeric score for each pair. The normal aligner elsewhere in the project handles one pair at a time on the CPU and does produce full alignments. This module exists for tasks that need the same query compared against thousands of sequences, where only the score matters.

## Why one alignment is sequential but many are parallel

A single pairwise dynamic-programming (DP) alignment fills a table `S` where the value in cell `(i, j)` depends on three neighbors that must already be filled:

```
S[i][j] = max(
    S[i-1][j-1] + substitution(query[i], target[j]),   # diagonal (match/mismatch)
    S[i-1][j]   - gap,                                  # up   (gap in target)
    S[i][j-1]   - gap,                                  # left (gap in query)
)
```

Because each cell reads its up, left, and diagonal neighbors, you cannot compute `(i, j)` until those three are done. That data dependency makes filling **one** table an inherently sequential, step-by-step process. There is no way to correctly compute the whole table in a single parallel shot, because later cells literally need earlier cells' results.

However, **different pairwise alignments are completely independent of each other.** Aligning the query to target #1 has nothing to do with aligning it to target #2. So if you have thousands of independent alignments to run, you can run all of them side by side. That is the parallelism this module exploits: not parallelism *within* one alignment's cell-to-cell dependency, but parallelism *across* many alignments (the batch dimension `B`) and, on top of that, across the independent cells of a single anti-diagonal (explained next). PyTorch tensors let a GPU do the same arithmetic on all `B` alignments simultaneously.

The module's docstring states the motivating use case directly: some tasks need the same query aligned against many other sequences, for example the significance test, which aligns against hundreds of shuffled sequences. Those are all independent, so they run in parallel on the GPU.

## The anti-diagonal wavefront

Even inside a single alignment table, not every cell depends on every other cell. Look again at the dependencies of cell `(i, j)`: its inputs are `(i-1, j-1)`, `(i-1, j)`, and `(i, j-1)`. Every one of those neighbors has a smaller value of `i + j` than `(i, j)` does. If we group cells by the sum `d = i + j` (called an **anti-diagonal**), then:

- All cells on anti-diagonal `d` depend only on cells with sum `d-1` and `d-2`.
- Two cells on the *same* anti-diagonal `d` never depend on each other, because neither is a neighbor of the other.

This means all cells on one anti-diagonal are mutually independent and can be computed at the same time. So instead of filling the table cell by cell (which would be `n * m` sequential steps for an `n` by `m` table), we fill it **one anti-diagonal at a time**. Each anti-diagonal becomes a single vectorized GPU operation that updates a whole strip of cells at once. This diagonal-by-diagonal sweep is called a **wavefront**: the computation moves across the table like a wave, and there are only `n + m - 1` waves total instead of `n * m` individual cell updates.

Combine the two levels of parallelism and each wavefront step is doing arithmetic on `B` alignments times `L` cells (where `L` is the length of that anti-diagonal) all in one tensor operation.

## How `batched_global_scores` and `_run_batch` work

The public entry point is:

```python
def batched_global_scores(query, targets, matrix, gap, alphabet, device=None):
```

- `query` is the single sequence (a string).
- `targets` is a list of sequences that **must all be the same length** `m` (the code reads `m = len(targets[0])` and assumes the rest match).
- `matrix` is a substitution-score dictionary of the form `{(a, b): score}` (for example BLOSUM62 or PAM250 entries).
- `gap` is the linear gap penalty (a positive number; the code subtracts it).
- `alphabet` is the set of allowed characters, and its order defines the index of each character.
- `device` is optional; if `None`, it is chosen automatically by `best_device()`.

It returns a plain Python list of scores, one per target. It wraps the real worker `_run_batch` in a `try/except` for the CPU fallback (described later), then delegates all real work to `_run_batch`.

### Encoding sequences and the substitution matrix as tensors

Before any DP runs, `_run_batch` converts everything into integer indices and float tensors so the GPU can index and add efficiently.

**The substitution matrix becomes an `[A, A]` lookup tensor.** `build_matrix_tensor(matrix, alphabet)` builds:

```python
idx = {c: i for i, c in enumerate(alphabet)}   # char -> integer index
A = len(alphabet)
M = torch.zeros((A, A), dtype=torch.float32)
for a in alphabet:
    for b in alphabet:
        M[idx[a], idx[b]] = matrix.get((a, b), matrix.get((a, "X"), 0))
```

So `M[i, j]` is the score of aligning the character with index `i` against the character with index `j`. If a specific `(a, b)` pair is missing from the dictionary, it falls back to the `(a, "X")` entry (the score against the "any/unknown" symbol), and if that is also missing, it falls back to `0`. This tensor is moved to the device as `SM = M.to(device)`.

**Sequences become integer-index tensors.** `encode(seq, idx)` maps each character to its integer index, and any character not in the alphabet is mapped to the index of `"X"` (or `0` if `"X"` is not present):

```python
def encode(seq, idx):
    x = idx.get("X", 0)
    return [idx.get(c, x) for c in seq]
```

In `_run_batch`:

```python
qi = torch.tensor(encode(query, idx), device=device)                 # shape [n]
ti = torch.tensor([encode(t, idx) for t in targets], device=device)  # shape [B, m]
```

where `n = len(query)`, `m = len(targets[0])`, and `B = len(targets)`. `qi` holds the query's indices; `ti` holds all `B` targets' indices stacked into one 2-D tensor.

### The score tensor and its initialization

The DP tables for all `B` alignments live in one 3-D tensor:

```python
S = torch.empty((B, n + 1, m + 1), device=device, dtype=torch.float32)
```

`S[b, i, j]` is the best score for alignment `b` using the first `i` characters of the query and the first `j` characters of target `b`. The `+1` in each dimension is the standard extra row and column for the empty-prefix case.

The boundaries are the pure-gap rows and columns, filled for **all** batches at once using broadcasting:

```python
# first row of every table: aligning i=0 query chars against j target chars = j gaps
S[:, 0, :] = -gap * torch.arange(m + 1, ...)
# first column of every table: aligning i query chars against j=0 target chars = i gaps
S[:, :, 0] = -gap * torch.arange(n + 1, ...).unsqueeze(0)
```

`S[:, 0, j] = -gap * j` and `S[:, i, 0] = -gap * i`, exactly the Needleman-Wunsch boundary conditions, applied identically to every alignment in the batch.

### Updating one anti-diagonal per step

The main loop walks the anti-diagonal index `d` from `2` up to `n + m` (inclusive). For each `d` it works out which cells `(i, j)` with `i + j = d` are actually inside the table:

```python
for d in range(2, n + m + 1):
    i_lo = max(1, d - m)      # smallest valid row on this diagonal
    i_hi = min(n, d - 1)      # largest valid row on this diagonal
    if i_lo > i_hi:
        continue              # this diagonal has no interior cells
    i_idx = torch.arange(i_lo, i_hi + 1, device=device)   # [L] the row numbers
    j_idx = d - i_idx                                      # [L] matching columns
```

`i_idx` and `j_idx` are the coordinate lists of every interior cell on this one anti-diagonal (length `L`), guaranteed to stay within `1..n` for rows and `1..m` for columns by the `i_lo`/`i_hi` clamps.

Then it gathers the substitution scores and the three neighbor values, and takes the max, for **all `B` batches and all `L` cells in one shot**:

```python
qrow = qi[i_idx - 1]                            # [L]   query char index for each cell
tcol = ti[:, j_idx - 1]                         # [B, L] target char index per batch/cell
sub  = SM[qrow.unsqueeze(0).expand(B, -1), tcol]# [B, L] substitution score per batch/cell

diag = S[:, i_idx - 1, j_idx - 1] + sub         # [B, L] diagonal (match/mismatch)
up   = S[:, i_idx - 1, j_idx]     - gap         # [B, L] gap in target
left = S[:, i_idx,     j_idx - 1] - gap         # [B, L] gap in query

S[:, i_idx, j_idx] = torch.maximum(torch.maximum(diag, up), left)
```

A few details worth noting:

- The `- 1` on `i_idx` and `j_idx` when indexing `qi` and `ti` converts from the DP table's 1-based prefix coordinates to the 0-based sequence positions.
- `sub` is a **2-D gather** into the `[A, A]` matrix tensor: `qrow` (broadcast to `[B, L]` via `unsqueeze(0).expand(B, -1)`) supplies the row indices and `tcol` supplies the column indices, so `sub[b, k]` is the score of aligning the query character at cell `k` against batch `b`'s target character at cell `k`.
- `diag`, `up`, and `left` reproduce exactly the three-way Needleman-Wunsch recurrence, but as whole-tensor operations. Because all cells on anti-diagonal `d` only read cells on diagonals `d-1` and `d-2` (already written on earlier iterations), writing them all at once is correct.

After the loop, `sync(device)` waits for the GPU to finish (it calls `torch.xpu.synchronize()` or `torch.cuda.synchronize()` on those backends, and does nothing otherwise), and the final scores are read out of the bottom-right corner of every table:

```python
return S[:, n, m].tolist()
```

`S[b, n, m]` is the completed global-alignment score for the query against target `b`, and `.tolist()` moves the results back to the CPU as an ordinary Python list.

## Device selection and the CPU fallback

`best_device()` chooses the best available accelerator by trying backends in a fixed order, and each check is wrapped in its own `try/except` so that a backend that is not installed or errors on the availability check is simply skipped rather than crashing the program:

```
1. Intel GPU        -> "xpu"   (torch.xpu.is_available())
2. NVIDIA / AMD     -> "cuda"  (torch.cuda.is_available(); AMD works via a ROCm build of torch)
3. Apple Silicon    -> "mps"   (torch.backends.mps.is_available())
4. none of the above-> "cpu"
```

So the priority is: Intel `xpu` first, then `cuda` (which covers both NVIDIA GPUs and AMD GPUs through a ROCm build of PyTorch), then Apple `mps`, and finally plain `cpu` if nothing else is available.

There is a second layer of safety in `batched_global_scores`. Even after a device is chosen, the actual GPU run can still fail (out of memory, an unsupported operation, a driver problem, etc.). The `try/except` handles that by retrying the whole batch on the CPU:

```python
try:
    return _run_batch(query, targets, matrix, gap, alphabet, device)
except Exception:
    if device != "cpu":
        return _run_batch(query, targets, matrix, gap, alphabet, "cpu")
    raise
```

If the run was already on the CPU, there is nothing to fall back to, so the exception is re-raised. Otherwise it transparently redoes the computation on the CPU and returns the same result. The net effect is that this code always produces an answer on any machine, GPU or not.

## Measured result

On an Intel Arc GPU (the `xpu` backend), the batched wavefront runs about **28x faster than the CPU** for this workload, and it produces the **exact same scores** as the CPU path. The scores are identical because the arithmetic is the same integer-derived recurrence; the GPU only changes how many cells are evaluated at once, not what each cell computes.

## Complexity and the scores-only limitation

Let `n` be the query length, `m` the (shared) target length, and `B` the number of targets in the batch.

- **Time complexity:** `O(B * n * m)` total arithmetic, the same total work as running `B` separate `n` by `m` DP fills. The win is not fewer operations but parallelism: there are only `n + m - 1` wavefront steps (loop iterations), and each step does `O(B * L)` work in parallel on the GPU, where `L` is that diagonal's length.
- **Space complexity:** `O(B * n * m)` for the full score tensor `S`, because this implementation keeps every table entry for every batch element. (A scores-only aligner could in principle keep just the last two anti-diagonals, but this code stores the whole `S` tensor.) The substitution tensor adds `O(A^2)` and the encoded sequences add `O(n + B*m)`, both small next to `S`.

**Scores only, no traceback.** This module computes the final alignment score for each pair and nothing else. Because it never stores traceback pointers and reads out only `S[:, n, m]`, it cannot reconstruct which characters aligned to which, and cannot show gaps or the aligned strings. That is by design: the tasks that use it (such as the significance test aligning against many shuffled sequences) only need the numeric score distribution, so skipping traceback saves both time and memory. Any task that needs an actual alignment must use the project's normal CPU aligner instead.

## Constants defined at the bottom

Two alphabet constants are provided for callers to pass as `alphabet`:

```python
PROTEIN_ALPHABET = "ARNDCQEGHILKMFPSTWYVBZX*"   # matches the BLOSUM62 / PAM250 matrices
DNA_ALPHABET     = "ACGTN"
```

The protein alphabet's character set and ordering line up with the standard BLOSUM62 and PAM250 scoring matrices (including the `B`, `Z`, `X` ambiguity codes and the `*` stop symbol), and the DNA alphabet includes `N` as the ambiguous-base symbol.

---

# Part 4. The on-device AI models


This project ships two neural networks that are trained offline in Python (PyTorch), exported to the ONNX format, and then run inside the WinUI app through Windows ML. Neither model is downloaded pre-trained. Both are trained from scratch on real biological data. This section documents exactly what each script does.

The two build scripts live in `C:\Users\thomas\Desktop\Alignment\app\`:

- `build_ss_model.py` builds the secondary-structure predictor and exports `ss_model.onnx`.
- `build_plm.py` trains the protein language model and saves the PyTorch weights `plm.pt`.
- `export_plm.py` loads `plm.pt`, wraps it as an embedding model, and exports `plm_embed.onnx` plus a precomputed family database `plm_db.tsv`.

---

## Model 1: Secondary-structure predictor (1D-CNN)

### What it does

Given a protein's amino-acid sequence (a string over the 20-letter alphabet `ACDEFGHIKLMNPQRSTVWY`), the model predicts, for every residue, one of three secondary-structure classes:

- helix (alpha helix),
- sheet (beta strand),
- coil (everything else, i.e. loops and turns).

This 3-class problem is the classic "Q3" secondary-structure prediction task. Predicting structure from sequence alone is a genuine (if small) machine-learning problem, because the same amino acid can be in a helix in one protein and a sheet in another. The model has to use the local sequence context to decide.

### Training data (this is real, not synthetic)

The training data is downloaded live from the RCSB Protein Data Bank (PDB), which holds experimentally determined 3D protein structures.

1. **Choosing which structures to download.** `get_pdb_ids(rows=400)` sends a JSON query to the RCSB search API (`https://search.rcsb.org/rcsbsearch/v2/query`). The query asks for entries that are:
   - `"Protein (only)"` (no DNA/RNA, no complexes),
   - between 60 and 250 residues long (`entity_poly.rcsb_sample_sequence_length` range 60 to 250),
   - solved at better than 2.0 Angstrom resolution (`rcsb_entry_info.resolution_combined` less than `2.0`), which means high-quality structures.

   It returns up to 400 candidate PDB identifiers.

2. **Reading sequence and structure off the coordinates.** For each ID, `get_seq_ss(pid)` downloads the raw `.pdb` coordinate file from `https://files.rcsb.org/download/<id>.pdb` (caching it under the `ss_cache` folder so re-runs are fast). Then, using the `biotite` library:
   - `pdb.PDBFile.read(...).get_structure(model=1)` parses the atoms of the first model.
   - `struc.filter_amino_acids(st)` keeps only amino-acid atoms, and the code then keeps a single chain (`st.chain_id == st.chain_id[0]`) so the labels line up with one sequence.
   - `struc.annotate_sse(st)` is the key step: biotite looks at the actual 3D coordinates and assigns each residue a secondary-structure letter. This is a computational reading of the experimental geometry, not a guess.
   - `struc.get_residues(st)` gives the residue names, which are converted from three-letter codes to one-letter codes via the `THREE2ONE` dictionary to build the sequence string.

3. **Labels.** biotite returns `a` (alpha helix), `b` (beta sheet), or `c` (coil). The code maps these with `m = {"a": 0, "b": 1, "c": 2}` to integer labels, defaulting anything unknown to `2` (coil). Sequences shorter than 40 residues, or where the label count does not match the residue count, are discarded (`return None`).

4. **Dataset size.** The main loop collects structures until it has 240 usable proteins (`if len(data) >= 240: break`). It then shuffles with a fixed seed (`np.random.RandomState(0)`) and splits 85% train / 15% test (`split = int(0.85 * len(data))`). So training uses roughly 240 real proteins totaling tens of thousands of labeled residues.

### Architecture (the `SSNet` class)

The network is a 1D convolutional neural network (1D-CNN). Convolutions are a natural fit here because secondary structure depends on a residue's local neighborhood, and a convolution slides a fixed-width window along the sequence.

```python
class SSNet(nn.Module):
    def __init__(self, vocab=X+1, emb=32, hid=64, k=11):
        self.emb = nn.Embedding(vocab, emb)          # 21 -> 32
        self.c1  = nn.Conv1d(emb, hid, k, padding=k//2)   # width-11 window
        self.c2  = nn.Conv1d(hid, hid, k, padding=k//2)   # width-11 window
        self.c3  = nn.Conv1d(hid, hid, 3, padding=1)      # width-3 window
        self.out = nn.Conv1d(hid, 3, 1)              # per-residue 3-way head

    def forward(self, x):
        e = self.emb(x).transpose(1, 2)   # [B,L,32] -> [B,32,L]
        h = torch.relu(self.c1(e))
        h = torch.relu(self.c2(h))
        h = torch.relu(self.c3(h))
        return self.out(h).transpose(1, 2)  # [B,3,L] -> [B,L,3]
```

Layer by layer:

- **Embedding** (`nn.Embedding(vocab=21, emb=32)`). `X = len(AA) = 20`, and `vocab = X + 1 = 21` because unknown residues are encoded as index 20 by `enc` (`AAI.get(c, X)`). Each amino acid becomes a learned 32-dimensional vector. The `transpose(1, 2)` reshapes from `[batch, length, 32]` into `[batch, 32, length]` because `Conv1d` expects channels before length.
- **Conv1d c1 and c2** each use kernel size `k = 11` with `padding = k // 2 = 5`. A width-11 window means each output position sees 5 residues on each side. `padding=5` keeps the output the same length as the input (this is "same" padding). Two stacked width-11 layers give an effective receptive field of about 21 residues, which is a realistic span for a helix or strand.
- **Conv1d c3** uses kernel size 3 (`padding=1`) to mix the features locally one more time.
- **`ReLU`** activations (`torch.relu`) add nonlinearity between conv layers.
- **Output head** `out = nn.Conv1d(hid, 3, 1)` is a width-1 convolution, which is just a per-position linear layer, producing 3 numbers (logits) per residue. The final `transpose(1, 2)` puts it back to `[batch, length, 3]`.

The output for a length-L protein is an `L x 3` array of logits. `argmax` over the last dimension gives the predicted class at each residue.

### Training

```python
model = SSNet()
opt   = torch.optim.Adam(model.parameters(), lr=1e-3)
lossf = nn.CrossEntropyLoss()
for epoch in range(40):
    for s, lab in train:
        x = torch.tensor([enc(s)], dtype=torch.long)   # batch of 1 protein
        logits = model(x)[0]                            # [L, 3]
        loss = lossf(logits, torch.tensor(lab))         # per-residue CE
        opt.zero_grad(); loss.backward(); opt.step()
```

- **Loss.** `nn.CrossEntropyLoss()` is standard multi-class classification loss. It compares the 3 logits at each residue against that residue's true class (0/1/2) and is summed/averaged over all residues. This is exactly the right loss for "pick one of three labels per position."
- **Optimizer.** Adam with learning rate `1e-3`.
- **Schedule.** 40 epochs. One protein at a time (batch size 1), with the training list reshuffled each epoch (`rs.shuffle(train)`).
- **Evaluation.** Every 8th epoch (and the last) it runs on the held-out test set and reports **Q3 accuracy**, the percentage of residues whose predicted class matches the true class:

  ```python
  pred = model(...)[0].argmax(1).numpy()
  cor += int((pred == lab).sum()); tot_res += len(lab)
  # test Q3 accuracy = 100 * cor / tot_res
  ```

  A well-trained model of this kind reports about **85% Q3 accuracy**, which is in the ballpark of classic sequence-only predictors. (Random guessing among 3 classes would be around 33%, and always guessing the most common class, coil, would do noticeably worse than 85%.)

### ONNX export

```python
torch.onnx.export(model, torch.zeros(1, 40, dtype=torch.long), "ss_model.onnx",
                  input_names=["seq"], output_names=["logits"],
                  dynamic_axes={"seq": {1: "L"}, "logits": {1: "L"}},
                  opset_version=13, dynamo=False)
```

- The model is traced with a dummy input of shape `[1, 40]` (a batch of one 40-long sequence).
- **`dynamic_axes` makes the length axis flexible.** By declaring axis 1 as a dynamic dimension named `"L"` on both `seq` and `logits`, the exported ONNX graph accepts sequences of any length, not just 40. This works cleanly here because a plain stack of `Conv1d` layers is length-agnostic (a convolution just slides further along a longer input).
- `opset_version=13`, and `dynamo=False` uses the stable tracing exporter.
- After export the script sanity-checks the ONNX file with `onnxruntime` on the CPU and asserts the ONNX output equals the PyTorch output (`ONNX == torch`), then prints an example prediction string over the alphabet `"HEC"` (Helix / shEet / Coil).

### Complexity

Let `L` be the sequence length. Each convolution touches every position once and does a constant amount of work per position (fixed kernel width and channel count), so inference is **O(L)** time and **O(L)** space. This is very cheap, which is why the app can predict structure instantly per residue.

### Where it is used in the project

The exported `ss_model.onnx` is loaded by the WinUI app and run through Windows ML (on the GPU or NPU when available). The app feeds it an aligned or single sequence and displays the per-residue helix/sheet/coil prediction, giving a structural annotation track alongside the sequence-alignment results.

---

## Model 2: Protein language model (transformer encoder)

### What it does, and what an "embedding" is

The protein language model (PLM) is a small transformer trained the way BERT is trained: masked-language modeling (MLM). It reads protein sequences and learns to fill in amino acids that have been hidden. In doing so it learns the statistical patterns of real proteins.

The payoff is the **embedding**. After training, the encoder can turn any protein into a single fixed-length vector. An embedding is that vector: a numeric summary of the whole protein. The useful property is that proteins in the same family land close together in this vector space, even when their raw sequences have diverged so far that a direct alignment would not detect the relationship. This is the "structure/function similarity that plain alignment misses" showcase of the project.

### Training data

Both the training corpus and the validation families come from UniProt via its REST streaming API (`https://rest.uniprot.org/uniprotkb/stream`).

- **Corpus.** `fetch("reviewed:true AND length:[40 TO 200]", 16000)` streams FASTA for up to **16,000** manually reviewed (Swiss-Prot) proteins of length 40 to 200. The `fetch` parser keeps a sequence only if it is 40 to 200 residues long and contains only the 20 standard amino acids (`set(cur) <= set(AA)`).
- **Labeled families for validation.** It separately downloads about 120 proteins each from five named families (`globin`, `cytochrome c`, `protein kinase`, `histone`, `ribonuclease`) using `protein_name:` queries. These labels are used only to check clustering, never for training.
- **Split.** The corpus is split 97% train / 3% validation (`split = int(0.97 * len(corpus))`).

### Tokenization and vocabulary

```python
AA = "ACDEFGHIKLMNPQRSTVWY"
TOK = {c: i for i, c in enumerate(AA)}   # A->0 ... Y->19
PAD, MASK, UNK = 20, 21, 22
VOCAB = 23
MAXLEN = 200
```

The 20 amino acids are tokens 0 to 19. Three special tokens are added: `PAD` (20) fills shorter sequences up to a common length, `MASK` (21) is the placeholder shown to the model where an amino acid was hidden, and `UNK` (22) covers any unexpected character. `encode(s)` maps a string to these integer ids.

### Architecture (the `PLM` class)

```python
class PLM(nn.Module):
    def __init__(self, d=192, nhead=6, layers=4, ff=512):
        self.emb  = nn.Embedding(VOCAB, d)          # token embedding, 23 -> 192
        self.pos  = nn.Embedding(MAXLEN + 1, d)     # positional embedding, 201 -> 192
        layer     = nn.TransformerEncoderLayer(d, nhead, ff,
                                               batch_first=True, dropout=0.1)
        self.enc  = nn.TransformerEncoder(layer, layers)   # 4 stacked layers
        self.head = nn.Linear(d, VOCAB)             # MLM head, 192 -> 23

    def hidden(self, x, pad):
        L = x.size(1)
        h = self.emb(x) + self.pos(torch.arange(L, device=x.device))[None]
        return self.enc(h, src_key_padding_mask=pad)

    def forward(self, x, pad):                       # MLM logits
        return self.head(self.hidden(x, pad))

    def embed(self, x, pad):                         # pooled embedding
        h = self.hidden(x, pad)
        w = (~pad).unsqueeze(-1).float()
        return (h * w).sum(1) / w.sum(1).clamp(min=1)
```

- **Token + positional embeddings.** Each residue token becomes a 192-dim vector (`self.emb`), and its position in the sequence becomes another 192-dim vector (`self.pos`). They are added, so the model knows both what each residue is and where it sits. `MAXLEN + 1 = 201` positional slots are allocated.
- **Transformer encoder.** `nn.TransformerEncoderLayer` with model dimension `d = 192`, `nhead = 6` attention heads, feed-forward width `ff = 512`, and `dropout = 0.1`. `batch_first=True` means tensors are shaped `[batch, length, features]`. `nn.TransformerEncoder(layer, layers=4)` stacks 4 identical layers. Self-attention lets every residue look at every other residue in the sequence, which is what a plain CNN cannot do.
- **`src_key_padding_mask=pad`.** The boolean `pad` mask (True at padded positions) tells attention to ignore padding, so padded slots do not pollute the real residues.
- **MLM head.** `nn.Linear(d, VOCAB)` turns each residue's 192-dim hidden state into 23 logits (a prediction over the vocabulary) for the training task.
- **Parameter count.** The script prints `sum(p.numel() ...) / 1e6` and reports about **1.4M parameters**, deliberately small so it trains in minutes and runs on a laptop NPU.

**What an embedding is, concretely** (the `embed` method): run the encoder, then average the hidden vectors over the real (non-pad) residues. `w = (~pad)` is 1 on real positions and 0 on padding, so `(h * w).sum(1) / w.sum(1)` is a masked mean. That single 192-dim vector is the protein's embedding. In the ONNX export it is additionally L2-normalized (see below).

### Training (masked-language modeling)

```python
model = PLM().to(D)
opt   = torch.optim.AdamW(model.parameters(), lr=3e-4)
lossf = nn.CrossEntropyLoss(ignore_index=-100)
...
real = ~pad
m = (torch.rand_like(x, dtype=torch.float) < 0.15) & real   # mask 15% of real residues
target = torch.where(m, x, torch.full_like(x, -100))         # only masked positions count
xin    = torch.where(m, torch.full_like(x, MASK), x)         # show MASK to the model
logits = model(xin, pad)
loss   = lossf(logits.reshape(-1, VOCAB), target.reshape(-1))
```

- **The task.** For each batch, 15% of the real (non-padded) residues are selected. Those positions are replaced with the `MASK` token in the input `xin`, and the model must predict the original amino acid there.
- **Loss.** `nn.CrossEntropyLoss(ignore_index=-100)`. The trick is that `target` is set to the true residue only at masked positions and to `-100` everywhere else; `ignore_index=-100` makes the loss skip all the unmasked positions, so the model is graded only on the residues it had to guess.
- **Batching.** `batchify` groups sequences into batches of 32, pads each batch to the longest sequence in it, and builds the matching `pad` mask.
- **Optimizer.** AdamW, learning rate `3e-4`.
- **Time-boxed training.** Instead of a fixed number of epochs, training runs for a wall-clock budget: `TIME_BUDGET = 720` seconds (12 minutes), looping over epochs until the budget is spent. It prefers the GPU: `device()` returns `xpu` (Intel Arc), else `cuda`, else `cpu`.
- **Validation during training.** After each epoch it reports masked-residue accuracy on the validation split (percent of masked positions predicted correctly), noting that random chance is `100/20 = 5%`. Getting well above 5% shows the model has learned real amino-acid statistics.

### Validation: do family embeddings cluster?

The real test of the embedding is downstream. Using the five labeled families, the script embeds each protein, standardizes the features, and checks clustering two ways:

1. **A PCA plot.** It computes a 2D PCA (via `np.linalg.svd`) and scatter-plots the families by color, saving `../figures/fig11_plm_embeddings.png`. Same-family proteins should visibly group together.
2. **Nearest-neighbour same-family rate.** For every protein it finds the single closest other protein by squared Euclidean distance in embedding space and checks whether that neighbour is in the same family:

   ```python
   for i in range(len(P)):
       best, bj = inf, -1
       for j in range(len(P)):
           if i == j: continue
           dsq = float(((E[i] - E[j]) ** 2).sum())
           if dsq < best: best, bj = dsq, j
       if labs[bj] == labs[i]: same += 1
   # rate = 100 * same / len(P), chance ~ 100 // len(FAMS) = 20%
   ```

   With five families, chance is about 20%. A trained model reports about **90% nearest-neighbour same-family**, meaning the embedding space really does organize proteins by family. Finally `torch.save(model.state_dict(), "plm.pt")` writes the weights.

### ONNX export as an embedder (`export_plm.py`)

`export_plm.py` reloads `plm.pt` and wraps it in an `Embedder` module whose output is exactly the vector the app wants: mean-pooled over real residues and then L2-normalized.

```python
class Embedder(nn.Module):
    def forward(self, x, mask):                 # x:[1,200] int64, mask:[1,200] float (1 = pad)
        padb = mask > 0.5
        pos  = torch.arange(FIXED, device=x.device)      # FIXED = 200
        h = self.plm.emb(x) + self.plm.pos(pos)[None]
        h = self.plm.enc(h, src_key_padding_mask=padb)
        w = (1.0 - mask).unsqueeze(-1)
        v = (h * w).sum(dim=1) / w.sum(dim=1).clamp(min=1.0)
        return v / (v.norm(dim=1, keepdim=True) + 1e-6)   # L2-normalized
```

- **L2 normalization** (`v / v.norm(...)`) scales every embedding to unit length. That makes comparing two embeddings by cosine similarity equivalent to a simple dot product, which is what the app uses to rank matches.
- **Fixed length 200 with an explicit mask.** `pad_seq` truncates or pads every sequence to exactly `FIXED = 200`, and builds a float `mask` that is `0.0` on real residues and `1.0` on padding. So the ONNX model always takes two inputs of shape `[1, 200]`: `seq` (int64) and `mask` (float32).
- **`enable_nested_tensor=False`.** The export-time `PLM` builds its `TransformerEncoder` with this flag off. Nested-tensor fast paths do dynamic, data-dependent reshaping that does not trace cleanly to a static ONNX graph, so it is disabled to keep the exported graph simple and deterministic.

```python
dx = torch.zeros(1, FIXED, dtype=torch.long)
dm = torch.zeros(1, FIXED, dtype=torch.float32)
torch.onnx.export(emb, (dx, dm), "plm_embed.onnx",
                  input_names=["seq", "mask"], output_names=["embedding"],
                  opset_version=17, dynamo=False)
```

Note there are **no `dynamic_axes` here**. The length is frozen at 200.

#### Why dynamic length fails for a traced TransformerEncoder (and why the CNN was fine)

The secondary-structure CNN could keep a dynamic length axis because a stack of convolutions is genuinely length-agnostic: the same kernel slides over however many positions exist, and no operation's shape depends on the specific value of `L`.

A traced `TransformerEncoder` is different. `torch.onnx.export` with `dynamo=False` works by tracing: it runs the model once on the dummy input and records the concrete operations that happened. Self-attention internally builds an `L x L` attention matrix and, with padding masks and the nested-tensor/fused-attention paths, performs reshapes and broadcasts whose recorded shapes are tied to the concrete `L` used during tracing. Those baked-in shapes do not generalize, so feeding a different length into a graph traced at one length leads to shape-mismatch errors or silently wrong results. The robust fix, taken here, is to pick one fixed length (200), always pad or truncate to it, and pass an explicit padding mask so the padded positions are ignored in attention and excluded from the mean pool. The output is identical to what a variable-length model would produce for the real residues, but the graph is now static and exports reliably. The script verifies this by checking `ONNX == torch embedding` with `np.allclose(..., atol=1e-4)` and prints the embedding dimension (192).

### The precomputed family database

After export, `export_plm.py` fetches about 40 proteins each from eight families (the original five plus `insulin`, `lysozyme`, `myoglobin`), embeds each one, and writes `../SequenceAlignerApp/plm_db.tsv`. Each line is `family<TAB>name<TAB>comma-separated 192-dim embedding` (formatted `%.4f`). The app loads this table so that when the user embeds a new protein it can compare against a ready-made set of known family vectors without re-running any downloads.

### Complexity

Let `L` be the sequence length (fixed at 200 for the exported model) and `d = 192` the model width. Self-attention compares every position with every other position, so a transformer layer is **O(L^2 * d)** time and **O(L^2)** space for the attention matrix, times the 4 layers. Because `L` is capped at 200 and the model is only about 1.4M parameters, a single embedding is still fast enough to run interactively on the NPU/GPU.

### Where it is used in the project

`plm_embed.onnx` runs inside the WinUI app through Windows ML (NPU or GPU). The app pads a query protein to length 200, computes its embedding, and compares it (by dot product, since embeddings are L2-normalized) against the entries in `plm_db.tsv` to find the nearest family. This is the project's demonstration that a learned embedding can recognize protein relationships that raw sequence alignment alone would miss.

---

## Both models run in the app via Windows ML

Both `ss_model.onnx` and `plm_embed.onnx` are standard ONNX files. The Python side only trains and exports them; at runtime the WinUI app hands them to Windows ML, which executes them on the best available accelerator (NPU or GPU, falling back to CPU). Training happens once on the developer machine; the app itself only performs fast inference on the exported graphs.


---

# Part 5. The Windows app


The Windows app is a WinUI 3 (Windows App SDK) desktop program in the `SequenceAlignerApp` namespace. It wraps the same alignment mathematics used in the project's Python notebook, but reimplements it in C# so it can run natively, add a rich interactive UI, draw a dot plot, show a real 3D protein structure, and run two on-device neural networks on the GPU or NPU. The relevant files are:

- `C:\Users\thomas\Desktop\Alignment\SequenceAlignerApp\Aligner.cs` (the pure algorithm engine, no UI)
- `C:\Users\thomas\Desktop\Alignment\SequenceAlignerApp\MainPage.xaml.cs` (all UI logic and event handlers)
- `C:\Users\thomas\Desktop\Alignment\SequenceAlignerApp\SsModel.cs` (secondary-structure neural net)
- `C:\Users\thomas\Desktop\Alignment\SequenceAlignerApp\EmbeddingModel.cs` (protein language model plus a small database)

---

## 1. `Aligner.cs`: the C# port of the Python engine

The file's own summary comment states its purpose directly: "C# port of the Python alignment engine: global and local alignment with linear or affine gaps, DNA or protein scoring, plus readouts and translation. Same algorithms the notebook uses." Everything in `Aligner` is `static`, so it is a stateless library of functions that the UI calls. The design goal is that for the same inputs, `Aligner.cs` produces the exact same alignments and scores as the Python code, which is why it duplicates the recurrences line for line rather than using any C# alignment library.

### 1.1 How it is structured

- **A scoring abstraction.** `public delegate double Scorer(char a, char b)` is a function type: given two symbols, return the score of pairing them. Two factories build scorers:
  - `Dna(double match, double mismatch)` returns `(a, b) => a == b ? match : -mismatch`. Note the sign: `match` is added as is, and `mismatch` is subtracted, so the slider value for mismatch is a positive penalty magnitude.
  - `FromMatrix(Dictionary<(char,char),int> m)` looks up a substitution-matrix entry. If the exact pair `(a, b)` is missing, it substitutes `'X'` (the "any residue" row/column) for whichever letter is unknown, so an out-of-alphabet character still scores instead of crashing.
- **A result record.** `public sealed record AlignResult(string Row1, string Row2, double Score, int Start1=0, int End1=0, int Start2=0, int End2=0)` holds the two aligned rows (with `'-'` for gaps), the score, and, for local alignment, the start/end coordinates in each original sequence.
- **A numerical tolerance.** `private const double Eps = 1e-9`. Because scores are `double`, traceback compares floating-point values with `Math.Abs(x - y) < Eps` instead of `==`.
- **The algorithms** (each detailed below): `Global`, `Local`, `GlobalAffine`, `SemiGlobal`, `Banded`, plus the shared linear-gap `Traceback` helper.
- **Utilities:** `ReverseComplement`, `PercentIdentity`, `GapStats`, `AlignmentScore`, `Translate`, and a private `Rev` (reverses a `StringBuilder`, used because traceback appends characters from the end of the alignment toward the front).
- **Embedded data:** the standard genetic-code codon table built by `BuildCodons`, and the `Blosum62` / `Pam250` substitution matrices parsed at load time by `ParseMatrix` from the plain-text tables `Blosum62Text` and `Pam250Text` stored at the bottom of the file.

### 1.2 Global alignment, linear gap (Needleman-Wunsch): `Global`

This is the classic Needleman-Wunsch algorithm with a constant (linear) gap penalty. `S` is an `(n+1) x (m+1)` score table. The first row and column are initialized to accumulating gap penalties, which forces the whole of both sequences to be aligned end to end.

```
S[i,0] = -i * gap          (first column: sequence 1 vs all gaps)
S[0,j] = -j * gap          (first row:    sequence 2 vs all gaps)

for i in 1..n, j in 1..m:
    diag = S[i-1,j-1] + sc(s1[i-1], s2[j-1])   # match/mismatch
    up   = S[i-1,j]   - gap                     # gap in sequence 2
    left = S[i,j-1]   - gap                     # gap in sequence 1
    S[i,j] = max(diag, up, left)
```

The final score is `S[n, m]`. The alignment is recovered by `Traceback(..., local: false)` starting from the bottom-right corner `(n, m)`.

- **Time:** O(n*m) (every cell filled once).
- **Space:** O(n*m) (the full `S` table is kept for traceback).
- **Used in:** `RunAlign` when mode is global, gaps are linear, and the sequences are not "huge."

### 1.3 Local alignment, linear gap (Smith-Waterman): `Local`

Smith-Waterman finds the single best-scoring subalignment. It differs from `Global` in two ways: every cell is floored at 0 (`Math.Max(0, ...)`), so a path can restart anywhere, and traceback begins at the highest-scoring cell rather than the corner.

```
for i in 1..n, j in 1..m:
    diag = S[i-1,j-1] + sc(s1[i-1], s2[j-1])
    up   = S[i-1,j]   - gap
    left = S[i,j-1]   - gap
    S[i,j] = max(0, diag, up, left)
    track (best, bi, bj) as the largest S seen so far
```

Traceback runs from `(bi, bj)` and stops as soon as it reaches a cell whose value is `<= Eps` (that is where the local segment began). The stop condition is the `!(local && S[i,j] <= Eps)` clause in the `while` loop. `Local` returns the best score and the local coordinates via `AlignResult(..., si, bi, sj, bj)`, where `si, sj` are the segment start indices returned by traceback and `bi, bj` are its end indices.

- **Time:** O(n*m). **Space:** O(n*m).
- **Used in:** `RunAlign` when mode index is 1 (local). If the user also selected affine gaps, the UI keeps linear here and shows the note "Affine gaps apply to global alignment only, so linear was used here."

### 1.4 The shared linear-gap traceback: `Traceback`

`Global` and `Local` share one helper. Starting from `(i, j)` it walks backward, at each step recomputing which of the three moves produced `S[i,j]` and taking the first that matches within `Eps`:

```
while (i>0 or j>0) and not (local and S[i,j] <= Eps):
    if  S[i,j] == S[i-1,j-1] + sc(s1[i-1], s2[j-1]):  diagonal (align both), i--, j--
    elif S[i,j] == S[i-1,j] - gap:                    up   (gap in seq 2),  i--
    else:                                             left (gap in seq 1),  j--
```

Because characters are appended from the end forward, the two `StringBuilder`s are reversed by `Rev` before returning. The helper also returns the final `(i, j)`, which for local alignment is the segment start.

### 1.5 Global affine gaps (Gotoh): `GlobalAffine`

An affine gap costs `open` for the first gap position and `extend` for each additional one, which models real biological indels (one long gap is cheaper than many short ones). Gotoh's algorithm uses three parallel matrices:

- `M[i,j]` best score ending with `s1[i-1]` aligned to `s2[j-1]` (a match/mismatch),
- `Ix[i,j]` best score ending with a gap in sequence 2 (a character of `s1` over a `'-'`),
- `Iy[i,j]` best score ending with a gap in sequence 1.

All cells start at negative infinity except `M[0,0] = 0`; the first column of `Ix` and first row of `Iy` are seeded with `-open - (k-1)*extend` so a leading run of gaps is charged one open plus extends. The recurrences are:

```
sub = sc(s1[i-1], s2[j-1])
M[i,j]  = max(M[i-1,j-1], Ix[i-1,j-1], Iy[i-1,j-1]) + sub
Ix[i,j] = max(M[i-1,j] - open,  Ix[i-1,j] - extend, Iy[i-1,j] - open)
Iy[i,j] = max(M[i,j-1] - open,  Iy[i,j-1] - extend, Ix[i,j-1] - open)
```

The code's own comment explains the third term in each gap recurrence: a gap can open either from a match or "right after a gap in the other sequence," which is why `Ix` can come from `Iy - open` and vice versa. The final score is `max(M[n,m], Ix[n,m], Iy[n,m])`, and traceback starts in whichever matrix (state `"M"`, `"Ix"`, or `"Iy"`) held that maximum, then follows the same max-argument logic backward, switching states as it goes.

- **Time:** O(n*m). **Space:** O(n*m) times three matrices.
- **Used in:** `RunAlign` only when mode is global and the affine gap model is selected (and the input is not huge).

### 1.6 Semi-global alignment (free end gaps): `SemiGlobal`

Semi-global (also called overlap or "glocal") alignment charges no penalty for gaps at the very start or end, which is right when one sequence is a fragment of, or overlaps, the other. It fills `S` with the same recurrence as `Global` but leaves the first row and column at 0 (no leading-gap penalty), and it takes the best score not from the corner but from anywhere along the last row or last column:

```
fill S[i,j] = max(diag, up, left)          # first row/col stay 0 -> free leading gaps
best over the last row S[n, *] and last column S[*, m]  -> (bi, bj)
```

Traceback runs from `(bi, bj)` until it reaches row 0 or column 0. The leftover leading pieces (`s1[0:si]`, `s2[0:sj]`) and trailing pieces (`s1[bi:]`, `s2[bj:]`) are then padded with `'-'` on the opposite row and stitched onto the aligned middle, so the returned rows still show the full sequences with the free end gaps made explicit.

- **Time:** O(n*m). **Space:** O(n*m).
- **Used in:** `RunAlign` when mode index is 2. Affine + semi-global falls back to linear with the same explanatory note as local.

### 1.7 Banded global alignment: `Banded`

For long, similar sequences a full O(n*m) table would use far too much memory (aligning two ~3000-character sequences means about 9 million cells per matrix). `Banded` computes only the cells within `band` of the main diagonal, so memory is O(n * band). Its comment notes it is "Exact when the best path stays in the band," which holds for similar sequences whose optimal alignment does not stray far from the diagonal.

Key details:

- **Automatic band width.** If `band < 0` (the default), it is set to `Math.Abs(n - m) + Math.Max(64, Math.Max(n, m) / 20)`: enough to cover the length difference plus a margin of at least 64, or 5% of the longer sequence, whichever is larger.
- **Rolling rows.** Only two score rows (`prev`, `cur`) of width `W = 2*band + 1` are kept and swapped each iteration. Column `c` in the band maps to sequence-2 index `j = i - band + c`. Cells with `j < 0` or `j > m` are set to negative infinity and marked direction `3` ("invalid").
- **Direction array for traceback.** A separate `byte[]` `dir` of size `(n+1) * W` records the chosen move per cell (`0` diag, `1` up, `2` left, `3` invalid). This is the one large allocation, but it is bytes, not doubles, and only `band`-wide.

Per-cell recurrence (offsets adjusted because neighbors live at shifted band columns):

```
diag = prev[c]     + sc(s1[i-1], s2[j-1])
up   = prev[c+1]   - gap        # if in range, else -inf
left = cur[c-1]    - gap        # if in range, else -inf
cur[c] = max(diag, up, left)    # record the winning direction in dir[]
```

Traceback walks from `(n, m)` reading `dir`, converting each `(ii, jj)` back to a band column `c = jj - (ii - band)`, with fallbacks if a cell falls outside the band (it still consumes a diagonal, up, or left move so the walk always terminates).

- **Time:** O(n * band). **Space:** O(n * band) for `dir`, O(band) for the score rows.
- **Used in:** `RunAlign` whenever `cells = s1.Length * s2.Length` exceeds `9_000_000`. At that size the app force-uses banded global even if the user picked local, semi-global, or affine, and shows a note saying those were skipped for size.

### 1.8 Utility functions

- **`ReverseComplement(dna)`** walks the string from the end, mapping `A<->T`, `C<->G`, `N->N`, and `U->A` (so RNA input works), leaving unknown characters unchanged. Input is upper-cased per character.
- **`PercentIdentity(r1, r2)`** counts columns where neither row is a gap (`aligned`) and how many of those match (`matches`), returning `100 * matches / aligned` (0 if nothing aligned). Gap columns are ignored, so this is identity over aligned positions.
- **`GapStats(r1, r2)`** scans for runs where either row has `'-'` and returns the number of gap runs plus a list of their lengths.
- **`AlignmentScore(r1, r2, sc, gap)`** re-scores an existing pair of rows: each gap column subtracts `gap`, each aligned column adds `sc(...)`. (This is a verification helper; the score shown in the UI comes from the DP table itself.)
- **`Translate(dna)`** upper-cases, converts `U` to `T`, then reads codons three at a time. Unknown codons become `'X'`; a stop codon (`'*'`) ends translation. The codon table is built by `BuildCodons` from the amino-acid string `"FFLLSSSS...GGGG"` in standard `TCAG x TCAG x TCAG` order (64 codons).
- **`ParseMatrix`** reads the embedded BLOSUM62/PAM250 text into a `Dictionary<(char,char),int>`, which is what `FromMatrix` looks up. Both matrices include the extended columns `B Z X *`.

---

## 2. `MainPage.xaml.cs`: the UI and every feature

`MainPage` is the single page holding all interaction. Fields worth knowing: four brushes for coloring (`GreenB` match `#2E7D46`, `RedB` mismatch `#C43D2E`, `GrayB` gap `#8A8A8A`, `SheetB` beta-sheet orange `#D98A26`); `MaxViewColumns = 1500` (the on-screen alignment is capped, though Save/Copy get the full text); a `_generation` counter used to discard stale async results; `_plainText` (the current alignment as text for copy/save); `_protein` (whether the current alignment is protein); and `_lastResult` (the most recent `AlignResult`, used by the 3D, secondary-structure, and AI features).

### 2.1 Startup and input cleaning

The constructor sets every dropdown to index 0, loads the two toy DNA samples (`SampleData.ToyDna1/2`) into the input boxes, marks `_ready = true`, and calls `RunAlign` so the app opens already showing a result. `Clean(s)` normalizes user input: it splits on newlines, skips blank lines and FASTA header lines (those starting with `>`), strips whitespace, and upper-cases every remaining character. So the boxes accept raw sequence, pasted FASTA, or multi-line text.

### 2.2 The parameter controls

`UpdateVisibility` wires the controls together: the affine sliders panel (`AffinePanel`) is only shown when the gap model is affine; the substitution-matrix dropdown is enabled only for protein; and the DNA match/mismatch sliders are enabled only for DNA. `Params_Changed` (dropdowns) and `Slider_Changed` (sliders) both re-run the alignment live; `Slider_Changed` additionally updates the numeric labels ("Match reward", "Mismatch penalty", "Gap penalty", "Gap open", "Gap extend"). So every parameter change re-aligns immediately, which is `_ready`-guarded so nothing fires during construction.

The controls map to the engine like this:

- **Type** (`TypeBox`): index 0 DNA, 1 protein. Protein uses `FromMatrix`; DNA uses `Dna(match, mismatch)`.
- **Mode** (`ModeBox`): 0 global, 1 local, 2 semi-global.
- **Matrix** (`MatrixBox`): 0 BLOSUM62, 1 PAM250 (protein only).
- **Gap model** (`GapModelBox`): 0 linear, 1 affine.
- **Sliders**: `MatchSlider`, `MismatchSlider`, `GapSlider`, `OpenSlider`, `ExtendSlider`.

### 2.3 Samples, translate, swap, reverse complement, file load

- **`Sample_Changed`** loads named pairs from `SampleData` and sets the type automatically: toy DNA, human vs gorilla HBB (protein), human vs zebrafish HBB, SARS vs SARS-CoV-2 spike DNA, the same spikes translated to protein (via `Aligner.Translate`), the two whole SARS genomes, and human vs yeast cytochrome c.
- **`Translate_Click`** cleans both boxes, replaces each with `Aligner.Translate(...)`, switches Type to protein, and re-aligns. This is the "Translate DNA to protein" action referenced by the 3D, SS, and AI features.
- **`Swap_Click`** swaps the two boxes' text and re-aligns.
- **`RevComp_Click`** replaces sequence 2 with its reverse complement (`Aligner.ReverseComplement`) and re-aligns, which is how you catch alignments that only appear on the opposite strand.
- **`LoadFile1_Click` / `LoadFile2_Click`** use `PickFastaAsync` (a `FileOpenPicker` filtered to `.fasta`, `.fa`, `.txt`, and `*`) to read a file's text into a box, then re-align.

### 2.4 Running an alignment: `RunAlign` and the auto-switch

`RunAlign` is `async void`. It increments `_generation` and captures the value as `gen`; when the background work finishes it checks `gen != _generation` and drops the result if a newer request has started. This prevents an old, slow alignment from overwriting a newer one while the user drags a slider.

It cleans both boxes, reads all parameters, builds the `Scorer`, and picks the algorithm on a background thread (`Task.Run`). The size logic:

```
cells = s1.Length * s2.Length
huge  = cells > 9_000_000     # about 3000x3000; full tables would use too much memory

if huge:                         Aligner.Banded(...)      # forced, with an explanatory note
elif affine and mode == global:  Aligner.GlobalAffine(...)
elif mode == local:              Aligner.Local(...)       # note if affine was requested
elif mode == semiglobal:         Aligner.SemiGlobal(...)  # note if affine was requested
else:                            Aligner.Global(...)
```

Status text is set to "Long sequences: running a fast banded alignment..." when huge, or "Aligning..." for a mid-size job (over 2,000,000 cells). On completion it stores `_protein`, calls `Render`, and shows any `note` (for example, that banded skipped the local/affine choice).

### 2.5 Rendering: readouts, stats, colored view, plain text

`Render` sets `_lastResult` and drives everything downstream:

- **Readouts:** `ScoreText` (the score to two decimals), `IdentityText` (`PercentIdentity`), and `GapsText` (gap count plus a bracketed list of run lengths via `Format`, which caps the display at 12 lengths and appends `,...`).
- **`UpdateStats`** writes a one-line summary: the two ungapped lengths, GC content for DNA (`Gc` counts G and C over length, omitted for protein), the number of aligned columns, and the match/mismatch counts.
- **`BuildColoredView(r1, r2)`** renders the alignment into `AlignmentView` as RichTextBlock paragraphs, 60 columns per block, in three rows: sequence 1, a match line, and sequence 2. Each column is colored green for a match, red for a mismatch, gray for a gap; the middle "match line" character is `'|'` for a match, `'.'` for a mismatch, and a space for a gap. It draws at most `MaxViewColumns` (1500) columns and, if the alignment is longer, appends a gray note telling the user to Save or Copy for the full alignment.
- **`BuildPlainText`** builds the same three-row layout as plain text (used by Copy and Save), prefixed with the local coordinates when local, and suffixed with score, percent identity, and gap stats.

### 2.6 Copy and Save

`Copy_Click` puts `_plainText` on the clipboard via a `DataPackage`. `Save_Click` opens a `FileSavePicker` (default name `alignment.txt`) and writes `_plainText`. Both use the full alignment text, not the 1500-column on-screen cap.

### 2.7 The dot plot: `DrawDotPlot`

A dot plot marks every position where a short window of sequence 1 matches a window of sequence 2; a diagonal line of dots means the two sequences are similar there. This one is drawn straight to a bitmap (comment: "so it works even for whole genomes (no skipping)"). It works on the ungapped sequences (`Row1.Replace("-","")` and the same for Row2).

- **Window size `k`** adapts to length to avoid a noisy plot: protein uses `k = 3` over 1500 residues else `2`; DNA uses `11` over 4000, `6` over 800, else `4`.
- **Method:** it indexes every length-`k` substring of sequence 2 into a `Dictionary<string, List<int>>`, then for each length-`k` window of sequence 1 looks up matching positions and paints a green pixel at `(x, y)`, where the sequence indices are scaled to the 340x340 image.
- **Output:** the pixel buffer (white BGRA background, green `#2E7D46` hits) is written into a `WriteableBitmap` and assigned to `DotImage.Source`.

### 2.8 The BLOSUM62 / PAM250 viewer

`ShowMatrix_Click` renders the selected substitution matrix (`Aligner.Pam250` or `Aligner.Blosum62`) as a monospaced, padded grid over the 20 standard amino acids and shows it in a scrollable `ContentDialog`. This lets the student see the actual scores driving a protein alignment.

---

## 3. The 3D structure viewer

The 3D viewer is a `WebView2` control (`Viewer3D`) that hosts the JavaScript library **3Dmol.js** (loaded from `https://3Dmol.org/build/3Dmol-min.js`). The structure it shows is colored by the current alignment: residues that differ between the two sequences are painted red on an otherwise light-gray model, so you literally see where the mutations are on the fold.

### 3.1 Two structure sources: `Load3D_Click`

`StructureSourceBox` chooses between predicting a structure and downloading a known one:

- **Predict (ESMFold).** This path requires a protein alignment (`_protein`); otherwise it tells the user to translate first or switch to a PDB ID. It takes sequence 2 with gaps removed and requires 10 to 400 residues (ESMFold's practical range for this app). It calls `FoldAsync(seq2)`, which POSTs the raw sequence to `https://api.esmatlas.com/foldSequence/v1/pdb/` (the ESMFold web service) with a 90-second HTTP timeout, and accepts the reply only if it contains `ATOM` records (real PDB text). The returned PDB text is stored in `_last3dPdb`, base64-encoded, and handed to 3Dmol via `viewer.addModel(atob('...'),'pdb')`.
- **From PDB ID.** Reads `PdbBox.Text` (defaulting to `6VXX`, a SARS-CoV-2 spike structure), and loads it directly with `$3Dmol.download('pdb:<id>', ...)`. Downloaded structures are not saved locally, so `_last3dPdb` is cleared in this branch.

### 3.2 Coloring by the alignment

Before loading, `Load3D_Click` walks the alignment to find the differing residues in sequence 2's own numbering: it counts only non-gap positions of `Row2` (that is the residue index `resi`), and marks a residue as changed when `Row1` has a gap there or the two letters differ. That list becomes the JavaScript array `changed`, and the viewer paints those residues red.

### 3.3 The generated page: `Build3DHtml`

`Build3DHtml` assembles a self-contained HTML string (loaded with `Viewer3D.NavigateToString`) containing the `changed` array, a 3Dmol viewer, and helper functions:

- `styleFor(rep, color)` returns the style object for the chosen representation.
- `applyStyle(rep)` styles the whole model light gray, then re-styles the `changed` residues red.
- `setSpin(on)` toggles spinning about the y axis.
- `ready()` applies the current representation, zooms to fit, and renders. It is called once the model has loaded (the `loadJs` fragment differs for predict vs download).

### 3.4 Representation, spin, and save

- **`Rep_Changed`** (dropdown `RepBox`) switches representation live by calling `window.applyStyle('cartoon'|'stick'|'sphere')` through `ExecuteScriptAsync`. `RepName()` maps index 1 to `"stick"`, 2 to `"sphere"`, and anything else to `"cartoon"`.
- **`Spin_Toggled`** (`SpinCheck`) calls `window.setSpin(true|false)`.
- **`SaveStructure_Click`** saves `_last3dPdb` to a `.pdb` file via a `FileSavePicker` (default name `structure`). It works only for a predicted structure, since downloaded structures are not held in memory; if none is loaded it says so.

All three guard on `Viewer3D.CoreWebView2 == null` and wrap the script calls in try/catch, and the whole load path catches failures with a message that the WebView2 runtime may be missing.

---

## 4. The two on-device AI features (Windows ML)

Both models are ONNX networks executed with **Windows ML** (`Windows.AI.MachineLearning`). Each tries the fast hardware first and falls back to the processor with the exact same pattern: it iterates over `(LearningModelDeviceKind.DirectXHighPerformance, "GPU / NPU")` then `(LearningModelDeviceKind.Cpu, "CPU")`, creating a `LearningModelSession` on the first device that succeeds and recording which one in a `Provider` string. `DirectXHighPerformance` targets the discrete/integrated GPU or NPU; if session creation throws, the loop simply tries CPU. If the `.onnx` file is not next to the executable (`AppContext.BaseDirectory`), the model stays unloaded and `Ready` is false, and the UI reports that gracefully. Both features are protein-only and app-only (they have no Python-notebook counterpart).

### 4.1 Secondary-structure prediction: `SsModel.cs`

`SsModel` is described in its comment as "a small neural network that predicts per-residue secondary structure (helix / sheet / coil) from a protein sequence," trained on real Protein Data Bank structures.

- **Encoding.** The amino-acid alphabet is `"ACDEFGHIKLMNPQRSTVWY"` (20 letters). `Predict(seq)` maps each residue to its index in that string; anything not found becomes index 20 (an out-of-alphabet token). The ids become an `int64` tensor of shape `[1, L]` bound to the input named `"seq"`.
- **Inference and decoding.** `_session.Evaluate(binding, "0")` produces a float output named `"logits"`, read as a flat row-major `[1 * L * 3]` vector. For each residue it takes the argmax over the 3 classes and appends `"HEC"[best]`, so the return is a string of `H` (helix), `E` (sheet), `C` (coil), one character per residue.
- **Where it is used.** `PredictSS_Click` requires a protein alignment, lazily constructs the model (`_ss ??= new SsModel()`), and runs `Predict` on sequence 2 (gaps removed) inside `Task.Run` so the UI stays responsive. It then draws sequence 2 in `SsView` with a colored secondary-structure track beneath it (60 columns per block), coloring `H` red, `E` orange (`SheetB`), and coil gray, and shows counts like "GPU / NPU: 12 helix, 4 sheet, 30 coil" where the prefix is the actual provider used.

### 4.2 Protein language model and similarity: `EmbeddingModel.cs`

`EmbeddingModel` is, per its comment, "a transformer protein language model that turns a protein into a 192-number embedding," trained from scratch on about 16,000 real proteins.

- **Fixed-length encoding.** Input is padded/truncated to `FIXED = 200`. For each position: if it is within the sequence, map the residue through the same `"ACDEFGHIKLMNPQRSTVWY"` alphabet (unknown residues become id 22) and set its mask to 0; otherwise it is padding, id 20, mask 1. Two tensors are bound, `"seq"` (`int64`, shape `[1, 200]`) and `"mask"` (`float`, shape `[1, 200]`).
- **Output.** `Evaluate` returns a float output named `"embedding"`, copied into a `float[]` (the 192-number vector).
- **Cosine similarity.** `EmbeddingModel.Cosine(a, b)` computes `dot / (|a| * |b| + 1e-6)`, a standard cosine with a tiny denominator guard.
- **The database.** `ProteinDatabase` loads `plm_db.tsv` from the executable's folder; each line is `Family <tab> Name <tab> comma-separated embedding`, parsed into `Entry(Family, Name, Emb)` records.
- **Where it is used.** `Analyze_Click` requires a protein alignment, lazily builds `_emb` and `_db`, embeds both aligned sequences (gaps removed) on a background thread, and reports: the AI cosine similarity of sequence 1 vs sequence 2 as a percentage alongside the alignment percent identity for comparison; then, if the database loaded, the top 5 most similar database proteins to sequence 2 (by cosine), and a predicted family chosen by majority vote among those top 5. The status line again names the provider actually used ("Running the transformer on the GPU / NPU..." then "Done on the ...").

---

## 5. How the pieces fit together

`Aligner.cs` is the trusted, UI-free core that mirrors the Python notebook (same five DP algorithms, same matrices, same codon table), so results are comparable across the two implementations. `MainPage.xaml.cs` is the interactive shell: it cleans input, chooses the right algorithm (including the automatic banded switch for large inputs and the generation-guarded live re-run), renders colored alignments, statistics, and a bitmap dot plot, and drives the three "extra" features. Those extras layer on top of a finished alignment: the 3D viewer visualizes sequence 2's fold with mutations highlighted (predicted by ESMFold or fetched by PDB ID through WebView2 + 3Dmol.js), and the two Windows ML models (`SsModel`, `EmbeddingModel`) run on the GPU/NPU with CPU fallback to add secondary-structure prediction and transformer-based similarity/family classification.