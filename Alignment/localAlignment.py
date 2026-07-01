from localAlignmentScore import local_score_table

type Alignment = list[str]

# Backtracking compares floating-point scores, so we allow a tiny tolerance.
EPSILON = 1e-9


def local_alignment(str1: str, str2: str, match: float, mismatch: float, gap: float) -> tuple[Alignment, int, int, int, int]:
    """
    Returns a highest-scoring local alignment of the two strings, with the start and end indices of the aligned substring in each string.

    Parameters:
    - str1: str
    - str2: str
    - match: float
    - mismatch: float
    - gap: float

    Returns:
    - tuple[Alignment, int, int, int, int]: the alignment, then start1, end1, start2, end2
      (indices are 0-based and half-open, so the aligned region of str1 is str1[start1:end1])
    """
    s = local_score_table(str1, str2, match, mismatch, gap)

    # Find the highest-scoring cell anywhere in the table. That is where the
    # best local alignment ends.
    best = 0.0
    end_r = 0
    end_c = 0
    for i in range(len(str1) + 1):
        for j in range(len(str2) + 1):
            if s[i][j] > best:
                best = s[i][j]
                end_r = i
                end_c = j

    row1 = ""  # aligned version of the str1 substring
    row2 = ""  # aligned version of the str2 substring

    r = end_r
    c = end_c

    # Walk backwards until we hit a cell whose score is 0 (the alignment starts there).
    while s[r][c] > EPSILON:
        if r > 0 and c > 0:
            if str1[r - 1] == str2[c - 1]:
                diag = s[r - 1][c - 1] + match
            else:
                diag = s[r - 1][c - 1] - mismatch
        else:
            diag = None

        up = s[r - 1][c] - gap if r > 0 else None

        if diag is not None and abs(s[r][c] - diag) < EPSILON:
            row1 = str1[r - 1] + row1
            row2 = str2[c - 1] + row2
            r -= 1
            c -= 1
        elif up is not None and abs(s[r][c] - up) < EPSILON:
            row1 = str1[r - 1] + row1
            row2 = "-" + row2
            r -= 1
        else:
            row1 = "-" + row1
            row2 = str2[c - 1] + row2
            c -= 1

    start_r = r
    start_c = c

    return [row1, row2], start_r, end_r, start_c, end_c


def main() -> None:
    pass


if __name__ == "__main__":
    main()
