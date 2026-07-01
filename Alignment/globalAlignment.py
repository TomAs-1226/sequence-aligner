from globalAlignmentScore import global_score_table

type Alignment = list[str]

# Backtracking compares floating-point scores, so we allow a tiny tolerance.
EPSILON = 1e-9


def global_alignment(str1: str, str2: str, match: float, mismatch: float, gap: float) -> Alignment:
    """
    Returns a highest-scoring global alignment of the two strings as two equal-length rows (using '-' for indels).

    Parameters:
    - str1: str
    - str2: str
    - match: float
    - mismatch: float
    - gap: float

    Returns:
    - Alignment: the two rows of a maximum-scoring global alignment
    """
    s = global_score_table(str1, str2, match, mismatch, gap)

    row1 = ""  # aligned version of str1
    row2 = ""  # aligned version of str2

    r = len(str1)
    c = len(str2)

    # Walk backwards from the bottom-right corner to the top-left corner.
    while r > 0 or c > 0:
        if r > 0 and c > 0:
            if str1[r - 1] == str2[c - 1]:
                diag = s[r - 1][c - 1] + match
            else:
                diag = s[r - 1][c - 1] - mismatch
        else:
            diag = None

        up = s[r - 1][c] - gap if r > 0 else None
        # left = s[r][c - 1] - gap if c > 0 else None  (handled implicitly below)

        if diag is not None and abs(s[r][c] - diag) < EPSILON:
            # diagonal move: line up the two characters
            row1 = str1[r - 1] + row1
            row2 = str2[c - 1] + row2
            r -= 1
            c -= 1
        elif up is not None and abs(s[r][c] - up) < EPSILON:
            # up move: a character of str1 against a gap
            row1 = str1[r - 1] + row1
            row2 = "-" + row2
            r -= 1
        else:
            # left move: a gap against a character of str2
            row1 = "-" + row1
            row2 = str2[c - 1] + row2
            c -= 1

    return [row1, row2]


def main() -> None:
    pass


if __name__ == "__main__":
    main()
