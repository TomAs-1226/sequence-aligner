def local_score_table(str1: str, str2: str, match: float, mismatch: float, gap: float) -> list[list[float]]:
    """
    Returns the local alignment dynamic programming score table for the two strings (entries are floored at 0).

    Parameters:
    - str1: str
    - str2: str
    - match: float
    - mismatch: float
    - gap: float

    Returns:
    - list[list[float]]: the local alignment score table
    """
    num_rows = len(str1) + 1
    num_cols = len(str2) + 1

    # initialization: zero-filled table. For local alignment the first row and
    # first column stay 0, which is what lets an alignment start anywhere.
    s: list[list[float]] = []
    for i in range(num_rows):
        new_row: list[float] = []
        for j in range(num_cols):
            new_row.append(0.0)
        s.append(new_row)

    # fill the rest of the table, never letting a score drop below 0
    for row in range(1, num_rows):
        for col in range(1, num_cols):
            up = s[row - 1][col] - gap
            left = s[row][col - 1] - gap
            if str1[row - 1] == str2[col - 1]:
                diag = s[row - 1][col - 1] + match
            else:
                diag = s[row - 1][col - 1] - mismatch
            s[row][col] = max(0.0, up, left, diag)

    return s


def main() -> None:
    pass


if __name__ == "__main__":
    main()
