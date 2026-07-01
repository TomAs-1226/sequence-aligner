def global_score_table(str1: str, str2: str, match: float, mismatch: float, gap: float) -> list[list[float]]:
    """
    Returns the global alignment dynamic programming score table for the two strings under the given match reward, mismatch penalty, and gap penalty.

    Parameters:
    - str1: str
    - str2: str
    - match: float
    - mismatch: float
    - gap: float

    Returns:
    - list[list[float]]: the global alignment score table
    """
    num_rows = len(str1) + 1
    num_cols = len(str2) + 1

    # initialization: zero-filled table
    s: list[list[float]] = []
    for i in range(num_rows):
        new_row: list[float] = []
        for j in range(num_cols):
            new_row.append(0.0)
        s.append(new_row)

    # first column and first row: every step is a gap, so the cost piles up
    for i in range(1, num_rows):
        s[i][0] = s[i - 1][0] - gap
    for j in range(1, num_cols):
        s[0][j] = s[0][j - 1] - gap

    # fill the rest of the table
    for row in range(1, num_rows):
        for col in range(1, num_cols):
            up = s[row - 1][col] - gap
            left = s[row][col - 1] - gap
            if str1[row - 1] == str2[col - 1]:
                diag = s[row - 1][col - 1] + match
            else:
                diag = s[row - 1][col - 1] - mismatch
            s[row][col] = max(up, left, diag)

    return s


def main() -> None:
    pass


if __name__ == "__main__":
    main()
