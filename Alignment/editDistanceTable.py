def edit_distance_table(str1: str, str2: str) -> list[list[int]]:
    """
    Returns the dynamic programming table for the edit distance between two
    strings. Entry (i, j) is the edit distance between the first i symbols of
    str1 and the first j symbols of str2.

    Parameters:
    - str1: str
    - str2: str

    Returns:
    - list[list[int]]: the edit distance table
    """
    if len(str1) == 0 or len(str2) == 0:
        raise ValueError("Error: empty string given to edit_distance_table.")

    num_rows = len(str1) + 1
    num_cols = len(str2) + 1

    scoring_matrix: list[list[int]] = []
    for i in range(num_rows):
        new_row: list[int] = []
        for j in range(num_cols):
            new_row.append(0)
        scoring_matrix.append(new_row)

    for j in range(num_cols):
        scoring_matrix[0][j] = j
    for i in range(num_rows):
        scoring_matrix[i][0] = i

    for row in range(1, num_rows):
        for col in range(1, num_cols):
            up = scoring_matrix[row - 1][col] + 1
            left = scoring_matrix[row][col - 1] + 1
            diag = scoring_matrix[row - 1][col - 1]
            if str1[row - 1] != str2[col - 1]:
                diag += 1
            scoring_matrix[row][col] = min(up, left, diag)

    return scoring_matrix


def main() -> None:
    pass


if __name__ == "__main__":
    main()
