def lcs_length(str1: str, str2: str) -> int:
    """
    Returns the length of a longest common subsequence of two strings.

    Parameters:
    - str1: str
    - str2: str

    Returns:
    - int: LCS length
    """
    if len(str1) == 0 or len(str2) == 0:
        raise ValueError("Error: empty string given to lcs_length.")

    scoring_matrix = lcs_score_matrix(str1, str2)
    return scoring_matrix[len(str1)][len(str2)]


def lcs_score_matrix(str1: str, str2: str) -> list[list[int]]:
    """
    Returns the dynamic programming scoring matrix for LCS.

    Parameters:
    - str1: str
    - str2: str

    Returns:
    - list[list[int]]: scoring matrix
    """
    if len(str1) == 0 or len(str2) == 0:
        raise ValueError("Error: empty string given to lcs_score_matrix.")

    num_rows = len(str1) + 1
    num_cols = len(str2) + 1

    # initialization
    scoring_matrix: list[list[int]] = []
    for i in range(num_rows):
        new_row: list[int] = []
        for j in range(num_cols):
            new_row.append(0)
        scoring_matrix.append(new_row)

    # set all the values of the matrix
    for row in range(1, num_rows):
        for col in range(1, num_cols):
            up = scoring_matrix[row - 1][col]
            left = scoring_matrix[row][col - 1]
            diag = scoring_matrix[row - 1][col - 1]
            if str1[row - 1] == str2[col - 1]:
                diag += 1
            scoring_matrix[row][col] = max(up, left, diag)

    return scoring_matrix


def main() -> None:
    pass


if __name__ == "__main__":
    main()
