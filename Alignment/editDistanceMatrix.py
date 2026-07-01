from editDistance import edit_distance


def edit_distance_matrix(patterns: list[str]) -> list[list[int]]:
    """
    Returns a matrix whose (i,j)th entry is the edit distance between
    the i-th and j-th strings in patterns.

    Parameters:
    - patterns: list[str]

    Returns:
    - list[list[int]]: pairwise edit distance matrix
    """
    num_rows = len(patterns)
    num_cols = len(patterns)
    mtx = initialize_integer_matrix(num_rows, num_cols)

    for i in range(num_rows):
        for j in range(i + 1, num_cols):
            d = edit_distance(patterns[i], patterns[j])
            mtx[i][j] = d
            mtx[j][i] = d

    return mtx


def initialize_integer_matrix(num_rows: int, num_cols: int) -> list[list[int]]:
    """
    Returns a zero-filled integer matrix with the given dimensions.

    Parameters:
    - num_rows: int
    - num_cols: int

    Returns:
    - list[list[int]]: zero-filled matrix
    """
    mtx: list[list[int]] = []
    for i in range(num_rows):
        row: list[int] = []
        for j in range(num_cols):
            row.append(0)
        mtx.append(row)
    return mtx


def main() -> None:
    pass


if __name__ == "__main__":
    main()
