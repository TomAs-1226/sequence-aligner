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
