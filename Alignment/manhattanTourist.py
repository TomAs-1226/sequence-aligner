def manhattan_tourist(n: int, m: int, down: list[list[int]], right: list[list[int]]) -> int:
    """
    Returns the maximum weight of any path from the source (0, 0) to the sink
    (n, m) in a weighted grid, moving only down or right at each step.

    Parameters:
    - n: int
    - m: int
    - down: list[list[int]]
    - right: list[list[int]]

    Returns:
    - int: the weight of a longest path from (0, 0) to (n, m)
    """
    s: list[list[int]] = []
    for i in range(n + 1):
        new_row: list[int] = []
        for j in range(m + 1):
            new_row.append(0)
        s.append(new_row)

    for i in range(1, n + 1):
        s[i][0] = s[i - 1][0] + down[i - 1][0]
    for j in range(1, m + 1):
        s[0][j] = s[0][j - 1] + right[0][j - 1]

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            from_above = s[i - 1][j] + down[i - 1][j]
            from_left = s[i][j - 1] + right[i][j - 1]
            s[i][j] = max(from_above, from_left)

    return s[n][m]


def main() -> None:
    pass


if __name__ == "__main__":
    main()
