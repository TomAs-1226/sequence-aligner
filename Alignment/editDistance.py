from editDistanceTable import edit_distance_table


def edit_distance(str1: str, str2: str) -> int:
    """
    Returns the edit distance between the two strings: the minimum number of
    substitutions, insertions, and deletions to turn one into the other.

    Parameters:
    - str1: str
    - str2: str

    Returns:
    - int: edit distance
    """
    if len(str1) == 0 or len(str2) == 0:
        raise ValueError("Error: empty string given to edit_distance.")

    scoring_matrix = edit_distance_table(str1, str2)
    return scoring_matrix[len(str1)][len(str2)]


def main() -> None:
    pass


if __name__ == "__main__":
    main()
