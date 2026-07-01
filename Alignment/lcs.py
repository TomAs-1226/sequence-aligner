from lcsLength import lcs_score_matrix


def longest_common_subsequence(str1: str, str2: str) -> str:
    """
    Returns a longest common subsequence of two strings.

    Parameters:
    - str1: str
    - str2: str

    Returns:
    - str: longest common subsequence
    """
    if len(str1) == 0 or len(str2) == 0:
        return ""

    scoring_matrix = lcs_score_matrix(str1, str2)

    r = len(str1)
    c = len(str2)

    lcs = ""

    while r != 0 and c != 0:
        up = scoring_matrix[r-1][c]
        left = scoring_matrix[r][c-1]
        diag = scoring_matrix[r-1][c-1]

        if str1[r-1] == str2[c-1]:
            diag += 1

        if scoring_matrix[r][c] == up:
            r -= 1
        elif scoring_matrix[r][c] == left:
            c -= 1
        elif scoring_matrix[r][c] == diag:
            if str1[r-1] == str2[c-1]:
                # match: add to LCS
                lcs = str1[r-1] + lcs
            r -= 1
            c -= 1
        else:
            raise ValueError("Error: bad scoring matrix state.")

    return lcs


def main() -> None:
    pass


if __name__ == "__main__":
    main()
