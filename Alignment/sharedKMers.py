def count_shared_kmers(str1: str, str2: str, k: int) -> int:
    """
    Returns the number of k-mers shared by the two strings, counting multiplicity (the minimum count of each shared k-mer).

    Parameters:
    - str1: str
    - str2: str
    - k: int

    Returns:
    - int: number of shared k-mers
    """
    if k <= 0 or k > len(str1) or k > len(str2):
        return 0

    counts1: dict[str, int] = {}
    for i in range(len(str1) - k + 1):
        kmer = str1[i:i + k]
        counts1[kmer] = counts1.get(kmer, 0) + 1

    counts2: dict[str, int] = {}
    for i in range(len(str2) - k + 1):
        kmer = str2[i:i + k]
        counts2[kmer] = counts2.get(kmer, 0) + 1

    shared = 0
    for kmer in counts1:
        if kmer in counts2:
            shared += min(counts1[kmer], counts2[kmer])

    return shared


def main() -> None:
    pass


if __name__ == "__main__":
    main()
