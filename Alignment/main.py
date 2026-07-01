from io_python import read_fasta_file, read_multi_fasta_file, write_alignment_to_fasta, write_local_alignment_to_fasta, write_distance_matrix_to_file, print_edit_distance
from localAlignment import local_alignment
from globalAlignment import global_alignment
from editDistance import edit_distance
from lcsLength import lcs_length
from editDistanceMatrix import edit_distance_matrix


def main() -> None:
    print("Sequence alignment!")
    sars_lcs_and_edit_distance()
    hemoglobin_edit_distance_matrix()


def sars_lcs_and_edit_distance() -> None:
    """Computes the LCS length and the edit distance between the SARS-CoV (2003) and SARS-CoV-2 (2020) genomes."""
    print("Reading in genomes.")
    sars1 = read_fasta_file("Data/Coronaviruses/SARS-CoV_genome.fasta")
    sars2 = read_fasta_file("Data/Coronaviruses/SARS-CoV-2_genome.fasta")

    print(f"SARS-CoV genome length:   {len(sars1)}")
    print(f"SARS-CoV-2 genome length: {len(sars2)}")

    print("Computing LCS length (this is a big computation on whole genomes).")
    shared = lcs_length(sars1, sars2)
    print(f"LCS length of the two genomes: {shared}")

    print("Computing edit distance (this is a big computation on whole genomes).")
    distance = edit_distance(sars1, sars2)
    print(f"Edit distance between the two genomes: {distance}")


def hemoglobin_edit_distance_matrix() -> None:
    """Builds the pairwise edit distance matrix for a multi-species hemoglobin protein dataset and writes it to a file."""
    names, sequences = read_multi_fasta_file("Data/Hemoglobin/hemoglobin_protein.fasta")
    print(f"Read {len(names)} hemoglobin sequences.")

    print("Building the pairwise edit distance matrix.")
    matrix = edit_distance_matrix(sequences)

    write_distance_matrix_to_file(names, matrix, "Output/hemoglobin_edit_distance_matrix.csv")
    print("Edit distance matrix written to Output/hemoglobin_edit_distance_matrix.csv.")


def hemoglobin_edit_distance() -> None:
    """Prints the edit distance between human hemoglobin and each of gorilla, cow, and zebrafish."""
    human = read_fasta_file("Data/Hemoglobin/Homo_sapiens_hemoglobin.fasta")
    gorilla = read_fasta_file("Data/Hemoglobin/Gorilla_gorilla_hemoglobin.fasta")
    cow = read_fasta_file("Data/Hemoglobin/Bos_taurus_hemoglobin.fasta")
    zebrafish = read_fasta_file("Data/Hemoglobin/Danio_rerio_hemoglobin.fasta")

    print_edit_distance("human", "gorilla", edit_distance(human, gorilla))
    print_edit_distance("human", "cow", edit_distance(human, cow))
    print_edit_distance("human", "zebrafish", edit_distance(human, zebrafish))


def sars_local_spike_alignment() -> None:
    """Computes and writes a local alignment of the SARS-CoV spike protein against the SARS-CoV-2 genome."""
    sars_spike = read_fasta_file("Data/Coronaviruses/SARS-CoV_genome_spike_protein.fasta")
    sars2 = read_fasta_file("Data/Coronaviruses/SARS-CoV-2_genome.fasta")

    match = 1.0
    mismatch = 1.0
    gap = 1.0

    spike_alignment, start1, end1, start2, end2 = local_alignment(sars_spike, sars2, match, mismatch, gap)
    write_local_alignment_to_fasta(spike_alignment, "Output/sars_genome_spike_alignment_GOOD.txt", start1, end1, start2, end2)


def sars_spike_alignment() -> None:
    """Computes and writes a global alignment of the SARS-CoV spike protein against the SARS-CoV-2 genome."""
    sars_spike = read_fasta_file("Data/Coronaviruses/SARS-CoV_genome_spike_protein.fasta")
    sars2 = read_fasta_file("Data/Coronaviruses/SARS-CoV-2_genome.fasta")

    match = 1.0
    mismatch = 1.0
    gap = 2.0

    spike_alignment = global_alignment(sars_spike, sars2, match, mismatch, gap)
    write_alignment_to_fasta(spike_alignment, "Output/sars_genome_spike_alignment.txt")


def sars_alignment() -> None:
    """Computes and writes a global alignment of the two SARS genomes."""
    print("Reading in genomes.")
    sars1 = read_fasta_file("Data/Coronaviruses/SARS-CoV_genome.fasta")
    sars2 = read_fasta_file("Data/Coronaviruses/SARS-CoV-2_genome.fasta")

    print("Genomes read. Running global alignment.")
    match = 1.0
    mismatch = 1.0
    gap = 2.0

    sars_aln = global_alignment(sars1, sars2, match, mismatch, gap)
    print("Global alignment run. Printing alignment to file.")
    write_alignment_to_fasta(sars_aln, "Output/sars_genome_alignment.txt")


def hemoglobin() -> None:
    """Computes and writes global alignments of human hemoglobin against zebrafish, cow, and gorilla."""
    zebrafish = read_fasta_file("Data/Hemoglobin/Danio_rerio_hemoglobin.fasta")
    human = read_fasta_file("Data/Hemoglobin/Homo_sapiens_hemoglobin.fasta")
    cow = read_fasta_file("Data/Hemoglobin/Bos_taurus_hemoglobin.fasta")
    gorilla = read_fasta_file("Data/Hemoglobin/Gorilla_gorilla_hemoglobin.fasta")

    match = 1.0
    mismatch = 1.0
    gap = 3.0

    alignment1 = global_alignment(zebrafish, human, match, mismatch, gap)
    write_alignment_to_fasta(alignment1, "Output/zebrafish_human_hemoglobin.txt")

    alignment2 = global_alignment(cow, human, match, mismatch, gap)
    write_alignment_to_fasta(alignment2, "Output/cow_human_hemoglobin.txt")

    alignment3 = global_alignment(gorilla, human, match, mismatch, gap)
    write_alignment_to_fasta(alignment3, "Output/gorilla_human_hemoglobin.txt")


if __name__ == "__main__":
    main()
