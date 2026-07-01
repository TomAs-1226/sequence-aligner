type Alignment = list[str]


def read_fasta_file(filename: str) -> str:
    """
    Takes a file name with a single FASTA header and reads out all elements
    as a string other than header elements.

    Parameters:
    - filename: str

    Returns:
    - str: genome string
    """
    try:
        with open(filename, 'r') as file:
            genome = ""
            for line in file:
                current_line = line.strip()
                if len(current_line) > 0 and current_line[0] != '>':
                    genome += current_line
            return genome
    except IOError:
        raise Exception("Error: something went wrong with file open (probably you gave wrong filename).")


def read_multi_fasta_file(filename: str) -> tuple[list[str], list[str]]:
    """
    Reads a FASTA file containing multiple records and returns the record
    names along with their sequences.

    Parameters:
    - filename: str

    Returns:
    - tuple[list[str], list[str]]: (names, sequences) in file order
    """
    try:
        with open(filename, 'r') as file:
            names: list[str] = []
            sequences: list[str] = []
            current_sequence = ""
            for line in file:
                current_line = line.strip()
                if len(current_line) == 0:
                    continue
                if current_line[0] == '>':
                    if len(names) > 0:
                        sequences.append(current_sequence)
                    name = current_line[1:]
                    names.append(name)
                    current_sequence = ""
                else:
                    current_sequence += current_line
            if len(names) > 0:
                sequences.append(current_sequence)
            return names, sequences
    except IOError:
        raise Exception("Error: something went wrong with file open (probably you gave wrong filename).")


def write_distance_matrix_to_file(names: list[str], matrix: list[list[int]], filename: str) -> None:
    """
    Writes a labeled distance matrix to a CSV file. The first row and first
    column hold the record names, and entry (i, j) is matrix[i][j].

    Parameters:
    - names: list[str]
    - matrix: list[list[int]]
    - filename: str
    """
    try:
        with open(filename, 'w') as file:
            header_line = "," + ",".join(names)
            file.write(header_line + "\n")
            for i in range(len(names)):
                row_values: list[str] = []
                for j in range(len(names)):
                    row_values.append(str(matrix[i][j]))
                current_line = names[i] + "," + ",".join(row_values)
                file.write(current_line + "\n")
    except IOError as e:
        raise Exception(e)


def write_alignment_to_fasta(a: Alignment, filename: str) -> None:
    """
    Takes an alignment and a file name and writes the alignment to the file as a FASTA.
    Uses "string_1" and "string_2" as the headers.

    Parameters:
    - a: Alignment
    - filename: str
    """
    try:
        with open(filename, 'w') as file:
            file.write(">string_1\n")
            file.write(a[0] + "\n")

            new_symbols = [' '] * len(a[0])
            for i in range(len(a[0])):
                if a[0][i] == '-' or a[1][i] == '-':
                    new_symbols[i] = ' '
                elif a[0][i] == a[1][i]:
                    new_symbols[i] = '|'
                else:
                    # mismatch
                    new_symbols[i] = '.'

            file.write("".join(new_symbols) + "\n")
            file.write(a[1] + "\n")
            file.write(">string_2\n")
    except IOError as e:
        raise Exception(e)


def write_local_alignment_to_fasta(a: Alignment, filename: str, start1: int, end1: int, start2: int, end2: int) -> None:
    """
    Takes an alignment and a filename. Writes the alignment to the file as a FASTA
    with indices for local alignment.

    Parameters:
    - a: Alignment
    - filename: str
    - start1: int
    - end1: int
    - start2: int
    - end2: int
    """
    try:
        with open(filename, 'w') as file:
            file.write(f">string_1: {start1} to {end1}\n")
            file.write(a[0] + "\n")

            new_symbols = [' '] * len(a[0])
            for i in range(len(a[0])):
                if a[0][i] == '-' or a[1][i] == '-':
                    new_symbols[i] = ' '
                elif a[0][i] == a[1][i]:
                    new_symbols[i] = '|'
                else:
                    # mismatch
                    new_symbols[i] = '.'

            file.write("".join(new_symbols) + "\n")
            file.write(a[1] + "\n")
            file.write(f">string_2: {start2} to {end2}\n")
    except IOError as e:
        raise Exception(e)


def print_alignment(a: Alignment) -> None:
    """
    Prints the two rows of an alignment.

    Parameters:
    - a: Alignment
    """
    print(a[0])
    print(a[1])


def print_edit_distance(name_1: str, name_2: str, distance: int) -> None:
    """
    Prints the edit distance between two named sequences.

    Parameters:
    - name_1: str
    - name_2: str
    - distance: int
    """
    print(f"Edit distance between {name_1} and {name_2}: {distance}")


def main() -> None:
    pass


if __name__ == "__main__":
    main()
