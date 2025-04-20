def validate_fasta(filepath, filetype=None):
    if filetype and filetype != "text/fasta":
        return False, "INVALID_FILE_TYPE"
    with open(filepath) as file:
        sequence = ""
        sequence_name = ""
        sequence_count = 0
        sequence_names = []

        for line in file:
            line = line.strip()

            if not line:
                continue
            if line.startswith(">"):
                if line in sequence_names:
                    return False, "DUPLICATE_SEQUENCE_NAME"
                if (
                    sequence_name
                    and line != sequence_name
                    and not any(char.isalpha() for char in sequence)
                ):
                    return False, "ZERO_LENGTH_SEQUENCE"
                sequence = ""
                sequence_name = line
                sequence_names.append(sequence_name)
                sequence_count += 1
                continue

            sequence += line.strip()

        if sequence_count < 2:
            return False, "NOT_ENOUGH_SEQUENCES"

    return True, "VALID"
