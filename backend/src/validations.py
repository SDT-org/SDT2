def validate_fasta(filepath):
    with open(filepath) as file:
        sequence = ""
        sequence_count = 0

        for line in file:
            line = line.strip()

            if not line:
                continue
            if line.startswith(">"):
                sequence = ""
                sequence_count += 1
                continue
            sequence += line.strip()

            if len(sequence) > 50000:
                return False, "SEQUENCE_TOO_LONG"

        if sequence_count < 2:
            return False, "NOT_ENOUGH_SEQUENCES"

    return True, "VALID"
