def validate_fasta(filepath):
    with open(filepath) as file:
        sequence = ""
        for line in file:
            line = line.strip()

            if not line:
                continue
            if line.startswith(">"):
                sequence = ""
                continue

            sequence += line.strip()
            if len(sequence) > 50000:
                return False, "SEQUENCE_TOO_LONG"
    return True, "VALID"
