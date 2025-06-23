export default {
  INVALID_FILE_TYPE:
    "This file is not compatible with SDT2. Please select a file with an extension of .fasta or .csv.",
  SEQUENCE_TOO_LONG:
    "This file contains one or more sequences that are too long (> 50000 characters.) Attempting to process this file could cause system instability.",
  NOT_ENOUGH_SEQUENCES: "File must contain 2 or more sequences.",
  DUPLICATE_SEQUENCE_NAME:
    "This file contains duplicate sequence names. Sequence names must be unique.",
  ZERO_LENGTH_SEQUENCE:
    "This file contains a zero length sequence. Sequences cannot be empty.",
  AMINO_ACID_NOT_SUPPORTED_LZANI:
    "LZ-ANI only supports nucleotide ANI calculations at this time.",
  PARASAIL_PERFORMANCE_WARNING:
    "Using Parasail for high sequence count datasets (>500 sequences) or datasets with long sequences (>30,000 nt) may cause performance issues. The LZ-ANI aligner is better suited for these cases.", // TODO: remove harded-coded values and show in performance warning?
  INVALID_FASTA_FORMAT: "The file is not a valid FASTA format.",
};
