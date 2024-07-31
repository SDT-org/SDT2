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
};
