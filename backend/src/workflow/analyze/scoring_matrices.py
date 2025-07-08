SIMPLE_MATRICES = {
    "simple_1_-1": {
        "name": "Simple (1/-1)",
        "match": 1,
        "mismatch": -1,
        "description": "Basic identity matrix"
    },
    "simple_2_-1": {
        "name": "Simple (2/-1)",
        "match": 2,
        "mismatch": -1,
        "description": "2:1 ratio match/mismatch"
    },
    "blast_5_-4": {
        "name": "BLAST (5/-4)", 
        "match": 5,
        "mismatch": -4,
        "description": "NCBI BLAST default for nucleotide comparison"
    },
    "megablast_2_-3": {
        "name": "MegaBLAST (2/-3)",
        "match": 2,
        "mismatch": -3,
        "description": "NCBI MegaBLAST for highly similar sequences"
    },
    "vsearch_2_-4": {
        "name": "VSEARCH (2/-4)",
        "match": 2,
        "mismatch": -4,
        "description": "VSEARCH default for sequence clustering"
    }
}

PAM_MATRICES = {
    "pam30": "PAM30 - For closely related sequences",
    "pam70": "PAM70 - For moderately related sequences", 
    "pam120": "PAM120 - For distantly related sequences",
    "pam250": "PAM250 - For very distant sequences"
}

BLOSUM_MATRICES = {
    "blosum45": "BLOSUM45 - For distant sequences",
    "blosum62": "BLOSUM62 - General purpose (default)",
    "blosum80": "BLOSUM80 - For closely related sequences",
    "blosum90": "BLOSUM90 - For very similar sequences"
}

NUCLEOTIDE_MATRICES = {
    "dnafull": "DNAFULL - Handles ambiguous nucleotide codes",
    "nuc44": "NUC44 - Simple 4x4 nucleotide matrix"
}

ALL_MATRICES = {
    "simple": SIMPLE_MATRICES,
    "nucleotide": NUCLEOTIDE_MATRICES,
    "pam": PAM_MATRICES,
    "blosum": BLOSUM_MATRICES
}
