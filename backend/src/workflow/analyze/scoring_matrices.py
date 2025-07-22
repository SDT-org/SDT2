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
    "simple_3_-2": {
        "name": "Simple (3/-2)",
        "match": 3,
        "mismatch": -2,
        "description": "3:2 ratio match/mismatch"
    },
    "simple_10_-5": {
        "name": "Simple (10/-5)",
        "match": 10,
        "mismatch": -5,
        "description": "10:5 ratio match/mismatch"
    },
    "blast_5_-4": {
        "name": "BLAST (5/-4)", 
        "match": 5,
        "mismatch": -4,
        "description": "NCBI BLAST"
    },
    "megablast_2_-3": {
        "name": "MegaBLAST (2/-3)",
        "match": 2,
        "mismatch": -3,
        "description": "NCBI MegaBLAST"
    },

}

PAM_MATRICES = {
    "pam30": "PAM30",
    "pam70": "PAM70", 
    "pam120": "PAM120",
    "pam250": "PAM250"
}

BLOSUM_MATRICES = {
    "blosum45": "BLOSUM45",
    "blosum62": "BLOSUM62",
    "blosum80": "BLOSUM80",
    "blosum90": "BLOSUM90"
}

NUCLEOTIDE_MATRICES = {
    "dnafull": "DNAFULL",
    "nuc44": "NUC44"
}

ALL_MATRICES = {
    "simple": SIMPLE_MATRICES,
    "nucleotide": NUCLEOTIDE_MATRICES,
    "pam": PAM_MATRICES,
    "blosum": BLOSUM_MATRICES
}
