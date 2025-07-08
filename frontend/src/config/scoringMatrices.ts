export interface ScoringMatrix {
  id: string;
  name: string;
  type: "nucleotide" | "protein" | "universal";
  description?: string;
}

export const scoringMatrices: ScoringMatrix[] = [
  // Simple matrices
  {
    id: "simple_1_-1",
    name: "Simple (1/-1)",
    type: "universal",
    description: "Basic identity matrix",
  },
  {
    id: "simple_2_-1",
    name: "Simple (2/-1)",
    type: "universal",
    description: "2:1 ratio match/mismatch",
  },
  {
    id: "blast_5_-4",
    name: "BLAST (5/-4)",
    type: "nucleotide",
    description: "NCBI BLAST default for nucleotide comparison",
  },
  {
    id: "megablast_2_-3",
    name: "MegaBLAST (2/-3)",
    type: "nucleotide",
    description: "NCBI MegaBLAST for highly similar sequences",
  },
  {
    id: "vsearch_2_-4",
    name: "VSEARCH (2/-4)",
    type: "nucleotide",
    description: "VSEARCH default for sequence clustering",
  },
  // Nucleotide matrices
  {
    id: "dnafull",
    name: "DNAFULL",
    type: "nucleotide",
    description: "Handles ambiguous nucleotide codes",
  },
  {
    id: "nuc44",
    name: "NUC44",
    type: "nucleotide",
    description: "Simple 4x4 nucleotide matrix",
  },
  // BLOSUM matrices
  {
    id: "blosum45",
    name: "BLOSUM45",
    type: "protein",
    description: "For distant sequences",
  },
  {
    id: "blosum62",
    name: "BLOSUM62",
    type: "protein",
    description: "General purpose (default)",
  },
  {
    id: "blosum80",
    name: "BLOSUM80",
    type: "protein",
    description: "For closely related sequences",
  },
  {
    id: "blosum90",
    name: "BLOSUM90",
    type: "protein",
    description: "For very similar sequences",
  },
  // PAM matrices
  {
    id: "pam30",
    name: "PAM30",
    type: "protein",
    description: "For closely related sequences",
  },
  {
    id: "pam70",
    name: "PAM70",
    type: "protein",
    description: "For moderately related sequences",
  },
  {
    id: "pam120",
    name: "PAM120",
    type: "protein",
    description: "For distantly related sequences",
  },
  {
    id: "pam250",
    name: "PAM250",
    type: "protein",
    description: "For very distant sequences",
  },
];

export const getRecommendedMatrix = (isAminoAcid: boolean): string => {
  return isAminoAcid ? "blosum62" : "blast_5_-4";
};

export const getMatricesByType = (
  type: "nucleotide" | "protein" | "universal",
) => {
  return scoringMatrices.filter(
    (m) => m.type === type || m.type === "universal",
  );
};
