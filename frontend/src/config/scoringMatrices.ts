export interface ScoringMatrix {
  id: string;
  name: string;
  type: "nucleotide" | "protein" | "universal";
}

export const scoringMatrices: ScoringMatrix[] = [
  { id: "simple_1_-1", name: "Simple (1/-1)", type: "universal" },
  { id: "simple_2_-1", name: "Simple (2/-1)", type: "universal" },
  { id: "simple_3_-2", name: "Simple (3/-2)", type: "universal" },
  { id: "simple_10_-5", name: "Simple (10/-5)", type: "universal" },
  { id: "blast_5_-4", name: "BLAST (5/-4)", type: "nucleotide" },
  { id: "megablast_2_-3", name: "MegaBLAST (2/-3)", type: "nucleotide" },
  { id: "vsearch_2_-4", name: "VSEARCH (2/-4)", type: "nucleotide" },
  { id: "dnafull", name: "DNAFULL", type: "nucleotide" },
  { id: "nuc44", name: "NUC44", type: "nucleotide" },
  { id: "blosum45", name: "BLOSUM45", type: "protein" },
  { id: "blosum62", name: "BLOSUM62", type: "protein" },
  { id: "blosum80", name: "BLOSUM80", type: "protein" },
  { id: "blosum90", name: "BLOSUM90", type: "protein" },
  { id: "pam30", name: "PAM30", type: "protein" },
  { id: "pam70", name: "PAM70", type: "protein" },
  { id: "pam120", name: "PAM120", type: "protein" },
  { id: "pam250", name: "PAM250", type: "protein" },
];

export const getRecommendedMatrix = (isAminoAcid: boolean): string => {
  return isAminoAcid ? "blosum62" : "simple_1_-1";
};

export const getMatricesByType = (
  type: "nucleotide" | "protein" | "universal",
) => {
  return scoringMatrices.filter(
    (m) => m.type === type || m.type === "universal",
  );
};
