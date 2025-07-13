/**
 * LZ-ANI Parameter Legend:
 *
 * aw (Anchor Width) - Window size for approximate matching. Larger values allow more mismatches.
 * am (Anchor Mismatch) - Maximum number of mismatches allowed in the anchor window.
 * mal (Min Anchor Length) - Minimum length of exact match anchors. Shorter = more sensitive.
 * msl (Min Seed Length) - Minimum length of seed matches. Shorter = more sensitive.
 * mrd (Max Reference Distance) - Maximum distance between matches in reference sequence.
 * mqd (Max Query Distance) - Maximum distance between matches in query sequence.
 * reg (Region) - Minimum length of aligned regions to keep. Smaller = keeps shorter alignments.
 * ar (Anchor Ratio) - Minimum ratio of matched nucleotides in anchor extensions.
 *
 * General principles:
 * - Smaller mal/msl = more sensitive (finds more alignments)
 * - Larger aw/am = more tolerant of mismatches
 * - Larger mrd/mqd = allows more distant matches
 * - Smaller reg = keeps shorter alignments
 */

export interface LzaniPreset {
  id: string;
  name: string;
  description: string;
  settings: {
    aw: number;
    am: number;
    mal: number;
    msl: number;
    mrd: number;
    mqd: number;
    reg: number;
    ar: number;
  };
}

export const lzaniPresets: LzaniPreset[] = [
  {
    id: "high_sensitivity",
    name: "High Sensitivity",
    description: "SNP detection",
    settings: {
      aw: 11,
      am: 3,
      mal: 9,
      msl: 5,
      mrd: 20,
      mqd: 20,
      reg: 25,
      ar: 2,
    },
  },
  {
    id: "balanced",
    name: "Balanced",
    description: "Default",
    settings: {
      aw: 15,
      am: 7,
      mal: 11,
      msl: 7,
      mrd: 40,
      mqd: 40,
      reg: 35,
      ar: 3,
    },
  },
  {
    id: "large_genomes",
    name: "Large Genomes",
    description: "Speed optimized",
    settings: {
      aw: 20,
      am: 10,
      mal: 15,
      msl: 10,
      mrd: 60,
      mqd: 60,
      reg: 50,
      ar: 4,
    },
  },
  {
    id: "rearrangements",
    name: "Rearrangements",
    description: "Structural variants",
    settings: {
      aw: 15,
      am: 7,
      mal: 11,
      msl: 7,
      mrd: 100,
      mqd: 100,
      reg: 35,
      ar: 3,
    },
  },
  {
    id: "high_accuracy",
    name: "High Accuracy",
    description: "Conservative",
    settings: {
      aw: 13,
      am: 5,
      mal: 13,
      msl: 8,
      mrd: 30,
      mqd: 30,
      reg: 40,
      ar: 3,
    },
  },
  {
    id: "fragmented",
    name: "Fragmented",
    description: "Contigs/MAGs",
    settings: {
      aw: 15,
      am: 8,
      mal: 10,
      msl: 6,
      mrd: 50,
      mqd: 50,
      reg: 30,
      ar: 2,
    },
  },
  {
    id: "ultra_fast",
    name: "Ultra-Fast",
    description: "Quick screening",
    settings: {
      aw: 25,
      am: 12,
      mal: 20,
      msl: 12,
      mrd: 80,
      mqd: 80,
      reg: 60,
      ar: 5,
    },
  },
];

export const getPresetById = (id: string): LzaniPreset | undefined => {
  return lzaniPresets.find((preset) => preset.id === id);
};

export const getDefaultPreset = (): string => {
  return "balanced";
};
