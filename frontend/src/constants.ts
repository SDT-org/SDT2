export const plotFontMonospace = {
  family: '"Roboto Mono", monospace',
  weight: 400,
};

export const plotFontSansSerif = {
  family: '"Inter", sans-serif',
  weight: 400,
};

export const plotFont = plotFontMonospace;

export const reorderMethods = {
  single: {
    name: "Nearest Neighbor",
    description: "Single Linkage method",
  },
  complete: {
    name: "Farthest Neighbor",
    description: "Complete Linkage method",
  },
  average: {
    name: "UPGMA",
    description: "Average Linkage method",
  },
  weighted: {
    name: "WPGMA",
    description: "Weighted Linkage method",
  },
  centroid: {
    name: "UPGMC",
    description: "Centroid Linkage method",
  },
  median: {
    name: "WPGMC",
    description: "Median Linkage method",
  },
  ward: {
    name: "Minimum Variance",
    description: "Ward's Linkage method",
  },
};

export const reorderMethodKeys = Object.keys(
  reorderMethods,
) as (keyof typeof reorderMethods)[];
