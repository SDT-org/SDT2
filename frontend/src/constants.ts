export const plotFontMonospace = {
  family: '"Roboto Mono", monospace',
  weight: 400,
};

export const plotFontSansSerif = {
  family: '"Inter", sans-serif',
  weight: 400,
};

export const plotFont = plotFontMonospace;

export const saveableImageFormats = {
  svg: "SVG",
  png: "PNG",
  jpeg: "JPEG",
};

export const reorderMethods = {
  single: {
    name: "Single",
    description: "Nearest Neighbor",
  },
  complete: {
    name: "Complete",
    description: "Farthest Neighbor",
  },
  average: {
    name: "Average",
    description: "UPGMA",
  },
  weighted: {
    name: "Weighted",
    description: "WPGMA",
  },
  centroid: {
    name: "Centroid",
    description: "UPGMC",
  },
  median: {
    name: "Median",
    description: "WPGMC",
  },
  ward: {
    name: "Ward",
    description: "Minimum Variance",
  },
};

export const reorderMethodKeys = Object.keys(
  reorderMethods,
) as (keyof typeof reorderMethods)[];
