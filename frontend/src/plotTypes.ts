import { colorScales } from "./colorScales";

export type Colorscale = keyof typeof colorScales;

export enum LineColor {
  White = "white",
  Black = "black",
  Red = "red",
  Blue = "blue",
  Green = "green",
}
export enum DistributionPlots {
  Bar = "Bar Plot",
  Scatter = "Scatter PLot",
  Box = "Box Plot",
  Violin = "Violin Plot",
  Contour = "Contour Plot",
}
export interface HeatmapSettings {
  colorscale: Colorscale | "Discrete";
  reverse: boolean;
  vmax: number;
  vmin: number;
  cellspace: number;
  annotation: boolean;
  annotation_font_size: number;
  annotation_rounding: 0 | 1 | 2;
  annotation_alpha: string;
  color: string;
  showscale: boolean;
  cbar_shrink: number;
  cbar_pad: string;
  cbar_aspect: string;
  axis_labels: boolean;
  axlabel_xrotation: number;
  axlabel_xfontsize: number;
  axlabel_yrotation: number;
  axlabel_yfontsize: number;
  cutoff_1: number;
  cutoff_2: number;
}
export type HeatmapData = string[][];

export type GetDataResponse = {
  data: string[][];
  metadata: {
    minVal: number;
    maxVal: number;
  };
  gc_stats: number[];
  length_stats: number[];
  identity_scores: [string, string, number][];
};

export type DistributionData = Omit<
  GetDataResponse,
  "data" | "identity_scores"
> & {
  raw_mat: number[];
  identity_combos: [string, string][];
};
