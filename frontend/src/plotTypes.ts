import { colorScales } from "./colorScales";

export type Colorscale = keyof typeof colorScales;

interface DataToPlot {
  x: number[];
  y: number[];
  gcContent:number[];
  length:number[];
}
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

export interface DistributionSettings {
  showHistogram: boolean;
  showLinePlot: boolean;
  lineColor: string;
  lineWidth: number;
  lineShape: string;
  barColor: string;
  barlineColor: string;
  barOutlineWidth: number;
  markerSymbol: string;
  markerSize: number;
  markerColor: string;
  plotTitle: string;
  showGrid: boolean;
  showLine: boolean;
  showZeroLine: boolean;
  showScatterPlot: boolean;
  showTickLabels: boolean;
  showAxisLabels: boolean;
}
export type DistributionData = {
  x: number[];
  y: number[];
  x2: number[];
  y2: number[];
  raw_mat:number[];
  round_mat:number[];
  gc_stats: number[];
  len_stats: number[];
  hover: string[];
};
