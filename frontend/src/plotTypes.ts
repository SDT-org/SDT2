import { colorScales } from "./colorScales";

export type Colorscale = keyof typeof colorScales;

export enum LineColor {
  White = "white",
  Black = "black",
  Red = "red",
  Blue = "blue",
  Green = "green",
}

export interface HeatmapSettings {
  colorscale: Colorscale;
  reverse: boolean;
  vmax: number;
  vmin: number;
  cellspace: number;
  annotation: boolean;
  annotation_font_size: number;
  annotation_rounding: 0 | 1 | 2;
  annotation_alpha: string;
  color: string;
  showscale:boolean,
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

export interface HistogramSettings {
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
export type HistogramData = {
  x: number[];
  y: number[];
  x2: number[];
  y2: number[];
  hover: string[];
};
