import React from "react";
import { type ColorString, Colors } from "./colors";

export type Visualization = "histogram" | "violin" | "raincloud";
type DataSet = number[];
export type DataSets = {
  scores: DataSet;
  gc: DataSet;
  length: DataSet;
};

type VisualizationBase = {
  plotTitle: string;
  lineColor: ColorString;
  lineWidth: number;
  showGrid: boolean;
  showTickLabels: boolean;
  showAxisLabels: boolean;
};

export type DistributionState = {
  visualization: Visualization;
  dataSet: keyof DataSets;
  histogram: VisualizationBase & {
    barColor: ColorString;
    binSize: number;
    histOutlineWidth: number;
    histnorm: "probability" | "percent";
    showHistogram: boolean;
    showLine: boolean;
    histlineColor: ColorString;
  };
  raincloud: VisualizationBase & {
    bandwidth: number;
    fillColor: ColorString;
    jitter: number;
    markerColor: ColorString;
    markerSize: number;
    pointOpacity: number;
    pointPos: number;
    points: "all" | "outliers" | "suspectedoutliers" | false;
    showAxisLines: boolean;
    showPoints: boolean;
    showZeroLine: boolean;
    violinOpacity: number;
    plotOrientation: "horizontal" | "vertical";
  };
  violin: VisualizationBase & {
    bandwidth: number;
    boxOpacity: number;
    boxWidth: number;
    boxfillColor: ColorString;
    boxlineColor: ColorString;
    boxlineWidth: number;
    fillColor: ColorString;
    jitter: number;
    markerColor: ColorString;
    markerSize: number;
    pointOpacity: number;
    pointOrientation: "Violin" | "Box";
    pointPos: number;
    points: "all" | "outliers" | "suspectedoutliers" | false;
    showAxisLines: boolean;
    showBox: boolean;
    showMeanline: boolean;
    showPoints: boolean;
    showViolin: boolean;
    showZeroLine: boolean;
    violinOpacity: number;
    whiskerWidth: number;
    plotOrientation: "horizontal" | "vertical";
  };
};

const visualizationDefaults = {
  plotTitle: "Distribution of Percent Identities",
  lineColor: Colors.Tomato,
  lineWidth: 3,
  showAxisLabels: true,
  showGrid: true,
  showTickLabels: true,
};

const initialDistributionState: DistributionState = {
  visualization: "histogram",
  dataSet: "scores",
  histogram: {
    ...visualizationDefaults,
    barColor: Colors.LightBlue,
    binSize: 1,
    histOutlineWidth: 1,
    histlineColor: Colors.Tomato,
    histnorm: "probability",
    showHistogram: true,
    showLine: true,
  },
  raincloud: {
    ...visualizationDefaults,
    bandwidth: 8,
    jitter: 0.5,
    markerColor: Colors.Tomato,
    markerSize: 7,
    plotOrientation: "horizontal",
    pointOpacity: 0.5,
    pointPos: -1.5,
    points: "all",
    showAxisLines: true,
    showPoints: true,
    showZeroLine: false,
    violinOpacity: 0.5,
    fillColor: Colors.LightBlue,
  },
  violin: {
    ...visualizationDefaults,
    bandwidth: 5,
    boxOpacity: 0.5,
    boxWidth: 0.95,
    boxfillColor: Colors.LightBlue,
    boxlineColor: Colors.Tomato,
    boxlineWidth: 3,
    fillColor: Colors.LightBlue,
    jitter: 0.5,
    markerColor: Colors.Tomato,
    markerSize: 7,
    plotOrientation: "vertical",
    pointOpacity: 0.5,
    pointOrientation: "Violin",
    pointPos: 0,
    points: "all",
    showAxisLines: true,
    showBox: true,
    showMeanline: true,
    showPoints: true,
    showViolin: true,
    showZeroLine: false,
    violinOpacity: 0.5,
    whiskerWidth: 0.2,
  },
};

export const useDistributionState = () => {
  const key = "distribution-state";
  const saved = localStorage.getItem(key);
  const [state, setState] = React.useState<DistributionState>(
    saved ? JSON.parse(saved) : initialDistributionState,
  );

  React.useEffect(() => {
    localStorage.setItem(key, JSON.stringify(state));
  }, [state]);

  const updateVisualization =
    (key: DistributionState["visualization"]) =>
    (values: Partial<DistributionState["histogram"]>) =>
      setState((prev) => ({
        ...prev,
        [key]: { ...prev[key], ...values },
      }));

  return {
    state,
    setState,
    updateHistogram: updateVisualization("histogram"),
    updateRaincloud: updateVisualization("raincloud"),
    updateViolin: updateVisualization("violin"),
  };
};
