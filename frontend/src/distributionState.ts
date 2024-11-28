import React from "react";
import { ColorOption } from "./colors";

export type Visualization = "histogram" | "violin" | "raincloud";
type DataSet = number[];
export type DataSets = {
  scores: DataSet;
  gc: DataSet;
  length: DataSet;
};

type VisualizationBase = {
  plotTitle: string;
  lineColor: ColorOption;
  lineWidth: number;
  showGrid: boolean;
  showTickLabels: boolean;
  showAxisLabels: boolean;
};

export type DistributionState = {
  visualization: Visualization;
  dataSet: keyof DataSets;
  histogram: VisualizationBase & {
    barColor: ColorOption;
    binSize: number;
    histOutlineWidth: number;
    histnorm: "probability" | "percent";
    showHistogram: boolean;
    showLine: boolean;
    histlineColor: ColorOption;
  };
  raincloud: VisualizationBase & {
    bandwidth: number;
    fillColor: ColorOption;
    jitter: number;
    markerColor: ColorOption;
    markerSize: number;
    pointOpacity: number;
    pointPos: number;
    points: "all" | "outliers" | "suspectedoutliers" | "None";
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
    boxfillColor: ColorOption;
    boxlineColor: ColorOption;
    boxlineWidth: number;
    fillColor: ColorOption;
    jitter: number;
    markerColor: ColorOption;
    markerSize: number;
    pointOpacity: number;
    pointOrientation: "Violin" | "Box" | "None";
    pointPos: number;
    points: "all" | "outliers" | "suspectedoutliers" | "None";
    showAxisLines: boolean;
    showBox: boolean;
    showMeanline: "Violin" | "Box" | "None";
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
  lineColor: ColorOption.Tomato,
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
    barColor: ColorOption.Light_Blue,
    binSize: 1,
    histOutlineWidth: 1,
    histlineColor: ColorOption.Tomato,
    histnorm: "probability",
    showHistogram: true,
    showLine: true,
  },
  raincloud: {
    ...visualizationDefaults,
    bandwidth: 8,
    jitter: 0.5,
    markerColor: ColorOption.Tomato,
    markerSize: 7,
    plotOrientation: "horizontal",
    pointOpacity: 0.5,
    pointPos: -1.5,
    points: "all",
    showAxisLines: true,
    showPoints: true,
    showZeroLine: false,
    violinOpacity: 0.5,
    fillColor: ColorOption.Light_Blue,
  },
  violin: {
    ...visualizationDefaults,
    bandwidth: 5,
    boxOpacity: 0.5,
    boxWidth: 0.95,
    boxfillColor: ColorOption.Light_Blue,
    boxlineColor: ColorOption.Tomato,
    boxlineWidth: 3,
    fillColor: ColorOption.Light_Blue,
    jitter: 0.5,
    markerColor: ColorOption.Tomato,
    markerSize: 7,
    plotOrientation: "vertical",
    pointOpacity: 0.5,
    pointOrientation: "Violin",
    pointPos: 0,
    points: "all",
    showAxisLines: true,
    showBox: true,
    showMeanline: "Violin",
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
  }, [state, setState]);

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
