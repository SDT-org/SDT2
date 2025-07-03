import type { PlotData } from "plotly.js";
import Plotly from "plotly.js-cartesian-dist-min";
import React from "react";
import createPlotlyComponent from "react-plotly.js/factory";
import { plotFontMonospace, plotFontSansSerif } from "../constants";
import type { DataSets, DistributionState } from "../distributionState";
import { useExportPlotly } from "../hooks/useExportPlotly";
import {
  useRelayoutHideSubtitle,
  useRelayoutUpdateTitles,
} from "../hooks/useRelayoutUpdateTitles";
import type { DistributionData, MetaData } from "../plotTypes";
import { HistogramSidebar } from "./HistogramSidebar";

const Plot = createPlotlyComponent(Plotly);

export const Histogram = ({
  data,
  dataSets,
  dataSetKey,
  metaData,
  sidebarComponent,
  settings,
  updateSettings,
}: {
  data: DistributionData | undefined;
  dataSets: DataSets;
  dataSetKey: keyof DataSets;
  metaData?: MetaData;
  sidebarComponent?: React.ReactNode;
  settings: DistributionState["histogram"];
  updateSettings: React.Dispatch<Partial<DistributionState["histogram"]>>;
}) => {
  if (!data) {
    return (
      <div className="app-main centered">
        <h1>No data available</h1>
      </div>
    );
  }

  const dataSet = dataSets[dataSetKey];

  const updateTitles = useRelayoutUpdateTitles(updateSettings);
  useRelayoutHideSubtitle(!settings.showTitles);
  const exportPlotly = useExportPlotly("distribution_histogram");

  const axisTitle = React.useMemo(() => {
    if (dataSetKey === "scores") {
      // Only use metadata labels for scores
      const labels = {
        lzani: "Average Nucleotide Identity",
        parasail: "Percent Pairwise Identity",
        default: "Percent Pairwise Identity", // fallback
      };
      const analysisMethod = metaData?.run?.analysis_method || "default";
      return labels[analysisMethod] || labels.default;
    }
    if (dataSetKey === "gc") {
      return "GC Percentage";
    }
    if (dataSetKey === "length") {
      return "Length (nt)";
    }
    return "";
  }, [dataSetKey, metaData]);

  const histogramTrace = React.useMemo(() => {
    const orientation = settings.plotOrientation;
    return {
      type: "histogram",
      ...(orientation === "vertical"
        ? { x: dataSet, xbins: { size: settings.binSize } }
        : { y: dataSet, ybins: { size: settings.binSize } }),
      histnorm: "percent",
      line: {
        width: 0,
        color: settings.histlineColor,
      },
      marker: {
        color: settings.binColor,
        line: {
          width: settings.histOutlineWidth,
          color: settings.histlineColor,
        },
      },
      name: "Histogram",
      hovertemplate:
        orientation === "vertical"
          ? "%{y:.2f}%<extra></extra>"
          : "%{x:.2f}%<extra></extra>",
    } as Partial<PlotData>;
  }, [settings, dataSet]);
  return (
    <>
      <div className="app-main">
        <Plot
          data={[histogramTrace]}
          layout={{
            ...(settings.showTitles
              ? {
                  title: {
                    text: settings.title,
                    pad: {
                      t: 100,
                      r: 0,
                      b: 0,
                      l: 0,
                    },
                  },
                }
              : {}),
            font: {
              ...(settings.titleFont === "Monospace"
                ? plotFontMonospace
                : plotFontSansSerif),
              size: 14,
            },
            uirevision: "true",
            bargap: settings.barGap,
            xaxis: {
              title: {
                text: settings.showAxisLabels
                  ? settings.plotOrientation === "vertical"
                    ? axisTitle
                    : "Distribution"
                  : "",
                font: { family: "sans-serif", size: 12 },
              },
              tickfont: { family: "sans-serif", size: 11 },
              side: "bottom",
              linecolor: settings.histlineColor,
              linewidth: settings.histOutlineWidth,
              rangemode:
                settings.plotOrientation === "vertical" ? "normal" : "tozero",
              fixedrange: true,
              zeroline: false,
              showgrid: settings.showGrid,
              showticklabels: settings.showTickLabels,
              showline: settings.showAxisLines,
              tickmode: "auto",
              autotick: true,
            },
            yaxis: {
              title: {
                text: settings.showAxisLabels
                  ? settings.plotOrientation === "horizontal"
                    ? axisTitle
                    : "Distribution"
                  : "",
                font: { family: "sans-serif", size: 12 },
              },
              tickfont: { family: "sans-serif", size: 11 },
              side: "left",
              linecolor: settings.histlineColor,
              linewidth: settings.histOutlineWidth,
              rangemode:
                settings.plotOrientation === "vertical" ? "tozero" : "normal",
              fixedrange: true,
              zeroline: false,
              showgrid: settings.showGrid,
              showticklabels: settings.showTickLabels,
              showline: settings.showAxisLines,

              tickmode: "auto",
              autotick: true,
            },
            dragmode: "pan",
            barmode: "overlay",
            margin: { l: 50, r: 50, t: 50, b: 50 },
          }}
          onRelayout={updateTitles}
          onInitialized={exportPlotly}
          config={{
            responsive: true,
            displayModeBar: false,
            scrollZoom: true,
            displaylogo: false,
            editable: false,
          }}
          style={{ width: "100%", height: "100%" }}
        />
      </div>
      <HistogramSidebar
        settings={settings}
        updateSettings={updateSettings}
        sidebarComponent={sidebarComponent}
      />
    </>
  );
};
