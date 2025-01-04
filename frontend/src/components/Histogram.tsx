import type { PlotData } from "plotly.js";
import Plotly from "plotly.js-cartesian-dist-min";
import React from "react";
import createPlotlyComponent from "react-plotly.js/factory";
import { plotFontMonospace, plotFontSansSerif } from "../constants";
import type { DataSets, DistributionState } from "../distributionState";
import {
  useRelayoutHideSubtitle,
  useRelayoutUpdateTitles,
} from "../hooks/useRelayoutUpdateTitles";
import type { DistributionData } from "../plotTypes";
import { HistogramSidebar } from "./HistogramSidebar";

const Plot = createPlotlyComponent(Plotly);

export const Histogram = ({
  data,
  dataSets,
  dataSetKey,
  sidebarComponent,
  settings,
  updateSettings,
}: {
  data: DistributionData | undefined;
  dataSets: DataSets;
  dataSetKey: keyof DataSets;
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

  const histogramTrace = React.useMemo(() => {
    const orientation = settings.plotOrientation;
    return {
      type: "histogram",
      ...(orientation === "vertical"
        ? { x: dataSet, xbins: { size: settings.binSize } }
        : { y: dataSet, ybins: { size: settings.binSize } }),
      histnorm: "percent",
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
                    text:
                      settings.title +
                      (settings.subtitle
                        ? `<br><span style="font-size:0.8em;">${settings.subtitle}</span>`
                        : ""),
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
              // @ts-ignore
              weight: "bold",
            },
            uirevision: "true",
            xaxis: {
              ...(settings.showTitles
                ? {
                    title: {
                      text: settings.xtitle,
                      font: {
                        ...(settings.titleFont === "Monospace"
                          ? plotFontMonospace
                          : plotFontSansSerif),
                        //@ts-ignore
                        weight: "bold",
                      },
                    },
                    scaleratio: 1,
                  }
                : {}),
              side: "bottom",
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
              ...(settings.showTitles
                ? {
                    title: {
                      text: settings.ytitle,
                      font: {
                        ...(settings.titleFont === "Monospace"
                          ? plotFontMonospace
                          : plotFontSansSerif),
                        //@ts-ignore
                        weight: "bold",
                      },
                      pad: {
                        r: 15,
                      },
                    },
                  }
                : {}),
              side: "left",
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
