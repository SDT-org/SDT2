import type { PlotData } from "plotly.js";
import Plotly from "plotly.js-cartesian-dist-min";
import React from "react";
import createPlotlyComponent from "react-plotly.js/factory";

import { plotFontMonospace, plotFontSansSerif } from "../constants";
import type { DataSets, DistributionState } from "../distributionState";
import { arrayMinMax } from "../helpers";
import {
  useRelayoutHideSubtitle,
  useRelayoutUpdateTitles,
} from "../hooks/useRelayoutUpdateTitles";
import type { DistributionData, MetaData } from "../plotTypes";
import { RaincloudSidebar } from "./RaincloudSiderbar";
const Plot = createPlotlyComponent(Plotly);

export const Raincloud = ({
  data,
  dataSets,
  dataSetKey,
  settings,
  updateSettings,
  sidebarComponent,
}: {
  data: DistributionData | undefined;
  dataSets: DataSets;
  dataSetKey: keyof DataSets;
  metaData: MetaData;
  footer?: React.ReactNode;
  sidebarComponent?: React.ReactNode;
  settings: DistributionState["raincloud"];
  updateSettings: React.Dispatch<Partial<DistributionState["raincloud"]>>;
}) => {
  if (!data) {
    return (
      <div className="app-main centered">
        <h1>No data available</h1>
      </div>
    );
  }
  const dataSet = dataSets[dataSetKey];
  const [minDataValue, maxDataValue] = arrayMinMax(dataSet);
  const updateTitles = useRelayoutUpdateTitles(updateSettings);
  useRelayoutHideSubtitle(!settings.showTitles);
  const rainCloudTrace = React.useMemo(
    () =>
      ({
        type: "violin",
        name: "",
        x: dataSet,
        side: "positive",
        points: settings.showPoints ? settings.points : false,
        showPoints: true,
        line: {
          color: settings.lineColor,
          width: settings.lineWidth,
        },
        fillcolor: settings.fillColor,
        pointpos: settings.pointPos,
        jitter: settings.jitter,
        bandwidth: settings.bandwidth,
        marker: {
          visible: true,
          color: settings.markerColor,
          size: settings.markerSize,
        },
        meanline: {
          visible: settings.showMeanline,
        },
        hoveron: "points",
        hovertemplate: "%{text}<br> <br>Percent Identity: %{x}<extra></extra>",
        text: data.identity_combos.map(
          (ids) => `Seq 1: ${ids[0]}<br>Seq 2: ${ids[1]}`,
        ),
      }) as Partial<PlotData>,
    [data.identity_combos, dataSet, settings],
  );
  return (
    <>
      <div className="app-main">
        <Plot
          data={[rainCloudTrace]}
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
            uirevision: "true",
            font: {
              ...(settings.titleFont === "Monospace"
                ? plotFontMonospace
                : plotFontSansSerif),
              //@ts-ignore
              weight: "bold",
            },
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
              side: "left",
              rangemode: "tozero",
              fixedrange: true,
              zeroline: false,
              dtick: 5,
              showgrid: settings.showGrid,
              showline: settings.showAxisLines,
              showticklabels: settings.showTickLabels,
              range: [minDataValue - 20, maxDataValue + 20],
              tickmode: "array",
              tickfont: {
                ...plotFontMonospace,
              },
              ticktext: Array.from(
                { length: Math.ceil((maxDataValue + 20) / 5) + 1 },
                (_, i) => {
                  const value = i * 5;
                  return value <= 100 ? value.toString() : "";
                },
              ),
              tickvals: Array.from(
                { length: Math.ceil((maxDataValue + 20) / 5) + 1 },
                (_, i) => i * 5,
              ),
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
              fixedrange: true,
              dtick: 1,
              zeroline: false,
              showgrid: settings.showGrid,
              showline: settings.showAxisLines,
              showticklabels: settings.showTickLabels,
            },
            dragmode: "pan",
            barmode: "overlay",
            showlegend: false,
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
      <RaincloudSidebar
        settings={settings}
        updateSettings={updateSettings}
        sidebarComponent={sidebarComponent}
      />
    </>
  );
};
