import type { PlotData } from "plotly.js";
import Plotly from "plotly.js-cartesian-dist-min";
import React from "react";
import createPlotlyComponent from "react-plotly.js/factory";
import { plotFont } from "../constants";
import type { DataSets, DistributionState } from "../distributionState";
import { arrayMinMax } from "../helpers";
import {
  useRelayoutHideSubtitle,
  useRelayoutUpdateTitles,
} from "../hooks/useRelayoutUpdateTitles";
import type { DistributionData, MetaData } from "../plotTypes";
import { ViolinSidebar } from "./ViolinSidebar";

const Plot = createPlotlyComponent(Plotly);

export const Violin = ({
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
  metaData?: MetaData;
  footer?: React.ReactNode;
  sidebarComponent?: React.ReactNode;
  settings: DistributionState["violin"];
  updateSettings: React.Dispatch<Partial<DistributionState["violin"]>>;
}) => {
  if (!data) {
    return (
      <div className="app-main centered">
        <h1>No data available</h1>
      </div>
    );
  }

  const dataSet = dataSets[dataSetKey];
  console.log("dataSetKey:", dataSetKey);
  console.log("dataSet:", dataSet);
  console.log("data.identity_combos:", data.identity_combos);
  const [minDataValue, maxDataValue] = arrayMinMax(dataSet);
  const updateTitles = useRelayoutUpdateTitles(updateSettings);
  useRelayoutHideSubtitle(!settings.showTitles);

  const violinTrace = React.useMemo(
    () =>
      ({
        type: "violin",
        name: "",
        [settings.plotOrientation === "vertical" ? "y" : "x"]: dataSet,
        visible: settings.showViolin,
        line: {
          color: settings.lineColor,
          width: settings.lineWidth,
        },
        fillcolor: settings.fillColor,
        meanline: {
          visible: settings.showMeanline,
        },
        points:
          settings.showPoints &&
          settings.pointOrientation === "Violin" &&
          settings.points,
        pointpos: settings.pointPos,
        jitter: settings.jitter,
        bandwidth: settings.bandwidth,
        marker: {
          visible: settings.showPoints,
          color: settings.markerColor,
          size: settings.markerSize,
        },
        hoveron: "points",
        scalemode: "width",
        hovertemplate: `%{text}<br>Percent Identity: %{${settings.plotOrientation === "vertical" ? "y" : "x"}}<extra></extra>`,
        text: data.identity_combos.map(
          (ids) => `Seq 1: ${ids[0]}<br>Seq 2: ${ids[1]}`,
        ),
      }) as Partial<PlotData>,
    [data, dataSet, settings],
  );
  const boxTrace = React.useMemo(
    () =>
      ({
        type: "box",
        name: "",
        [settings.plotOrientation === "vertical" ? "y" : "x"]: dataSet,
        boxmean: settings.showMeanline,
        boxpoints:
          settings.showPoints &&
          settings.pointOrientation === "Box" &&
          settings.points,
        pointpos: settings.pointPos,
        jitter: settings.jitter,
        line: {
          color: settings.boxlineColor,
          width: settings.boxlineWidth,
        },
        whiskerwidth: settings.whiskerWidth,
        marker: {
          visible: settings.showBox,
          color: settings.markerColor,
          size: settings.markerSize,
        },
        fillcolor: settings.boxfillColor,
        hovermode: "closest",
        boxgap: 1 - settings.boxWidth,
        hovertemplate:
          "Percent Identity: %{x}<br>Percent Identity: %{y}<extra></extra>",
      }) as Partial<PlotData>,
    [dataSet, settings],
  );

  return (
    <>
      <div className="app-main">
        <Plot
          data={[
            settings.showViolin ? violinTrace : {},
            settings.showBox ? boxTrace : {},
          ]}
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
            uirevision: settings.plotOrientation,
            font: plotFont,
            // @ts-ignore
            boxgap: 1 - settings.boxWidth,
            xaxis:
              settings.plotOrientation === "vertical"
                ? {
                    fixedrange: true,
                    zeroline: false,
                    showgrid: settings.showGrid,
                    showline: settings.showAxisLines,
                    showticklabels: settings.showTickLabels,
                  }
                : {
                    side: "left",
                    rangemode: "tozero",
                    fixedrange: true,
                    zeroline: false,
                    showgrid: settings.showGrid,
                    showline: settings.showAxisLines,
                    showticklabels: settings.showTickLabels,
                    tickmode: "auto",
                    autotick: true,
                    range: [minDataValue - 20, maxDataValue + 20],
                  },
            yaxis:
              settings.plotOrientation === "vertical"
                ? {
                    side: "left",
                    rangemode: "tozero",
                    fixedrange: true,
                    zeroline: false,
                    showgrid: settings.showGrid,
                    showline: settings.showAxisLines,
                    showticklabels: settings.showTickLabels,
                    tickmode: "auto",
                    autotick: true,
                    range: [minDataValue - 20, maxDataValue + 20],
                  }
                : {
                    fixedrange: true,
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
      <ViolinSidebar
        settings={settings}
        updateSettings={updateSettings}
        sidebarComponent={sidebarComponent}
      />
    </>
  );
};
