import type { PlotData } from "plotly.js";
import Plotly from "plotly.js-cartesian-dist-min";
import React from "react";
import createPlotlyComponent from "react-plotly.js/factory";
import { plotFontMonospace, plotFontSansSerif } from "../constants";
import type { DataSets, DistributionState } from "../distributionState";
import { arrayMinMax } from "../helpers";
import { useRelayoutHideSubtitle, useRelayoutUpdateTitles } from "../hooks/useRelayoutUpdateTitles";
import type { DistributionData, MetaData } from "../plotTypes";
import { ViolinSidebar } from "./ViolinSidebar";

const Plot = createPlotlyComponent(Plotly);

const getAxisConfig = (minDataValue: number, maxDataValue: number, isLength: boolean) => {
  if (isLength) {
    return {
      tickmode: 'linear',
      tickfont: {
        ...plotFontMonospace,
      },
      dtick: 10,
      range: [minDataValue - 20, maxDataValue + 20]
    };
  }
  
  return {
    tickmode: 'array',
    tickfont: {
      ...plotFontMonospace,
    },
    range: [minDataValue - 20, maxDataValue + 20],
    ticktext: Array.from(
      {length: Math.ceil((maxDataValue + 20) / 5) + 1},
      (_, i) => {
        const value = i * 5;
        return value <= 100 ? value.toString() : '';
      }
    ),
    tickvals: Array.from(
      {length: Math.ceil((maxDataValue + 20) / 5) + 1},
      (_, i) => i * 5
    )
  };
};

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
  const [minDataValue, maxDataValue] = arrayMinMax(dataSet);
  const updateTitles = useRelayoutUpdateTitles(updateSettings);
  useRelayoutHideSubtitle(!settings.showTitles);

  const scoresText = data.identity_combos.map(
    (ids) =>
      `Seq 1: ${ids[0]}<br>Seq 2: ${ids[1]} <br>Percent Identity: %{${settings.plotOrientation === "vertical" ? "y" : "x"}}`,
  );

  const gcLengthIndex = dataSetKey === "gc" ? 1 : 2;
  const gcLengthSuffix = dataSetKey === "gc" ? "%" : " nt";
  const gcLengthTitle = dataSetKey === "gc" ? "GC" : "Length";
  const gcLengthText = data.full_stats.map(
    (stats) =>
      `${stats[0]}: ${gcLengthTitle} = ${stats[gcLengthIndex]}${gcLengthSuffix}`,
  );

  const hoverText = dataSetKey === "scores" ? scoresText : gcLengthText;

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
        hovertemplate: "%{text}<extra></extra>",
        text: hoverText,
      }) as Partial<PlotData>,
    [dataSet, settings, hoverText],
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
                    text: settings.title,
                    pad: {
                      t: 0,
                      r: 0,
                      b: 0,
                      l: 0,
                    },
                    subtitle: {
                      text: settings.subtitle,
                    },
                  },
                }
              : {}),
            font: {
              ...(settings.titleFont === "Monospace"
                ? plotFontMonospace
                : plotFontSansSerif),
              weight: "bold",
            },
            uirevision: settings.plotOrientation,
            boxgap: 1 - settings.boxWidth,
            xaxis: {
              ...(settings.plotOrientation === "vertical"
                ? {
                    fixedrange: true,
                    zeroline: false,
                    showgrid: settings.showGrid,
                    showline: settings.showAxisLines,
                    showticklabels: settings.showTickLabels,
                  }
                : {
                    side: "left",
                    fixedrange: true,
                    zeroline: false,
                    showgrid: settings.showGrid,
                    showline: settings.showAxisLines,
                    showticklabels: settings.showTickLabels,
                    ...getAxisConfig(minDataValue, maxDataValue, dataSetKey === "length")
                  }
              ),
              ...(settings.showTitles
                ? {
                    title: {
                      text: settings.xtitle,
                      font: {
                        ...(settings.titleFont === "Monospace"
                          ? plotFontMonospace
                          : plotFontSansSerif),
                        weight: "bold",
                      },
                    },
                    scaleratio: 1,
                  }
                : {}),
            },
            yaxis: {
              ...(settings.plotOrientation === "vertical"
                ? {
                    side: "left",
                    fixedrange: true,
                    zeroline: false,
                    showgrid: settings.showGrid,
                    showline: settings.showAxisLines,
                    showticklabels: settings.showTickLabels,
                    ...getAxisConfig(minDataValue, maxDataValue, dataSetKey === "length")
                  }
                : {
                    fixedrange: true,
                    zeroline: false,
                    showgrid: settings.showGrid,
                    showline: settings.showAxisLines,
                    showticklabels: settings.showTickLabels,
                  }),
              ...(settings.showTitles
                ? {
                    title: {
                      text: settings.ytitle,
                      font: {
                        ...(settings.titleFont === "Monospace"
                          ? plotFontMonospace
                          : plotFontSansSerif),
                        weight: "bold",
                      },
                      pad: {
                        r: 15,
                      },
                    },
                  }
                : {}),
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