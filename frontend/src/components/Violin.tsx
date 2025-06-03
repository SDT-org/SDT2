import type { PlotData } from "plotly.js";
import Plotly from "plotly.js-cartesian-dist-min";
import React from "react";
import createPlotlyComponent from "react-plotly.js/factory";
import { plotFontMonospace, plotFontSansSerif } from "../constants";
import type { DataSets, DistributionState } from "../distributionState";
import { arrayMinMax } from "../helpers";
import { useExportPlotly } from "../hooks/useExportPlotly";
import {
  useRelayoutHideSubtitle,
  useRelayoutUpdateTitles,
} from "../hooks/useRelayoutUpdateTitles";
import type { DistributionData, MetaData } from "../plotTypes";
import { ViolinSidebar } from "./ViolinSidebar";

const Plot = createPlotlyComponent(Plotly);

const getAxisConfig = (
  minDataValue: number,
  maxDataValue: number,
  isLength: boolean,
) => {
  if (isLength) {
    return {
      tickmode: "linear" as const,
      tickfont: {
        family: plotFontMonospace.family,
        size: 11,
      },
      dtick: 10,
      range: [minDataValue - 20, maxDataValue + 20],
    };
  }

  return {
    tickmode: "array" as const,
    tickfont: {
      family: plotFontMonospace.family,
      size: 11,
    },
    range: [minDataValue - 20, maxDataValue + 20],
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

  const exportPlotly = useExportPlotly("distribution_violin");

  const axisTitle = React.useMemo(() => {
    if (dataSetKey === "scores") {
      return "Percent Pairwise Identity";
    }
    if (dataSetKey === "gc") {
      return "GC Percentage";
    }
    if (dataSetKey === "length") {
      return "Length (nt)";
    }
    return "";
  }, [dataSetKey]);

  const hoverData = {
    scores: {
      index: 0,
      suffix: "%",
      title: "Identity Score",
    },
    gc: {
      index: 1,
      suffix: "%",
      title: "GC Content",
    },
    length: {
      index: 2,
      suffix: "nt",
      title: "Sequence Length",
    },
  };

  const scoresText = data.identity_combos.map(
    (idIndexes) =>
      `Seq 1: ${data.ids[idIndexes[0]]}<br>Seq 2: ${data.ids[idIndexes[1]]}<br><br>${hoverData.scores.title}: `,
  );
  const { index, suffix, title } = hoverData[dataSetKey];

  const gcLengthText = data.full_stats.map(
    (stats) => `${stats[0]}:<br>${title} = ${stats[index]}${suffix}`,
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
          width: settings.showBox ? settings.boxlineWidth : settings.lineWidth,
          color: settings.lineColor,
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
        // we are going to rip this nonsense out for D3 ASAP
        hovertemplate:
          dataSetKey === "scores"
            ? `%{text}%{${settings.plotOrientation === "vertical" ? "y" : "x"}}${hoverData.scores.suffix}<extra></extra>`
            : "%{text}<extra></extra>",
        text: hoverText,
      }) as Partial<PlotData>,
    [dataSet, dataSetKey, settings, hoverText],
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
          visible: settings.showPoints,
          color: settings.markerColor,
          size: settings.markerSize,
        },
        fillcolor: settings.boxfillColor,
        width: settings.boxWidth,
        hoverinfo: "text",
        hovertemplate:
          dataSetKey === "scores"
            ? `%{text}%{${settings.plotOrientation === "vertical" ? "y" : "x"}}${hoverData.scores.suffix}<extra></extra>`
            : "%{text}<extra></extra>",
        text: hoverText,
        hoveron: "points",
      }) as Partial<PlotData>,
    [dataSet, dataSetKey, settings, hoverText],
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
                  },
                }
              : {}),
            font: {
              ...(settings.titleFont === "Monospace"
                ? plotFontMonospace
                : plotFontSansSerif),
              size: 14,
            },
            uirevision: settings.plotOrientation,
            xaxis: {
              title: {
                text: settings.showAxisLabels
                  ? settings.plotOrientation === "vertical"
                    ? "Distribution"
                    : axisTitle
                  : "",
                font: { family: "sans-serif", size: 12 },
              },
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
                    ...getAxisConfig(
                      minDataValue,
                      maxDataValue,
                      dataSetKey === "length",
                    ),
                  }),
            },
            yaxis: {
              title: {
                text: settings.showAxisLabels
                  ? settings.plotOrientation === "horizontal"
                    ? "Distribution"
                    : axisTitle
                  : "",
                font: { family: "sans-serif", size: 12 },
              },
              ...(settings.plotOrientation === "vertical"
                ? {
                    side: "left",
                    fixedrange: true,
                    zeroline: false,
                    showgrid: settings.showGrid,
                    showline: settings.showAxisLines,
                    showticklabels: settings.showTickLabels,
                    ...getAxisConfig(
                      minDataValue,
                      maxDataValue,
                      dataSetKey === "length",
                    ),
                  }
                : {
                    fixedrange: true,
                    zeroline: false,
                    showgrid: settings.showGrid,
                    showline: settings.showAxisLines,
                    showticklabels: settings.showTickLabels,
                  }),
            },
            dragmode: "pan",
            barmode: "overlay",
            showlegend: false,
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
      <ViolinSidebar
        settings={settings}
        updateSettings={updateSettings}
        sidebarComponent={sidebarComponent}
        dataSetLength={dataSet.length}
      />
    </>
  );
};
