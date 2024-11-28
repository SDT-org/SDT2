import React from "react";
import Plotly from "plotly.js-dist-min";
import createPlotlyComponent from "react-plotly.js/factory";
import { NumberInput } from "./NumberInput";
import { Layout, PlotData } from "plotly.js-dist-min";
import { DistributionData } from "../plotTypes";
import { formatTitle } from "../helpers";
import { ColorOptions } from "./ColorOptions";
import { DataSets, DistributionState } from "../distributionState";
import { ColorOption } from "../colors";

const Plot = createPlotlyComponent(Plotly);

export const Raincloud = ({
  data,
  dataSets,
  dataSetKey,
  footer,
  sidebarComponent,
  settings,
  updateSettings,
}: {
  data: DistributionData | undefined;
  dataSets: DataSets;
  dataSetKey: keyof DataSets;
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
  const minDataValue = Math.min(...dataSet);
  const maxDataValue = Math.max(...dataSet);

  const rainCloudTrace = React.useMemo(
    () =>
      ({
        type: "violin",
        name: "",
        x: dataSet,
        side: "negative",
        points: settings.points !== "None" ? settings.points : false,
        line: {
          color: settings.lineColor,
          width: settings.lineWidth,
        },
        fillcolor: settings.fillColor,
        opacity: settings.violinOpacity,
        pointpos: settings.pointPos,
        jitter: settings.jitter,
        bandwidth: settings.bandwidth,
        marker: {
          visible: true,
          color: settings.markerColor,
          size: settings.markerSize,
          opacity: settings.pointOpacity,
        },
        meanline: {
          visible: true,
        },
        hovertemplate: "%{text}<br> <br>Percent Identity: %{x}<extra></extra>",
        text: data.identity_combos.map(
          (ids) => `Seq 1: ${ids[0]}<br>Seq 2: ${ids[1]}`,
        ),
      }) as Partial<PlotData>,
    [data, dataSetKey, settings],
  );

  const layout = React.useMemo(() => {
    return {
      title: settings.plotTitle,
      xaxis: {
        side: "left",
        rangemode: "tozero",
        fixedrange: true,
        zeroline: false,
        showgrid: settings.showGrid,
        showline: settings.showAxisLines,
        showticklabels: settings.showTickLabels,
        title: settings.showAxisLabels
          ? "Percent Pairwise Identity"
          : undefined,
        range: [minDataValue - 20, maxDataValue + 20],
      },
      yaxis: {
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
    } as Partial<Layout>;
  }, [data, dataSetKey, settings]);

  return (
    <>
      <div className="app-sidebar">
        <div className="app-sidebar-toolbar">
          <div className="form">
            {sidebarComponent}
            <div className="group">
              <div className="field">
                <label className="header">Title</label>
                <input
                  type="text"
                  value={settings.plotTitle}
                  onChange={(e) =>
                    updateSettings({ plotTitle: e.target.value })
                  }
                />
              </div>
              <div className="row">
                <div className="col-2">
                  <div className="field">
                    <label htmlFor="showGrid">
                      <input
                        type="checkbox"
                        name="showGrid"
                        id="showGrid"
                        checked={settings.showGrid}
                        onChange={() =>
                          updateSettings({ showGrid: !settings.showGrid })
                        }
                      />
                      Grid
                    </label>
                  </div>
                  <div className="field">
                    <label htmlFor="showTickLabels">
                      <input
                        type="checkbox"
                        name="showTickLabels"
                        id="showTickLabels"
                        checked={settings.showTickLabels}
                        onChange={() =>
                          updateSettings({
                            showTickLabels: !settings.showTickLabels,
                          })
                        }
                      />
                      Tick Labels
                    </label>
                  </div>
                </div>
                <div className="col-2">
                  <div className="field">
                    <label htmlFor="showAxisLines">
                      <input
                        type="checkbox"
                        name="showAxisLines"
                        id="showAxisLines"
                        checked={settings.showAxisLines}
                        onChange={() =>
                          updateSettings({
                            showAxisLines: !settings.showAxisLines,
                          })
                        }
                      />
                      Axis Lines
                    </label>
                  </div>
                  <div className="field">
                    <label htmlFor="showAxisLabels">
                      <input
                        type="checkbox"
                        name="showAxisLabels"
                        id="showAxisLabels"
                        checked={settings.showAxisLabels}
                        onChange={() =>
                          updateSettings({
                            showAxisLabels: !settings.showAxisLabels,
                          })
                        }
                      />
                      Axis Title
                    </label>
                  </div>
                </div>
              </div>
            </div>
            <div className="group">
              <div className="row">
                <div className="col-2">
                  <div className="field">
                    <label htmlFor="fill-color">Cloud Fill Color</label>
                    <select
                      id="fill-color"
                      value={settings.fillColor}
                      onChange={(e) =>
                        updateSettings({
                          fillColor: e.target.value as ColorOption,
                        })
                      }
                    >
                      <ColorOptions />
                    </select>
                  </div>
                  <div className="field">
                    <NumberInput
                      label="Band Width"
                      field="bandWidth"
                      value={settings.bandwidth}
                      updateValue={updateSettings}
                      min={0}
                      max={20}
                      step={1}
                    />
                  </div>
                </div>
                <div className="col-2">
                  <div className="field">
                    <label htmlFor="line-color">Line Color</label>
                    <select
                      id="line-color"
                      value={settings.lineColor}
                      onChange={(e) =>
                        updateSettings({
                          lineColor: e.target.value as ColorOption,
                        })
                      }
                    >
                      <ColorOptions />
                    </select>
                  </div>
                  <div className="field">
                    <NumberInput
                      label="Line Width"
                      field="lineWidth"
                      value={settings.lineWidth}
                      updateValue={updateSettings}
                      min={0}
                      max={20}
                      step={1}
                    />
                  </div>
                </div>
              </div>
            </div>
            <div className="group">
              <div className="row">
                <div className="col-2">
                  <div className="field">
                    <NumberInput
                      label="Point Position"
                      field="pointPos"
                      value={settings.pointPos}
                      type="float"
                      updateValue={updateSettings}
                      min={-2}
                      max={-1}
                      step={0.1}
                    />
                  </div>
                  <div className="field">
                    <label htmlFor="points">Points</label>
                    <select
                      id="points"
                      value={settings.points}
                      onChange={(e) =>
                        updateSettings({
                          points: e.target.value as
                            | "all"
                            | "outliers"
                            | "suspectedoutliers"
                            | "None",
                        })
                      }
                    >
                      {["all", "outliers", "suspectedoutliers", "None"].map(
                        (value) => (
                          <option key={value} value={value}>
                            {value === "None" ? "None" : formatTitle(value)}
                          </option>
                        ),
                      )}
                    </select>
                  </div>
                </div>
                <div className="row">
                  <div className="col-2">
                    <div className="field">
                      <label htmlFor="markerColor">Point Color</label>
                      <select
                        id="markerColor"
                        value={settings.markerColor}
                        onChange={(e) =>
                          updateSettings({
                            markerColor: e.target.value as ColorOption,
                          })
                        }
                      >
                        <ColorOptions />
                      </select>
                    </div>
                    <div className="field">
                      <NumberInput
                        label="Point Size"
                        field="markerSize"
                        value={settings.markerSize}
                        updateValue={updateSettings}
                        min={0}
                        max={20}
                        step={1}
                      />
                    </div>
                  </div>
                  <div className="row">
                    <div className="col-2">
                      <div className="field">
                        <NumberInput
                          label="Point Opacity"
                          field="pointOpacity"
                          type="float"
                          value={settings.pointOpacity}
                          updateValue={updateSettings}
                          min={0}
                          max={1}
                          step={0.1}
                        />
                      </div>
                      <div className="field">
                        <NumberInput
                          label="Jitter"
                          field="jitter"
                          value={settings.jitter}
                          type="float"
                          updateValue={updateSettings}
                          min={0}
                          max={1}
                          step={0.1}
                        />
                      </div>
                    </div>
                  </div>
                </div>
              </div>
            </div>
          </div>
          <div className="app-sidebar-footer">{footer}</div>
        </div>
      </div>
      <div className="app-main">
        <Plot
          data={[rainCloudTrace]}
          layout={layout}
          config={{
            responsive: true,
            displayModeBar: false,
            scrollZoom: true,
            displaylogo: false,
          }}
          style={{ width: "100%", height: "100%" }}
        />
      </div>
    </>
  );
};
