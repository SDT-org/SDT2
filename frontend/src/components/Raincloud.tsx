import Plotly from "plotly.js-dist-min";
import type { Layout, PlotData } from "plotly.js-dist-min";
import React from "react";
import { Label, ToggleButton, ToggleButtonGroup } from "react-aria-components";
import createPlotlyComponent from "react-plotly.js/factory";
import type { ColorString } from "../colors";
import { plotFont } from "../constants";
import type { DataSets, DistributionState } from "../distributionState";
import { arrayMinMax } from "../helpers";
import type { DistributionData, MetaData } from "../plotTypes";
import { ColorPicker } from "./ColorPicker";
import { Select, SelectItem } from "./Select";
import { Slider } from "./Slider";
import { Tooltip } from "./Tooltip";

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
        hovertemplate: "%{text}<br> <br>Percent Identity: %{x}<extra></extra>",
        text: data.identity_combos.map(
          (ids) => `Seq 1: ${ids[0]}<br>Seq 2: ${ids[1]}`,
        ),
      }) as Partial<PlotData>,
    [data.identity_combos, dataSet, settings],
  );

  const layout = React.useMemo(() => {
    return {
      title: "",
      uirevision: "true",
      font: plotFont,
      xaxis: {
        side: "left",
        rangemode: "tozero",
        fixedrange: true,
        zeroline: false,
        dtick:5,
        showgrid: settings.showGrid,
        showline: settings.showAxisLines,
        showticklabels: settings.showTickLabels,
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
  }, [settings, minDataValue, maxDataValue]);

  return (
    <>
      <div className="app-sidebar">
        <div className="app-sidebar-toolbar">
          <div className="form">
            {sidebarComponent}

            <div className="group">
              <div className="drawer">
                <ToggleButtonGroup
                  data-icon-only
                  selectionMode="multiple"
                  selectedKeys={Object.keys(settings).filter(
                    (key) =>
                      [
                        "showGrid",
                        "showTickLabels",
                        "showAxisLines",
                        "showAxisLabels",
                        "makeEditable",
                        "showMeanline",
                      ].includes(key) && settings[key as keyof typeof settings],
                  )}
                  onSelectionChange={(value) =>
                    updateSettings({
                      showGrid: value.has("showGrid"),
                      showTickLabels: value.has("showTickLabels"),
                      showAxisLines: value.has("showAxisLines"),
                      showAxisLabels: value.has("showAxisLabels"),
                      makeEditable: value.has("makeEditable"),
                      showMeanline: value.has("showMeanline"),
                    })
                  }
                >
                  <Tooltip tooltip="Toggle grid">
                    <ToggleButton id="showGrid" aria-label="Toggle grid">
                      <svg
                        xmlns="http://www.w3.org/2000/svg"
                        viewBox="0 0 24 24"
                        aria-hidden="true"
                      >
                        <g
                          style={{
                            fill: "none",
                            stroke: "currentcolor",
                            strokeWidth: 2,
                            strokeLinecap: "round",
                            strokeLinejoin: "round",
                            strokeMiterlimit: 10,
                          }}
                        >
                          <path d="M9 3v18M15 3v18M3 15h18M21 9H3M19 21H5a2 2 0 0 1-2-2V5a2 2 0 0 1 2-2h14a2 2 0 0 1 2 2v14a2 2 0 0 1-2 2z" />
                        </g>
                      </svg>
                    </ToggleButton>
                  </Tooltip>
                  <Tooltip tooltip="Toggle axis lines">
                    <ToggleButton
                      id="showAxisLines"
                      aria-label="Toggle axis lines"
                    >
                      <svg
                        xmlns="http://www.w3.org/2000/svg"
                        viewBox="0 0 24 24"
                        aria-hidden="true"
                      >
                        <g
                          style={{
                            fill: "none",
                            stroke: "currentcolor",
                            strokeWidth: 2,
                            strokeLinecap: "round",
                            strokeLinejoin: "round",
                            strokeMiterlimit: 10,
                          }}
                        >
                          <path d="M3 5a2 2 0 0 1 2 -2h14a2 2 0 0 1 2 2v14a2 2 0 0 1 -2 2h-14a2 2 0 0 1 -2 -2v-14z" />
                          <path d="M3 10h18" />
                          <path d="M10 3v18" />
                        </g>
                      </svg>
                    </ToggleButton>
                  </Tooltip>
                  <Tooltip tooltip="Toggle axis values">
                    <ToggleButton
                      id="showTickLabels"
                      aria-label="Toggle axis values"
                    >
                      <svg
                        xmlns="http://www.w3.org/2000/svg"
                        viewBox="0 0 24 24"
                        aria-hidden="true"
                      >
                        <g
                          style={{
                            fill: "none",
                            stroke: "currentcolor",
                            strokeWidth: 2,
                            strokeLinecap: "round",
                            strokeLinejoin: "round",
                            strokeMiterlimit: 10,
                          }}
                        >
                          <path d="M18 22H6a4 4 0 0 1-4-4V6a4 4 0 0 1 4-4h12a4 4 0 0 1 4 4v12a4 4 0 0 1-4 4zM9 7v10M6 9l3-2" />
                          <path d="M15.5 17a2.5 2.5 0 0 1-2.5-2.5v-5a2.5 2.5 0 1 1 5 0v5a2.5 2.5 0 0 1-2.5 2.5z" />
                        </g>
                      </svg>
                    </ToggleButton>
                  </Tooltip>
                  <Tooltip tooltip="Toggle editable mode">
                    <ToggleButton
                      id="makeEditable"
                      aria-label="Toggle editable mode"
                    >
                      <svg
                        xmlns="http://www.w3.org/2000/svg"
                        viewBox="0 0 24 24"
                        aria-hidden="true"
                      >
                        <g
                          style={{
                            fill: "none",
                            stroke: "currentcolor",
                            strokeWidth: 2,
                            strokeLinecap: "round",
                            strokeLinejoin: "round",
                            strokeMiterlimit: 10,
                          }}
                        >
                          <path d="M14 2 L18 6 L7 17 H3 V13 Z" />
                          <path d="M3 22 L21 22" />
                        </g>
                      </svg>
                    </ToggleButton>
                  </Tooltip>
                  <Tooltip tooltip="Toggle mean line">
                    <ToggleButton
                      id="showMeanline"
                      aria-label="Toggle mean line"
                    >
                      <svg
                        xmlns="http://www.w3.org/2000/svg"
                        viewBox="0 0 24 24"
                        aria-hidden="true"
                      >
                        <g
                          style={{
                            fill: "none",
                            stroke: "currentcolor",
                            strokeWidth: 2,
                            strokeLinecap: "round",
                            strokeLinejoin: "round",
                            strokeMiterlimit: 10,
                          }}
                        >
                          <path d="M5 12h2" />
                          <path d="M17 12h2" />
                          <path d="M11 12h2" />
                        </g>
                      </svg>
                    </ToggleButton>
                  </Tooltip>
                </ToggleButtonGroup>
                <hr className="compact" />
                <div className="col-2 auto-onefr align-items-center color-slider-gap">
                  <ColorPicker
                    value={settings.fillColor}
                    onChange={(value) => {
                      updateSettings({
                        fillColor: value.toString() as ColorString,
                      });
                    }}
                  />
                  <Slider
                    label="Bandwidth"
                    defaultValue={settings.bandwidth}
                    onChange={(value) => updateSettings({ bandwidth: value })}
                    minValue={1}
                    maxValue={20}
                    step={1}
                  />
                </div>
                <hr className="compact" />
                <div className="col-2 auto-onefr align-items-center color-slider-gap">
                  <ColorPicker
                    value={settings.lineColor}
                    onChange={(value) =>
                      updateSettings({
                        lineColor: value.toString() as ColorString,
                      })
                    }
                  />
                  <Slider
                    label="Width"
                    value={settings.lineWidth}
                    onChange={(value) => updateSettings({ lineWidth: value })}
                    minValue={1}
                    maxValue={20}
                    step={1}
                  />
                </div>
              </div>
            </div>
            <div className="group">
              <h4 className="setting-header">Points</h4>
              <div
                className="drawer"
                data-hidden={!settings.showPoints}
                aria-hidden={!settings.showPoints}
              >
                <div className="col-2 auto-onefr align-items-center">
                  <Label htmlFor="points-type">Type</Label>
                  <Select
                    data-compact
                    id="points-type"
                    selectedKey={settings.points}
                    onSelectionChange={(value) => {
                      updateSettings({
                        points: value as typeof settings.points,
                      });
                    }}
                    items={Object.entries({
                      all: "All",
                      outliers: "Outliers",
                      suspectedoutliers: "Suspected Outliers",
                    }).map(([id, name]) => ({
                      id,
                      name,
                    }))}
                  >
                    {(item) => (
                      <SelectItem textValue={item.name}>{item.name}</SelectItem>
                    )}
                  </Select>
                </div>

                <Slider
                  label="Size"
                  value={settings.markerSize}
                  onChange={(value) => updateSettings({ markerSize: value })}
                  minValue={1}
                  maxValue={20}
                  step={1}
                />

                <Slider
                  label="Position"
                  value={settings.pointPos}
                  onChange={(value) => updateSettings({ pointPos: value })}
                  minValue={-2}
                  maxValue={2}
                  step={0.1}
                />

                <Slider
                  label="Jitter"
                  value={settings.jitter}
                  onChange={(value) => updateSettings({ jitter: value })}
                  minValue={0.1}
                  maxValue={1}
                  step={0.1}
                />

                <div className="col-2 onefr-auto small-color align-items-center">
                  <Label>Color</Label>
                  <ColorPicker
                    value={settings.markerColor}
                    onChange={(value) => {
                      updateSettings({
                        markerColor: value.toString() as ColorString,
                      });
                    }}
                  />
                </div>
              </div>
            </div>
          </div>
        </div>
        {footer ? <div className="app-sidebar-footer">{footer}</div> : null}
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
            editable: !!settings.makeEditable,
          }}
          style={{ width: "100%", height: "100%" }}
        />
      </div>
    </>
  );
};
