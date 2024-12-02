import Plotly from "plotly.js-dist-min";
import type { Layout, PlotData } from "plotly.js-dist-min";
import React from "react";
import { Label, ToggleButton, ToggleButtonGroup } from "react-aria-components";
import createPlotlyComponent from "react-plotly.js/factory";
import type { ColorString } from "../colors";
import type { DataSets, DistributionState } from "../distributionState";
import type { DistributionData } from "../plotTypes";
import { ColorPicker } from "./ColorPicker";
import { Slider } from "./Slider";
import { Switch } from "./Switch";
import { Tooltip } from "./Tooltip";

const Plot = createPlotlyComponent(Plotly);

export const Histogram = ({
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

  const histogramTrace = React.useMemo(
    () =>
      ({
        type: "histogram",
        x: dataSet,
        histnorm: "percent",
        marker: {
          color: settings.binColor,
          line: {
            width: settings.histOutlineWidth,
            color: settings.histlineColor,
          },
        },
        xbins: {
          size: settings.binSize,
        },
        name: "Histogram",
        hovertemplate:
          "Percent Identity: %{x}<br>Percentage: %{y}<extra></extra>",
      }) as Partial<PlotData>,
    [settings, dataSet],
  );

  const layout = React.useMemo(
    () =>
      ({
        title: "",
        uirevision: "true",
        xaxis: {
          title: settings.showAxisLabels ? "Percent Pairwise Identity" : "",
          side: "bottom",
          rangemode: "normal",
          fixedrange: true,
          zeroline: false,
          showgrid: settings.showGrid,
          showticklabels: settings.showTickLabels,
          showline: settings.showLine,
        },
        yaxis: {
          title: settings.showAxisLabels
            ? "Proportion of Pairwise Identities"
            : "",
          side: "left",
          rangemode: "tozero",
          fixedrange: true,
          zeroline: false,
          showgrid: settings.showGrid,
          showticklabels: settings.showTickLabels,
          showline: settings.showLine,
        },
        dragmode: "pan",
        barmode: "overlay",
        showlegend: true,
        margin: { l: 50, r: 50, t: 50, b: 50 },
      }) as Partial<Layout>,
    [settings],
  );

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
                        "showLine",
                        "showAxisLabels",
                      ].includes(key) && settings[key as keyof typeof settings],
                  )}
                  onSelectionChange={(value) =>
                    updateSettings({
                      showGrid: value.has("showGrid"),
                      showTickLabels: value.has("showTickLabels"),
                      showLine: value.has("showLine"),
                      showAxisLabels: value.has("showAxisLabels"),
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
                    <ToggleButton id="showLine" aria-label="Toggle axis lines">
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
                          <path d="M1 1v22h22" />
                          <path d="m7 17 5-6 5 1 6-7" />
                        </g>
                      </svg>
                    </ToggleButton>
                  </Tooltip>
                  <Tooltip tooltip="Toggle axis labels">
                    <ToggleButton
                      id="showAxisLabels"
                      aria-label="Toggle axis labels"
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
                          <path d="M18 22H6a4 4 0 0 1-4-4V6a4 4 0 0 1 4-4h12a4 4 0 0 1 4 4v12a4 4 0 0 1-4 4z" />
                          <path d="M8 17v-6a4 4 0 0 1 8 0v6M8 13h8" />
                        </g>
                      </svg>
                    </ToggleButton>
                  </Tooltip>
                  <Tooltip tooltip="Toggle axis title">
                    <ToggleButton
                      id="showTickLabels"
                      aria-label="Toggle axis title"
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
                </ToggleButtonGroup>
              </div>
            </div>
            <div className="group">
              <div className="row">
                <div className="col-2">
                  <Label className="header">Bins</Label>
                  <div className="auto-onefr">
                    <ColorPicker
                      value={settings.binColor}
                      onChange={(value) => {
                        updateSettings({
                          binColor: value.toString() as ColorString,
                        });
                      }}
                    />
                    <Slider
                      label="Width"
                      defaultValue={settings.binSize}
                      onChange={(value) => updateSettings({ binSize: value })}
                      minValue={0}
                      maxValue={5}
                      step={0.5}
                    />
                  </div>
                </div>
                <div className="col-2">
                  <Label className="header">Outline</Label>
                  <div className="auto-onefr">
                    <ColorPicker
                      value={settings.histlineColor}
                      onChange={(value) => {
                        updateSettings({
                          histlineColor: value.toString() as ColorString,
                        });
                      }}
                    />
                    <Slider
                      label="Width"
                      defaultValue={settings.histOutlineWidth}
                      onChange={(value) =>
                        updateSettings({ histOutlineWidth: value })
                      }
                      minValue={0}
                      maxValue={15}
                      step={1}
                    />
                  </div>
                </div>
              </div>
            </div>
          </div>
        </div>
        {footer ? <div className="app-sidebar-footer">{footer}</div> : null}
      </div>
      <div className="app-main">
        <Plot
          data={[histogramTrace]}
          layout={layout}
          config={{
            responsive: true,
            displayModeBar: false,
            scrollZoom: true,
            displaylogo: false,
            editable: true,
          }}
          style={{ width: "100%", height: "100%" }}
        />
      </div>
    </>
  );
};
