import type { PlotData } from "plotly.js";
import Plotly from "plotly.js-cartesian-dist-min";
import React from "react";
import {
  Input,
  Label,
  TextField,
  ToggleButton,
  ToggleButtonGroup,
} from "react-aria-components";
import createPlotlyComponent from "react-plotly.js/factory";
import type { ColorString } from "../colors";
import { plotFont } from "../constants";
import type { DataSets, DistributionState } from "../distributionState";
import {
  useRelayoutHideSubtitle,
  useRelayoutUpdateTitles,
} from "../hooks/useRelayoutUpdateTitles";
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
            font: plotFont,
            uirevision: "true",
            xaxis: {
              ...(settings.showTitles
                ? {
                    title: {
                      text: settings.xtitle,
                    },
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
      <div className="app-sidebar app-sidebar-right">
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
                        "showMeanline",
                      ].includes(key) && settings[key as keyof typeof settings],
                  )}
                  onSelectionChange={(value) =>
                    updateSettings({
                      showGrid: value.has("showGrid"),
                      showTickLabels: value.has("showTickLabels"),
                      showAxisLines: value.has("showAxisLines"),
                      showAxisLabels: value.has("showAxisLabels"),
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
                  <Tooltip tooltip="Toggle tick values">
                    <ToggleButton
                      id="showTickLabels"
                      aria-label="Toggle axis tick values"
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
                <Label>Orientation</Label>
                <ToggleButtonGroup
                  data-compact
                  selectionMode="single"
                  disallowEmptySelection={true}
                  selectedKeys={[settings.plotOrientation]}
                  onSelectionChange={(value) =>
                    updateSettings({
                      plotOrientation: value.values().next()
                        .value as typeof settings.plotOrientation,
                    })
                  }
                >
                  <ToggleButton id="vertical">Vertical</ToggleButton>
                  <ToggleButton id="horizontal">Horizontal</ToggleButton>
                </ToggleButtonGroup>
                <hr className="compact" />
                <div className="col-2 auto-onefr align-items-center color-slider-gap">
                  <ColorPicker
                    value={settings.binColor}
                    onChange={(value) => {
                      updateSettings({
                        binColor: value.toString() as ColorString,
                      });
                    }}
                  />
                  <Slider
                    label="Bin Width"
                    defaultValue={settings.binSize}
                    onChange={(value) => updateSettings({ binSize: value })}
                    minValue={0.5}
                    maxValue={5}
                    step={0.5}
                  />
                </div>
                <hr className="compact" />
                <div className="col-2 auto-onefr align-items-center color-slider-gap">
                  <ColorPicker
                    value={settings.histlineColor}
                    onChange={(value) => {
                      updateSettings({
                        histlineColor: value.toString() as ColorString,
                      });
                    }}
                  />
                  <Slider
                    label="Outline Width"
                    defaultValue={settings.histOutlineWidth}
                    onChange={(value) =>
                      updateSettings({ histOutlineWidth: value })
                    }
                    minValue={1}
                    maxValue={15}
                    step={1}
                  />
                </div>
              </div>
            </div>
            <div className="group">
              <Switch
                isSelected={settings.showTitles}
                onChange={(value) => {
                  updateSettings({
                    showTitles: value,
                  });
                }}
              >
                Plot Titles
              </Switch>
              <div
                className="drawer"
                data-hidden={!settings.showTitles}
                aria-hidden={!settings.showTitles}
              >
                <div className="field">
                  <TextField
                    onChange={(value) => updateSettings({ title: value })}
                    value={settings.title}
                  >
                    <Label>Title</Label>
                    <Input />
                  </TextField>
                </div>
                <div className="field">
                  <TextField
                    onChange={(value) => updateSettings({ subtitle: value })}
                    value={settings.subtitle}
                  >
                    <Label>Subtitle</Label>
                    <Input />
                  </TextField>
                </div>
                <div className="field">
                  <TextField
                    onChange={(value) => updateSettings({ xtitle: value })}
                    value={settings.xtitle}
                  >
                    <Label>X Axis Title</Label>
                    <Input />
                  </TextField>
                </div>
                <div className="field">
                  <TextField
                    onChange={(value) => updateSettings({ ytitle: value })}
                    value={settings.ytitle}
                  >
                    <Label>Y Axis Title</Label>
                    <Input />
                  </TextField>
                </div>
              </div>
            </div>
          </div>
        </div>
      </div>
    </>
  );
};
