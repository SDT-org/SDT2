import Plotly from "plotly.js-dist-min";
import type { Layout, PlotData } from "plotly.js-dist-min";
import React from "react";
import {
  Label,
  Radio,
  RadioGroup,
  ToggleButton,
  ToggleButtonGroup,
} from "react-aria-components";
import createPlotlyComponent from "react-plotly.js/factory";
import type { ColorString, Colors } from "../colors";
import type { DataSets, DistributionState } from "../distributionState";
import type { DistributionData } from "../plotTypes";
import { ColorOptions } from "./ColorOptions";
import { ColorPicker } from "./ColorPicker";
import { NumberInput } from "./NumberInput";
import { Slider } from "./Slider";
import { Switch } from "./Switch";
import { Tooltip } from "./Tooltip";

const Plot = createPlotlyComponent(Plotly);

export const Violin = ({
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
  const minDataValue = Math.min(...dataSet);
  const maxDataValue = Math.max(...dataSet);

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
          settings.showPoints && settings.pointOrientation === "Violin"
            ? settings.points === false
              ? false
              : settings.points
            : false,
        pointpos: settings.pointPos,
        jitter: settings.jitter,
        bandwidth: settings.bandwidth,
        marker: {
          visible: settings.showPoints,
          color: settings.markerColor,
          size: settings.markerSize,
          opacity: settings.pointOpacity,
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
        opacity: settings.boxOpacity,
        whiskerwidth: settings.whiskerWidth,
        marker: {
          visible: settings.showBox,
          color: settings.markerColor,
          opacity: settings.pointOpacity,
        },
        fillcolor: settings.boxfillColor,
        hovertemplate:
          "Percent Identity: %{x}<br>Percent Identity: %{y}<extra></extra>",
      }) as Partial<PlotData>,
    [dataSet, settings],
  );

  const layout = React.useMemo(() => {
    const isVertical = settings.plotOrientation === "vertical";

    return {
      title: settings.plotTitle,
      uirevision: "true",
      xaxis: isVertical
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
            title: settings.showAxisLabels
              ? "Percent Pairwise Identity"
              : undefined,
            range: [minDataValue - 20, maxDataValue + 20],
          },
      yaxis: isVertical
        ? {
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
      boxgap: settings.boxWidth,
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
                      ].includes(key) && settings[key as keyof typeof settings],
                  )}
                  onSelectionChange={(value) =>
                    updateSettings({
                      showGrid: value.has("showGrid"),
                      showTickLabels: value.has("showTickLabels"),
                      showAxisLines: value.has("showAxisLines"),
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
              </div>
            </div>
            <div className="group">
              <Switch
                isSelected={settings.showViolin}
                onChange={(value) => updateSettings({ showViolin: value })}
              >
                Violin
              </Switch>
              <div
                className="drawer"
                data-hidden={!settings.showViolin}
                aria-hidden={!settings.showViolin}
              >
                <Label className="header">Band</Label>
                <div className="field col-2 has-color">
                  <ColorPicker
                    value={settings.fillColor}
                    onChange={(value) => {
                      updateSettings({
                        fillColor: value.toString() as ColorString,
                      });
                    }}
                  />
                  <Slider
                    label="Width"
                    defaultValue={settings.bandwidth}
                    isDisabled={!settings.showViolin}
                    onChangeEnd={(value) =>
                      updateSettings({ bandwidth: value })
                    }
                    minValue={0}
                    maxValue={20}
                    step={1}
                  />
                </div>
                <Label className="header">Line</Label>
                <div className="field col-2 has-color">
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
                    isDisabled={!settings.showViolin}
                    onChange={(value) => updateSettings({ lineWidth: value })}
                    minValue={0}
                    maxValue={20}
                    step={1}
                  />
                </div>
              </div>
            </div>
            <div className="group">
              <Switch
                isSelected={settings.showBox}
                onChange={(value) => updateSettings({ showBox: value })}
              >
                Box
              </Switch>
              <div
                className="drawer"
                data-hidden={!settings.showBox}
                aria-hidden={!settings.showBox}
              >
                <div className="col-2">
                  <div className="field">
                    <label htmlFor="box-fill-color">Fill Color</label>
                    <select
                      id="box-fill-color"
                      value={settings.boxfillColor}
                      disabled={!settings.showBox}
                      onChange={(e) =>
                        updateSettings({
                          boxfillColor: e.target.value as ColorOption,
                        })
                      }
                    >
                      <ColorOptions />
                    </select>
                  </div>
                  <div className="field">
                    <NumberInput
                      label="Box Width"
                      field="boxWidth"
                      value={settings.boxWidth}
                      isDisabled={!settings.showBox}
                      updateValue={updateSettings}
                      min={0.5}
                      max={1}
                      step={0.05}
                      type="float"
                    />
                  </div>
                </div>
                <div className="col-2">
                  <div className="field">
                    <label htmlFor="box-line-color">Line Color</label>
                    <select
                      id="box-line-color"
                      value={settings.boxlineColor}
                      disabled={!settings.showBox}
                      onChange={(e) =>
                        updateSettings({
                          boxlineColor: e.target.value as ColorOption,
                        })
                      }
                    >
                      <ColorOptions />
                    </select>
                  </div>
                  <div className="field">
                    <NumberInput
                      label="Line Width"
                      field="boxlineWidth"
                      value={settings.boxlineWidth}
                      isDisabled={!settings.showBox}
                      updateValue={updateSettings}
                      min={0}
                      max={20}
                      step={1}
                    />
                  </div>
                </div>
                <div className="row">
                  <div className="col-2">
                    <NumberInput
                      label="Box Opacity"
                      field="boxOpacity"
                      type="float"
                      value={settings.boxOpacity}
                      isDisabled={!settings.showBox}
                      updateValue={updateSettings}
                      min={0}
                      max={1}
                      step={0.1}
                    />
                    <NumberInput
                      label="Whiskers"
                      field="whiskerWidth"
                      value={settings.whiskerWidth}
                      type="float"
                      isDisabled={!settings.showBox}
                      updateValue={updateSettings}
                      min={0}
                      max={1}
                      step={0.1}
                    />
                  </div>
                </div>
              </div>
            </div>
            <div className="group">
              <Switch
                isSelected={settings.showPoints}
                onChange={(value) => updateSettings({ showPoints: value })}
              >
                Points
              </Switch>
              <div
                className="drawer"
                data-hidden={!settings.showPoints}
                aria-hidden={!settings.showPoints}
              >
                <ToggleButtonGroup
                  data-compact
                  selectionMode="single"
                  disallowEmptySelection={true}
                  selectedKeys={[settings.pointOrientation]}
                  onSelectionChange={(value) =>
                    updateSettings({
                      pointOrientation: value.values().next()
                        .value as typeof settings.pointOrientation,
                    })
                  }
                >
                  <ToggleButton id="Violin">Violin</ToggleButton>
                  <ToggleButton id="Box" isDisabled={!settings.showBox}>
                    Box
                  </ToggleButton>
                </ToggleButtonGroup>

                <Slider
                  label="Position"
                  value={settings.pointPos}
                  isDisabled={!settings.showViolin}
                  onChange={(value) => updateSettings({ pointPos: value })}
                  minValue={-2}
                  maxValue={2}
                  step={0.1}
                />

                <div className="field">
                  <label htmlFor="points">Visibility</label>
                  <select
                    id="points"
                    value={settings.points.toString()}
                    disabled={!settings.showPoints}
                    onChange={(e) =>
                      updateSettings({
                        points: e.target.value as typeof settings.points,
                      })
                    }
                  >
                    {["all", "outliers", "suspectedoutliers", false].map(
                      (value) => (
                        <option key={value.toString()} value={value.toString()}>
                          {value === false ? "None" : value}
                        </option>
                      ),
                    )}
                  </select>
                </div>
                <div className="col-2">
                  <div className="field">
                    <label htmlFor="markerColor">Point Color</label>
                    <select
                      id="markerColor"
                      value={settings.markerColor}
                      disabled={!settings.showPoints}
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
                      isDisabled={!settings.showPoints}
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
                        isDisabled={!settings.showPoints}
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
                        isDisabled={!settings.showPoints}
                        updateValue={updateSettings}
                        min={0}
                        max={1}
                        step={0.1}
                      />
                    </div>
                  </div>
                </div>
              </div>

              <div className="group">
                <Switch
                  isSelected={settings.showMeanline}
                  onChange={(value) => updateSettings({ showMeanline: value })}
                >
                  Mean
                </Switch>
              </div>
            </div>
          </div>
          {footer ? <div className="app-sidebar-footer">{footer}</div> : null}
        </div>
      </div>
      <div className="app-main">
        <Plot
          data={[
            settings.showViolin ? violinTrace : {},
            settings.showBox ? boxTrace : {},
            // settings. ? scatterPlotTrace : {},
          ]}
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
