import Plotly from "plotly.js-dist-min";
import type { Layout, PlotData } from "plotly.js-dist-min";
import React from "react";
import { Label, Radio, RadioGroup } from "react-aria-components";
import createPlotlyComponent from "react-plotly.js/factory";
import type { ColorOption } from "../colors";
import type { DataSets, DistributionState } from "../distributionState";
import type { DistributionData } from "../plotTypes";
import { ColorOptions } from "./ColorOptions";
import { NumberInput } from "./NumberInput";
import { Slider } from "./Slider";
import { Switch } from "./Switch";

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
        opacity: settings.violinOpacity,
        meanline: {
          visible: settings.showMeanline === "Violin",
        },
        points:
          settings.pointOrientation === "Violin"
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
        boxmean: settings.showMeanline === "Box",
        boxpoints: settings.pointOrientation === "Box" && settings.points,
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
              <div className="field">
                <label htmlFor="plot-title" className="header">
                  Title
                </label>
                <input
                  id="plot-title"
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
                      Labels
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
              <RadioGroup
                value={settings.plotOrientation}
                onChange={(value) =>
                  updateSettings({
                    plotOrientation: value as typeof settings.plotOrientation,
                  })
                }
              >
                <Label>Plot Orientation</Label>
                <Radio value="vertical">Vertical</Radio>
                <Radio value="horizontal">Horizontal</Radio>
              </RadioGroup>
            </div>
            <div className="group">
              <Switch
                isSelected={settings.showViolin}
                onChange={(value) => updateSettings({ showViolin: value })}
              >
                Violin Plot
              </Switch>
              <div
                className="row"
                data-hidden={!settings.showViolin}
                aria-hidden={!settings.showViolin}
              >
                <div className="col-2">
                  <div className="field">
                    <label htmlFor="fill-color">Fill Color</label>
                    <select
                      id="fill-color"
                      value={settings.fillColor}
                      disabled={!settings.showViolin}
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
                      field="bandwidth"
                      value={settings.bandwidth}
                      isDisabled={!settings.showViolin}
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
                      disabled={!settings.showViolin}
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
                      isDisabled={!settings.showViolin}
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
              <Switch
                isSelected={settings.showBox}
                onChange={(value) => updateSettings({ showBox: value })}
              >
                Box Plot
              </Switch>
              <div
                className="row"
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
              <div className="field">
                <label htmlFor="show-points">Show Points</label>
                <div className="col-3">
                  <label>
                    <input
                      id="show-points"
                      type="radio"
                      name="pointOrientation"
                      value="Violin"
                      checked={settings.pointOrientation === "Violin"}
                      onChange={(e) =>
                        updateSettings({
                          pointOrientation: e.target
                            .value as typeof settings.pointOrientation,
                        })
                      }
                    />
                    Violin
                  </label>
                  <label>
                    <input
                      type="radio"
                      name="pointOrientation"
                      value="Box"
                      checked={settings.pointOrientation === "Box"}
                      onChange={(e) =>
                        updateSettings({
                          pointOrientation: e.target
                            .value as typeof settings.pointOrientation,
                        })
                      }
                    />
                    Box
                  </label>
                  <label>
                    <input
                      type="radio"
                      name="pointOrientation"
                      value="None"
                      checked={settings.pointOrientation === "None"}
                      onChange={(e) =>
                        updateSettings({
                          pointOrientation: e.target
                            .value as typeof settings.pointOrientation,
                        })
                      }
                    />
                    None
                  </label>
                </div>
              </div>
              <div className="field">
                <label htmlFor="show-mean">Show Mean</label>
                <div className="col-3">
                  <label>
                    <input
                      type="radio"
                      name="showMeanline"
                      value="Violin"
                      checked={settings.showMeanline === "Violin"}
                      onChange={(e) =>
                        updateSettings({
                          showMeanline: e.target
                            .value as typeof settings.showMeanline,
                        })
                      }
                    />
                    Violin
                  </label>
                  <label>
                    <input
                      type="radio"
                      name="showMeanline"
                      value="Box"
                      checked={settings.showMeanline === "Box"}
                      onChange={(e) =>
                        updateSettings({
                          showMeanline: e.target
                            .value as typeof settings.showMeanline,
                        })
                      }
                    />
                    Box
                  </label>
                  <label>
                    <input
                      type="radio"
                      name="showMeanline"
                      value="None"
                      checked={settings.showMeanline === "None"}
                      onChange={(e) =>
                        updateSettings({
                          showMeanline: e.target
                            .value as typeof settings.showMeanline,
                        })
                      }
                    />
                    None
                  </label>
                </div>
              </div>
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
                      max={2}
                      step={0.1}
                    />
                  </div>
                  <div className="field">
                    <label htmlFor="points">Points</label>
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
                          <option
                            key={value.toString()}
                            value={value.toString()}
                          >
                            {value === false ? "None" : value}
                          </option>
                        ),
                      )}
                    </select>
                  </div>
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
