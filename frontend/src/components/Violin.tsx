import React from "react";
import Plotly from "plotly.js-dist-min";
import createPlotlyComponent from "react-plotly.js/factory";
import { NumberInput } from "./NumberInput";
import { Layout, PlotData } from "plotly.js-dist-min";
import { DistributionData } from "../plotTypes";

const Plot = createPlotlyComponent(Plotly);

enum ColorOption {
  White = "white",
  Black = "black",
  Red = "tomato",
  Blue = "lightblue",
  Green = "lightgreen",
  Purple = "plum",
  Pink = "lightcoral",
  None = "none",
}

export const Violin = ({
  data,
  footer,
}: {
  data: DistributionData | undefined;
  footer?: React.ReactNode;
}) => {
  if (!data) {
    return (
      <div className="app-main centered">
        <h1>No data available</h1>
      </div>
    );
  }
  const minDataValue = Math.min(...data.raw_mat);
  const maxDataValue = Math.max(...data.raw_mat);

  const [settings, setSettings] = React.useState({
    plotTitle: "Distribution of Percent Identities",
    plotOrientation: "vertical", 
    fillColor: "lightblue",
    bandWidth: 8,
    lineColor: "tomato",
    lineWidth: 3,
    violinOpacity: 0.5,
    boxfillColor: "lightblue",
    boxWidth: 0.8,
    boxlineColor: "tomato",
    boxlineWidth: 3,
    boxOpacity: 0.5,
    whiskerWidth: 0.2,
    markerColor: "tomato",
    markerSize: 7,
    pointOrientation: "Violin",
    points: "all",
    pointPos: 0,
    showViolin: true,
    showBox: true,
    showPoints: true,
    showZeroLine: false,
    showGrid: true,
    showAxisLines: true,
    showTickLabels: true,
    showAxisLabels: true,
    bandwidth: 8,
    jitter: 0.5,
  });

  console.log(data);

  const updateSettings = (newState: Partial<typeof settings>) => {
    setSettings((previous) => {
      return {
        ...previous,
        ...newState,
      };
    });
  };
  const violinTrace = React.useMemo(
    () =>
      ({
        type: "violin",
        name: "",
        [settings.plotOrientation === "vertical" ? "y" : "x"]: data.raw_mat,
        visible: settings.showViolin,
        line: {
          color: settings.lineColor,
          width: settings.lineWidth,
        },
        fillcolor: settings.fillColor,
        opacity: settings.violinOpacity,
        points:
        settings.pointOrientation === "Violin"
          ? settings.points === "None"
            ? false
            : settings.points
          : false,
        pointpos: settings.pointPos,
        jitter: settings.jitter,
        bandwidth: settings.bandWidth,
        marker: {
          visible: settings.showPoints,
          color: settings.markerColor,
          size: settings.markerSize,
        },
        meanline: {
          visible: true,
        },
        hovertemplate:
          "Percent Identity: %{x}<br>Percent Identity: %{y}<extra></extra>",
      } as Partial<PlotData>),
    [data, settings]
  );
  const boxTrace = React.useMemo(
    () =>
      ({
        type: "box",
        name: "",
        [settings.plotOrientation === "vertical" ? "y" : "x"]: data.raw_mat,
        visible: settings.showBox,
        boxpoints:
        settings.pointOrientation === "Box"
          ? settings.points === "None"
            ? false
            : settings.points
          : false,
        pointpos: settings.pointPos,
        jitter: settings.jitter,
        line: {
          color: settings.boxlineColor,
          width: settings.boxlineWidth,
        },
        opacity: settings.boxOpacity,
        whiskerwidth:settings.whiskerWidth,
        marker: {
          visible: settings.showBox,
        },
        fillcolor:settings.boxfillColor,
        hovertemplate:
          "Percent Identity: %{x}<br>Percent Identity: %{y}<extra></extra>",
      } as Partial<PlotData>),
    [data, settings]
  );

    const layout = React.useMemo(() => {
      const isVertical = settings.plotOrientation === "vertical";
    
      return {
        title: settings.plotTitle,
        xaxis: isVertical
          ? {
              fixedrange: true,
              dtick: 1,
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
              dtick: 1,
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
    }, [data, settings]);
  

  return (
    <>
      <div className="app-sidebar">
        <div className="app-sidebar-toolbar">
          <div className="form">
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
              <div className="field">
                <label>Plot Orientation</label>
                <div style={{ display: "flex", gap: "100px" }}></div>
                <div>
                  <label>
                    <input
                      type="radio"
                      name="orientation"
                      value="vertical"
                      checked={settings.plotOrientation === "vertical"}
                      onChange={(e) => updateSettings({ plotOrientation: e.target.value })}
                    />
                    Vertical
                  </label>
                  <label>
                    <input
                      type="radio"
                      name="orientation"
                      value="horizontal"
                      checked={settings.plotOrientation === "horizontal"}
                      onChange={(e) => updateSettings({ plotOrientation: e.target.value })}
                    />
                    Horizontal
                  </label>
                </div>
              </div>
            </div>
            <div className="group">
              <div className="field">
                <label className="header">
                  <input
                    type="checkbox"
                    checked={settings.showViolin}
                    onChange={(e) =>
                      updateSettings({ showViolin: e.target.checked })
                    }
                  />
                  Violin Plot
                </label>
              </div>
              <div className="row">
                <div className="col-2">
                  <div className="field">
                    <label htmlFor="fill-color">Fill Color</label>
                    <select
                      id="fill-color"
                      value={settings.fillColor}
                      disabled={!settings.showViolin}
                      onChange={(e) =>
                        updateSettings({ fillColor: e.target.value })
                      }
                    >
                      {Object.entries(ColorOption).map(([key, value]) => (
                        <option key={key} value={value}>
                          {key}
                        </option>
                      ))}
                    </select>
                  </div>
                  <div className="field">
                    <NumberInput
                      label="Band Width"
                      field="bandWidth"
                      value={settings.bandWidth}
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
                        updateSettings({ lineColor: e.target.value })
                      }
                    >
                      {Object.entries(ColorOption).map(([key, value]) => (
                        <option key={key} value={value}>
                          {key}
                        </option>
                      ))}
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
              <div className="field">
                <label className="header">
                  <input
                    type="checkbox"
                    checked={settings.showBox}
                    onChange={(e) =>
                      updateSettings({ showBox: e.target.checked })
                    }
                  />
                  Box Plot
                </label>
              </div>
              <div className="row">
                <div className="col-2">
                  <div className="field">
                    <label htmlFor="box-fill-color">Fill Color</label>
                    <select
                      id="box-fill-color"
                      value={settings.boxfillColor}
                      disabled={!settings.showBox}
                      onChange={(e) =>
                        updateSettings({ boxfillColor: e.target.value })
                      }
                    >
                      {Object.entries(ColorOption).map(([key, value]) => (
                        <option key={key} value={value}>
                          {key}
                        </option>
                      ))}
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
              </div>
              <div className="row">
                <div className="col-2">
                  <div className="field">
                    <label htmlFor="box-line-color">Line Color</label>
                    <select
                      id="box-line-color"
                      value={settings.boxlineColor}
                      disabled={!settings.showBox}
                      onChange={(e) =>
                        updateSettings({ boxlineColor: e.target.value })
                      }
                    >
                      {Object.entries(ColorOption).map(([key, value]) => (
                        <option key={key} value={value}>
                          {key}
                        </option>
                      ))}
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
              </div>
              <div className="row">
                <div className="col-2">
                <div className="field">
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
                  </div>
                  <div className="field">
                    <NumberInput
                      label="Whisker Width"
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
                <label>Show Points</label>
                <div style={{ display: "flex", gap: "20px" }}>
                  <label>
                    <input
                      type="radio"
                      name="pointOrientation"
                      value="Violin"
                      checked={settings.pointOrientation === "Violin"}
                      onChange={(e) => updateSettings({ pointOrientation: e.target.value })}
                    />
                    Violin
                  </label>
                  <label>
                    <input
                      type="radio"
                      name="pointOrientation"
                      value="Box"
                      checked={settings.pointOrientation === "Box"}
                      onChange={(e) => updateSettings({ pointOrientation: e.target.value })}
                    />
                    Box
                  </label>
                  <label>
                    <input
                      type="radio"
                      name="pointOrientation"
                      value="None"
                      checked={settings.pointOrientation === "None"}
                      onChange={(e) => updateSettings({ pointOrientation: e.target.value })}
                    />
                    None
                  </label>
                </div>
              </div>

              <div className="row">
                <div className="col-2">
                  <div className="field">
                    <label htmlFor="pointpos">Point Position</label>
                    <select
                      id="pointpos"
                      value={settings.pointPos}
                      disabled={!settings.showPoints}
                      onChange={(e) =>
                        updateSettings({
                          pointPos: parseFloat(e.target.value),
                        })
                      }
                    >
                      {[
                        { label: "Far Left", value: -2 },
                        { label: "Left", value: -1 },
                        { label: "Mid", value: 0 },
                        { label: "Right", value: 1 },
                        { label: "Far Right", value: 2 },
                      ].map((option) => (
                        <option key={option.value} value={option.value}>
                          {option.label}
                        </option>
                      ))}
                    </select>
                  </div>
                  <div className="field">
                    <label htmlFor="points">Points</label>
                    <select
                      id="points"
                      value={settings.points}
                      disabled={!settings.showPoints}
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
                            {value === "None" ? "None" : value}
                          </option>
                        )
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
                        updateSettings({ markerColor: e.target.value })
                      }
                    >
                      {Object.entries(ColorOption).map(([key, value]) => (
                        <option key={key} value={value}>
                          {key}
                        </option>
                      ))}
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
              </div>
            </div>
          </div>
          <div className="app-sidebar-footer">{footer}</div>
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
