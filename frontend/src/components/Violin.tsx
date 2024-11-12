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

  const [settings, setSettings] = React.useState({
    lineColor: "tomato",
    lineWidth: 3,
    lineShape: "linear",
    boxlineColor: "tomato",
    boxlineWidth: 3,
    barOutlineWidth: 1,
    fillColor: "lightblue",
    markerColor: "tomato",
    markerSize: 7,
    points: "all",
    showViolin: true,
    showBox: true,
    showPoints: true,
    showZeroLine: false,
    showGrid: true,
    plotTitle: "Distribution of Percent Identities",
    showTickLabels: true,
    showAxisLabels: true,
    bandwidth: 8,
    jitter: .3
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
        type: 'violin',
        y: data.raw_mat, // x-y switches h-v
        line: {
          color: settings.lineColor,
          width: settings.lineWidth
        },
        points: settings.points,
        visible: settings.showViolin,
        pointpos: 0,
        fillcolor: settings.fillColor,
        opacity: 0.6,
        box: {
          visible: settings.showBox,
          line:{ 
            color: settings.boxlineColor,
            width: settings.boxlineWidth
          }
        },

        marker: {
        color:settings.markerColor,
        size:settings.markerSize
        },
        meanline: {
          visible: true 
        },
        bandwidth: 8,
        jitter: .5,
        hovertemplate: "Percent Identity: %{x}<br>Percent Identity: %{y}<extra></extra>",
      }) as Partial<PlotData>,
    [data, settings],
  );

  const layout = React.useMemo(
    () =>
      ({
        title: settings.plotTitle,
        xaxis: {
          fixedrange: true,
          dtick: 1,
          showticklabels: settings.showTickLabels,
        },
        yaxis: {
          side: "left",
          rangemode: "tozero",
          fixedrange: true,
          zeroline: false,
          showgrid: settings.showGrid,
          showticklabels: settings.showTickLabels,
        },

        dragmode: "pan",
        fixedrange: true,
        barmode: "overlay",
        showlegend: false,
        margin: { l: 50, r: 50, t: 50, b: 50 },
      }) as Partial<Layout>,
    [data, settings],
  );
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
                          updateSettings({
                            showGrid: !settings.showGrid,
                          })
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
              </div>
              <div className="col-2">
                <div className="field">
                  <label htmlFor="showBox">
                    <input
                      type="checkbox"
                      name="showBox"
                      id="showBox"
                      checked={settings.showBox}
                      onChange={() =>
                        updateSettings({
                          showBox: !settings.showBox,
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
            <div className="group">
              <div className="field">
                <label className="header">
                  <input
                    type="checkbox"
                    checked={settings.showViolin}
                    onChange={(e) =>
                      updateSettings({ showViolin: e.target.checked, })
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
                      onChange={(e) => updateSettings({ fillColor: e.target.value })}
                    >
                      {Object.entries(ColorOption).map(([key, value]) => (
                        <option key={key} value={value}>
                          {key}
                        </option>
                      ))}
                    </select>
                  </div>
                </div>
                    
                <div className="col-2">
                  <div className="field">
                    <label htmlFor="points">Points</label>
                    <select
                      id="points"
                      value={settings.points}
                      onChange={(e) =>
                        updateSettings({
                          points: e.target.value as "all" | "outliers" | "suspectedoutliers" | "false",
                        })
                      }
                    >
                      {["all", "outliers", "suspectedoutliers", "false"].map((value) => (
                        <option key={value} value={value}>
                          {value}
                        </option>
                      ))}
                    </select>
                  </div>
                </div>
              </div>
                    
              <div className="row">
                <div className="col-2">
                  <div className="field">
                    <label htmlFor="line-color">Line Color</label>
                    <select
                      id="line-color"
                      value={settings.lineColor}
                      disabled={!settings.showViolin}
                      onChange={(e) => updateSettings({ lineColor: e.target.value })}
                    >
                      {Object.entries(ColorOption).map(([key, value]) => (
                        <option key={key} value={value}>
                          {key}
                        </option>
                      ))}
                    </select>
                  </div>
                </div>
                    
                <div className="col-2">
                  <NumberInput
                    label="Line Width"
                    field="lineWidth"
                    value={settings.lineWidth}
                    disabled={!settings.showViolin}
                    updateValue={updateSettings}
                    min={0}
                    max={5}
                    step={1}
                  />
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
              <div className="col-2">
                <div className="field">
                  <label htmlFor="box-line-color">Color</label>
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
                <NumberInput
                  label="Width"
                  field="boxlineWidth"
                  value={settings.boxlineWidth}
                  isDisabled={!settings.showBox}
                  updateValue={updateSettings}
                  min={0}
                  max={10}
                  step={1}
                />
              </div>
            </div>
            <div className="group">
              <div className="field">
                <label className="header">
                  <input
                    type="checkbox"
                    checked={settings.showPoints}
                    onChange={(e) =>
                      updateSettings({ showPoints: e.target.checked })
                    }
                  />
                  Show Points
                </label>
              </div>
              <div className="col-2">
                <div className="field">
                  <label htmlFor="markerColor">Color</label>
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
                    label="Size"
                    field="markerSize"
                    value={settings.markerSize}
                    isDisabled={!settings.points}
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
      <div className="app-main">
        <Plot
          data={[
            settings.showViolin ? violinTrace : {},
            // settings.showBox ? linePlotTrace : {},
            // settings. ? scatterPlotTrace : {}, // Add a new trace for the scatter plot
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
