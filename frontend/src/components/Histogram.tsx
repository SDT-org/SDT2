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
enum LineOption {
  Bar = "hvh",
  Spline = "spline",
  Linear = "linear",
}
enum MarkerOption {
  Circle = "circle",
  Square = "square",
  Diamond = "diamond",
  Hexagram = "hexagram",
  TriangleUp = "triangle-up",
  TriangleDown = "triangle-down",
  Pentagon = "pentagon",
  CircleOpen = "circle-open",
  SquareOpen = "square-open",
  DiamondOpen = "diamond-open",
  HexagramOpen = "hexagram-open",
  TriangleUpOpen = "triangle-up-open",
  TriangleDownOpen = "triangle-down-open",
  PentagonOpen = "pentagon-open",
}

export const Histogram = ({
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
    barlineColor: "tomato",
    barOutlineWidth: 1,
    barColor: "lightblue",
    showHistogram: true,
    showLinePlot: false,
    showScatterPlot: false,
    markerSymbol: "square",
    markerColor: "tomato",
    markerSize: 7,
    showGrid: true,
    showLine: true,
    showZeroLine: true,
    plotTitle: "Distribution of Percent Identities",
    showTickLabels: true,
    showAxisLabels: true,
  });

  const updateSettings = (newState: Partial<typeof settings>) => {
    setSettings((previous) => {
      return {
        ...previous,
        ...newState,
      };
    });
  };

  const histogramTrace = React.useMemo(
    () =>
      ({
        type: "bar",
        x: data.x,
        y: data.y,
        marker: {
          color: settings.barColor,
          line: {
            width: settings.barOutlineWidth,
            color: settings.barlineColor,
          },
        },
        name: "Histogram",
        hovertemplate:
          "Percent Identity: %{x}<br>Proportion: %{y}<extra></extra>",
      }) as Partial<PlotData>,
    [data, settings],
  );

  const linePlotTrace = React.useMemo(
    () =>
      ({
        x: data.x,
        y: data.y,
        type: "scatter",
        line: {
          shape: settings.lineShape,
          color: settings.lineColor,
          width: settings.lineWidth,
        },
        name: "Line Plot",
        hovertemplate:
          "Percent Identity: %{x}<br>Proportion: %{y}<extra></extra>",
      }) as Partial<PlotData>,
    [data, settings],
  );

  const scatterPlotTrace = React.useMemo(
    () =>
      ({
        x: data.x,
        y: data.y,
        type: "scatter",
        mode: "markers",
        marker: {
          size: settings.markerSize,
          symbol: settings.markerSymbol,
          color: settings.markerColor,
        },
        name: "Line Plot",
        hovertemplate:
          "Percent Identity: %{x}<br>Proportion: %{y}<extra></extra>",
      }) as Partial<PlotData>,
    [data, settings],
  );

  const layout = React.useMemo(
    () =>
      ({
        title: settings.plotTitle,
        xaxis: {
          title: settings.showAxisLabels ? "Percent Pairwise Identity" : "",
          range: [
            Math.floor(Math.min(...data.x)),
            Math.ceil(Math.max(...data.x)),
          ],
          fixedrange: true,
          dtick: 1,
          showline: settings.showLine,
          zeroline: settings.showLine,
          showgrid: settings.showGrid,
          showticklabels: settings.showTickLabels,
        },
        yaxis: {
          title: settings.showAxisLabels
            ? "Proportion of Pairwise Identities"
            : "",
          side: "left",
          rangemode: "tozero",
          fixedrange: true,
          showline: settings.showLine,
          zeroline: settings.showLine,
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
                  <label htmlFor="showLine">
                    <input
                      type="checkbox"
                      name="showLine"
                      id="showLine"
                      checked={settings.showLine}
                      onChange={() =>
                        updateSettings({
                          showLine: !settings.showLine,
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
                    checked={settings.showHistogram}
                    onChange={(e) =>
                      updateSettings({ showHistogram: e.target.checked })
                    }
                  />
                  Bar Plot
                </label>
              </div>
              <div className="field">
                <label htmlFor="bar-color">Color</label>
                <select
                  id="bar-color"
                  value={settings.barColor}
                  disabled={!settings.showHistogram}
                  onChange={(e) => updateSettings({ barColor: e.target.value })}
                >
                  {Object.entries(ColorOption).map(([key, value]) => (
                    <option key={key} value={value}>
                      {key}
                    </option>
                  ))}
                </select>
              </div>
              <div className="col-2">
                <div className="field">
                  <label htmlFor="bar-line-color">Outline</label>
                  <select
                    id="bar-line-color"
                    value={settings.barlineColor}
                    disabled={!settings.showHistogram}
                    onChange={(e) =>
                      updateSettings({ barlineColor: e.target.value })
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
                  label="Outline Width"
                  field="barOutlineWidth"
                  value={settings.barOutlineWidth}
                  disabled={!settings.showHistogram}
                  updateValue={updateSettings}
                  min={0}
                  max={5}
                  step={1}
                />
              </div>
            </div>
            <div className="group">
              <div className="field">
                <label className="header">
                  <input
                    type="checkbox"
                    checked={settings.showLinePlot}
                    onChange={(e) =>
                      updateSettings({ showLinePlot: e.target.checked })
                    }
                  />
                  Line Plot
                </label>
              </div>

              <div className="field">
                <label htmlFor="line-shape">Shape</label>
                <select
                  id="line-shape"
                  value={settings.lineShape}
                  disabled={!settings.showLinePlot}
                  onChange={(e) =>
                    updateSettings({ lineShape: e.target.value })
                  }
                >
                  {Object.entries(LineOption).map(([key, value]) => (
                    <option key={key} value={value}>
                      {key}
                    </option>
                  ))}
                </select>
              </div>
              <div className="col-2">
                <div className="field">
                  <label htmlFor="line-color">Color</label>
                  <select
                    id="line-color"
                    value={settings.lineColor}
                    disabled={!settings.showLinePlot}
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
                <NumberInput
                  label="Width"
                  field="lineWidth"
                  value={settings.lineWidth}
                  isDisabled={!settings.showLinePlot}
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
                    checked={settings.showScatterPlot}
                    onChange={(e) =>
                      updateSettings({ showScatterPlot: e.target.checked })
                    }
                  />
                  Markers
                </label>
              </div>
              <div className="field">
                <label htmlFor="marker-symbol">Symbol</label>
                <select
                  id="marker-symbol"
                  value={settings.markerSymbol}
                  disabled={!settings.showScatterPlot}
                  onChange={(e) =>
                    updateSettings({ markerSymbol: e.target.value })
                  }
                >
                  {Object.entries(MarkerOption).map(([key, value]) => (
                    <option key={key} value={value}>
                      {key}
                    </option>
                  ))}
                </select>
              </div>
              <div className="col-2">
                <div className="field">
                  <label htmlFor="markerColor">Color</label>
                  <select
                    id="markerColor"
                    value={settings.markerColor}
                    disabled={!settings.showScatterPlot}
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
                    isDisabled={!settings.showScatterPlot}
                    updateValue={updateSettings}
                    min={0}
                    max={20}
                    step={1} //want to change to .5 but breaks
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
            settings.showHistogram ? histogramTrace : {},
            settings.showLinePlot ? linePlotTrace : {},
            settings.showScatterPlot ? scatterPlotTrace : {}, // Add a new trace for the scatter plot
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
