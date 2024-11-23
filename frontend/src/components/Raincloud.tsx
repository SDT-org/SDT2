import React from "react";
import Plotly from "plotly.js-dist-min";
import createPlotlyComponent from "react-plotly.js/factory";
import { NumberInput } from "./NumberInput";
import { Layout, PlotData } from "plotly.js-dist-min";
import { DistributionData } from "../plotTypes";
import { formatTitle } from "../helpers";

const Plot = createPlotlyComponent(Plotly);

enum ColorOption {
  White = "white",
  Black = "black",
  Tomato = "tomato",
  LightBLue = "lightblue",
  LightGreen = "lightgreen",
  Purple = "plum",
  Lightcoral = "lightcoral",
  Orange = "orange",
  Yellow = "gold",
  Cyan = "cyan",
  Teal = "teal",
  Magenta = "magenta",
  Brown = "saddlebrown",
  Lime = "lime",
  Coral = "coral",
  Turquoise = "turquoise",
  Indigo = "indigo",
  Violet = "violet",
  Lavender = "lavender",
  Peach = "peachpuff",
  SkyBlue = "skyblue",
  Olive = "olive",
  Tan = "tan",
  Salmon = "salmon",
  Maroon = "maroon",
  Navy = "navy",
  Khaki = "khaki",
  Periwinkle = "periwinkle",
  Mint = "mintcream",
  Azure = "azure",
  Chartreuse = "chartreuse",
  Goldrod = "goldenrod",
  SlateBlue = "slateblue",
  LightSeaGreen = "lightseagreen",
  DarkCyan = "darkcyan",
  RosyBrown = "rosybrown",
  PaleVioletRed = "palevioletred",
  DeepPink = "deeppink",
  DarkOrange = "darkorange",
  Crimson = "crimson",
  LightSalmon = "lightsalmon",
  Orchid = "orchid",
  Thistle = "thistle",
  DarkKhaki = "darkkhaki",
  LightCoral = "lightcoral",
  MediumOrchid = "mediumorchid",
  None = "rgba(0,0,0,0)", // No color (use transparent or ignore)
}

export const Raincloud = ({
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
    plotOrientation: "horizontal",
    fillColor: "lightblue",
    bandWidth: 8,
    lineColor: "tomato",
    lineWidth: 3,
    violinOpacity: 0.5,
    markerColor: "tomato",
    markerSize: 7,
    pointPos: -1.5,
    pointOpacity: 0.5,
    points: "all",
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
  const rainCloudTrace = React.useMemo(
    () =>
      ({
        type: "violin",
        name: "",
        x: data.raw_mat,
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
        bandwidth: settings.bandWidth,
        marker: {
          visible: true,
          color: settings.markerColor,
          size: settings.markerSize,
          opacity: settings.pointOpacity,
        },
        meanline: {
          visible: true,
        },
        hovertemplate:
          "Percent Identity: %{x}<br>Percent Identity: %{y}<extra></extra>",
      }) as Partial<PlotData>,
    [data, settings],
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
            <div className="row">
              <div className="col-2">
                <div className="field">
                  <label htmlFor="fill-color">Cloud Fill Color</label>
                  <select
                    id="fill-color"
                    value={settings.fillColor}
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
                    updateValue={updateSettings}
                    min={0}
                    max={20}
                    step={1}
                  />
                </div>
              </div>
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
