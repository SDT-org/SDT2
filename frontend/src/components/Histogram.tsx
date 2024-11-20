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
  Tomato = "tomato",
  Blue = "lightblue",
  Green = "lightgreen",
  Purple = "plum",
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
    histlineColor: "tomato",
    histOutlineWidth: 1,
    barColor: "lightblue",
    showHistogram: true,
    binSize: 1,
    showGrid: true,
    showLine: true,
    plotTitle: "Distribution of Percent Identities",
    showTickLabels: true,
    showAxisLabels: true,
    histnorm: 'percent',
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
        type: "histogram",
        x: data.raw_mat,
        histnorm: 'percent',
        marker: {
          colorscale: 'Viridis',
          color: settings.barColor,
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
          side: "bottom",
          rangemode: "normal",
          // fixedrange: true,
          // dtick: settings.showTickLabels ? 1 : undefined, 
          showline: settings.showLine,
          zeroline: settings.showLine,
          showgrid: settings.showGrid,
          showticklabels: settings.showTickLabels,
          // automargin: true,

        },
        yaxis: {
          title: settings.showAxisLabels
            ? "Proportion of Pairwise Identities"
            : "",
          side: "left",
          rangemode: "tozero",
          // fixedrange: true,
          showline: settings.showLine,
          zeroline: settings.showLine,
          showgrid: settings.showGrid,
          showticklabels: settings.showTickLabels,
          // automargin: true,

        },

        dragmode: "pan",
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
            <div className="row">
              <div className="col-2">
                <div className="field">
                  <label htmlFor="bin-color">Bin Color</label>
                  <select
                    id="bin-color"
                    value={settings.barColor}
                    onChange={(e) =>
                      updateSettings({ barColor: e.target.value })
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
                  label="Bin Size"
                  field="binSize"
                  value={settings.binSize}
                  updateValue={updateSettings}
                  type="float"
                  min={0.5}  
                  max={5}
                  step={0.5}
                />
              </div>
              <div className="col-2">
                <div className="field">
                  <label htmlFor="hist-line-color">Outline</label>
                  <select
                    id="hist-line-color"
                    value={settings.histlineColor}
                    onChange={(e) =>
                      updateSettings({ histlineColor: e.target.value })
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
                  field="histOutlineWidth"
                  value={settings.histOutlineWidth}
                  updateValue={updateSettings}
                  min={0}
                  max={15}
                  step={1}
                />
              </div>
            </div>
          </div>
        </div>
        <div className="app-sidebar-footer">{footer}</div>
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
          }}
          style={{ width: "100%", height: "100%" }}
        />
      </div>
    </>
  );
};
