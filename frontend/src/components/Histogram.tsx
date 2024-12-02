import Plotly from "plotly.js-dist-min";
import type { Layout, PlotData } from "plotly.js-dist-min";
import React from "react";
import createPlotlyComponent from "react-plotly.js/factory";
import type { Colors } from "../colors";
import type { DataSets, DistributionState } from "../distributionState";
import type { DistributionData } from "../plotTypes";
import { ColorOptions } from "./ColorOptions";
import { NumberInput } from "./NumberInput";

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
          "Percent Identity: %{x}<br>Percentage: %{y}<extra></extra>",
      }) as Partial<PlotData>,
    [settings, dataSet],
  );

  const layout = React.useMemo(
    () =>
      ({
        title: settings.plotTitle,
        uirevision: "true",
        xaxis: {
          title: settings.showAxisLabels ? "Percent Pairwise Identity" : "",
          side: "bottom",
          rangemode: "normal",
          fixedrange: true,
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
          fixedrange: true,
          showline: settings.showLine,
          zeroline: settings.showLine,
          showgrid: settings.showGrid,
          showticklabels: settings.showTickLabels,
          // automargin: true,
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
              <div className="row">
                <div className="col-2">
                  <div className="field">
                    <label htmlFor="bin-color">Bin Color</label>
                    <select
                      id="bin-color"
                      value={settings.barColor}
                      onChange={(e) =>
                        updateSettings({
                          barColor: e.target.value as ColorOption,
                        })
                      }
                    >
                      <ColorOptions />
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
                    <label htmlFor="hist-line-color">Outline Color</label>
                    <select
                      id="hist-line-color"
                      value={settings.histlineColor}
                      onChange={(e) =>
                        updateSettings({
                          histlineColor: e.target.value as ColorOption,
                        })
                      }
                    >
                      <ColorOptions />
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
            editable:true,
          }}
          style={{ width: "100%", height: "100%" }}
        />
      </div>
    </>
  );
};
