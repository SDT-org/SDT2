import React from "react";
import Plotly from "plotly.js-dist-min";
import createPlotlyComponent from "react-plotly.js/factory";
import { Colorscale, HeatmapData, HeatmapSettings } from "../plotTypes";
import { NumberInput } from "./NumberInput";
import { colorScales } from "../colorScales";
import tinycolor from "tinycolor2";
import { SaveImage } from "./SaveImage";

const Plot = createPlotlyComponent(Plotly);

export const Heatmap = ({
  settings,
  updateSettings,
  data,
  tickText,
}: {
  settings: HeatmapSettings;
  updateSettings: (_: Partial<HeatmapSettings>) => void;
  data: HeatmapData | undefined;
  tickText: string[];
}) => {
  // const [settings, setHeatmapSettings] = React.useState<HeatmapSettings>(
  //   {
  //     colorscale: "Portland",
  //     reverse: false,
  //     vmax: 100,
  //     vmin: 65,
  //     cellspace: 1,
  //     annotation: false,
  //     annotation_font_size: 10,
  //     annotation_rounding: 0,
  //     annotation_alpha: "0",
  //     color: "white",
  //     showscale: true,
  //     cbar_shrink: 1,
  //     cbar_aspect: "25",
  //     cbar_pad: "10",
  //     axis_labels: true,
  //     axlabel_xrotation: 270,
  //     axlabel_xfontsize: 12,
  //     axlabel_yrotation: 0,
  //     axlabel_yfontsize: 12,
  //   },
  // );

  // const updateSettings = (newState: Partial<HeatmapSettings>) => {
  //   setHeatmapSettings(previous => {
  //     return {
  //       ...previous,
  //       ...newState,
  //     };
  //   });
  // };

  const plotAnnotations = React.useMemo(() => {
    const x: number[] = [];
    const y: number[] = [];
    const text: string[] = [];
    const textColors: (string | null)[] = [];

    if (!data) {
      return;
    }

    const normalizedData = data
      .flat()
      .map(parseFloat)
      .filter((d) => d > 0);

    // https://stackoverflow.com/a/42623277
    const [dataMin, dataMax] = normalizedData.reduce(
      ([min, max], val) => [Math.min(min, val), Math.max(max, val)],
      [Number.POSITIVE_INFINITY, Number.NEGATIVE_INFINITY],
    );

    const dataDiff = dataMax - dataMin;

    data.forEach((row, rowIndex) => {
      row.forEach((datum, columnIndex) => {
        x.push(columnIndex);
        y.push(rowIndex);

        if (datum === null) {
          text.push(datum);
        } else {
          text.push(parseFloat(datum).toFixed(settings.annotation_rounding));
        }

        let parsedDatum = parseFloat(datum);
        if (Number.isInteger(parsedDatum)) {
          parsedDatum = parsedDatum / 100;
        }

        const normalizedDatum = Math.max(0, (parsedDatum - dataMin) / dataDiff);

        let scale: (typeof colorScales)[Colorscale] | [number, string][] =
          colorScales[settings.colorscale];

        if (settings.reverse) {
          scale = [...scale]
            .reverse()
            .map((data, i) => [(scale[i] || scale[0])[0], data[1]]);
        }

        const match = scale.reduce((previous, current) =>
          Math.abs(current[0] - normalizedDatum) <=
          Math.abs(previous[0] - normalizedDatum)
            ? current
            : previous,
        );

        if (datum === null) {
          textColors.push(null);
        } else if (tinycolor(match[1]).isLight()) {
          textColors.push("black");
        } else {
          textColors.push("white");
        }
      });
    });

    return {
      x,
      y,
      text,
      textColors,
    };
  }, [tickText, data, settings]);

  console.log(settings);

  return (
    <>
      <div className="app-sidebar">
        <div className="app-sidebar-toolbar">
          <div className="form">
            <div className="group">
              <h4>Heatmap</h4>
              <div className="col-2">
                <div className="field">
                  <label htmlFor="colorscale">Colorscale</label>
                  <select
                    id="colorscale"
                    value={settings.colorscale}
                    onChange={(e) =>
                      updateSettings({
                        colorscale: e.target.value as Colorscale,
                      })
                    }
                  >
                    {Object.keys(colorScales).map((value) => (
                      <option key={value} value={value}>
                        {value}
                      </option>
                    ))}
                  </select>
                  <div className="subfield">
                    <label htmlFor="reverse">
                      <input
                        type="checkbox"
                        name="reverse"
                        id="reverse"
                        defaultChecked={settings.reverse}
                        onChange={() =>
                          updateSettings({
                            reverse: !settings.reverse,
                          })
                        }
                      />
                      Reverse
                    </label>
                  </div>
                </div>
                <div className="field">
                  <NumberInput
                    label="Cell Spacing"
                    field="cellspace"
                    value={settings.cellspace}
                    updateValue={updateSettings}
                    min={0}
                    max={20}
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
                    name="annotation"
                    id="annotation"
                    checked={settings.annotation}
                    onChange={() =>
                      updateSettings({
                        annotation: !settings.annotation,
                      })
                    }
                  />
                  Percent Identities
                </label>
              </div>
              <div className="col-2">
                <div className="field">
                  <label htmlFor="round-vals">Precision</label>
                  <select
                    id="round-vals"
                    value={settings.annotation_rounding}
                    onChange={(event) =>
                      updateSettings({
                        annotation_rounding: parseInt(
                          event.target.value,
                        ) as HeatmapSettings["annotation_rounding"],
                      })
                    }
                    disabled={!settings.annotation}
                  >
                    <option value={0}>0</option>
                    <option value={1}>1</option>
                    <option value={2}>2</option>
                  </select>
                </div>
                <NumberInput
                  label="Font Size"
                  field="annotation_font_size"
                  value={settings.annotation_font_size}
                  updateValue={updateSettings}
                  min={1}
                  max={20}
                  step={1}
                  isDisabled={!settings.annotation}
                />
              </div>
            </div>
            <div className="group">
              <div className="field">
                <label className="header">
                  <input
                    type="checkbox"
                    name="axis_labels"
                    id="axis_labels"
                    checked={settings.axis_labels}
                    onChange={() =>
                      updateSettings({
                        axis_labels: !settings.axis_labels,
                      })
                    }
                  />
                  Axis Labels
                </label>
              </div>
              <div className="col-2">
                <NumberInput
                  label="X Font Size"
                  field="axlabel_xfontsize"
                  value={settings.axlabel_xfontsize}
                  updateValue={updateSettings}
                  min={1}
                  max={40}
                  step={1}
                />
                <NumberInput
                  label="X Rotation"
                  field="axlabel_xrotation"
                  value={settings.axlabel_xrotation}
                  updateValue={updateSettings}
                  min={0}
                  max={360}
                  step={10}
                />
                <NumberInput
                  label="Y Font Size"
                  field="axlabel_yfontsize"
                  value={settings.axlabel_yfontsize}
                  updateValue={updateSettings}
                  min={1}
                  max={40}
                  step={1}
                />
                <NumberInput
                  label="Y Rotation"
                  field="axlabel_yrotation"
                  value={settings.axlabel_yrotation}
                  updateValue={updateSettings}
                  min={0}
                  max={360}
                  step={10}
                />
              </div>
            </div>
            <div className="group">
              <div className="field">
                <label className="header">
                  <input
                    type="checkbox"
                    name="showscale"
                    id="showscale"
                    checked={settings.showscale}
                    onChange={() =>
                      updateSettings({
                        showscale: !settings.showscale,
                      })
                    }
                  />
                  Scale Bar
                </label>
              </div>
              <div className="range-group">
                <div className="field">
                  <label htmlFor="cbar-shrink">Size</label>
                  <input
                    type="range"
                    name="cbar-shrink"
                    id="cbar-shrink"
                    min="0.1"
                    max="1"
                    step=".1"
                    value={settings.cbar_shrink}
                    onChange={(e) =>
                      updateSettings({
                        cbar_shrink: parseFloat(e.target.value),
                      })
                    }
                  />
                </div>
                <div className="field">
                  <label htmlFor="cbar-aspect">Aspect</label>
                  <input
                    type="range"
                    name="cbar-aspect"
                    id="cbar-aspect"
                    min="1"
                    max="100"
                    step="1"
                    value={settings.cbar_aspect}
                    onChange={(e) =>
                      updateSettings({ cbar_aspect: e.target.value })
                    }
                  />
                </div>
                <div className="field">
                  <label htmlFor="cbar-pad">Pad</label>
                  <input
                    type="range"
                    name="cbar-pad"
                    id="cbar-pad"
                    min="0"
                    max="20"
                    step="1"
                    value={settings.cbar_pad}
                    onChange={(e) =>
                      updateSettings({ cbar_pad: e.target.value })
                    }
                  />
                </div>
              </div>
              <div className="col-2">
                <NumberInput
                  label="Min"
                  field="vmin"
                  value={settings.vmin}
                  updateValue={updateSettings}
                  min={1}
                  max={100} // need to add vmax min as upper lim
                  step={1}
                />
                <NumberInput
                  label="Max"
                  field="vmax"
                  value={settings.vmax}
                  updateValue={updateSettings}
                  min={1}
                  max={100}
                  step={1}
                />
              </div>
            </div>
          </div>
        </div>
        <div className="app-sidebar-footer">
          <div className="form">
            <SaveImage plotType="Heatmap" />
          </div>
        </div>
      </div>
      <div className="app-main">
        {data ? (
          <Plot
            id="plot"
            data={[
              {
                z: data,
                colorscale: settings.colorscale,
                reversescale: settings.reverse,
                type: "heatmap",
                hovertemplate:
                  "Seq 1: %{y}<br>Seq 2: %{x}<br>Percent Pairwise Identity: %{z}<extra></extra>", //remove trace:0 bubble and only display x,y, and
                // @ts-ignore
                hoverongaps: false, // disable  hover values for NaNs
                transpose: false,
                zsmooth: false,
                autosize: true,
                set_aspect: 3,
                showscale: settings.showscale,
                zmin: settings.vmin,
                zmax: settings.vmax,
                xgap: settings.cellspace, // Gap between columns of cells
                ygap: settings.cellspace,
                colorbar: {
                  len: settings.cbar_shrink,
                  thickness: parseInt(settings.cbar_aspect, 10),
                  xpad: parseInt(settings.cbar_pad, 10),
                },
                //add colorbar switch if w want with showscale: false
              },
              {
                type: "scattergl",
                mode: "text",
                // @ts-ignore
                scrollZoom: true,
                textfont: {
                  family: "Arial",
                  size: settings.annotation_font_size,
                  color: plotAnnotations?.textColors ?? "white",
                },
                hoverinfo: "skip",
                ...(settings.annotation ? plotAnnotations : {}),
              },
            ]}
            config={{
              responsive: true,
              displayModeBar: false,
              scrollZoom: true,
              modeBarButtonsToRemove: ["sendDataToCloud", "toImage"],
              displaylogo: false,
            }}
            layout={{
              margin: {
                l: 150,
                r: 100,
                b: 150,
                t: 100,
                pad: 5,
              },
              dragmode: "pan",
              hovermode: "closest",
              xaxis: {
                tickfont: {
                  family: "Arial, sans-serif",
                  size: settings.axlabel_xfontsize,
                  color: "black",
                },
                tickangle: settings.axlabel_xrotation,
                tickvals: [...Array(tickText.length).keys()],
                ticktext: tickText,
                showticklabels: settings.axis_labels,
                // @ts-ignore
                ticks: settings.axis_labels ? false : "",
                autorange: true,
                showline: false,
                zeroline: false,
                showgrid: false,
              },
              yaxis: {
                tickfont: {
                  family: "Arial, sans-serif",
                  size: settings.axlabel_yfontsize,
                  color: "black",
                },
                tickangle: settings.axlabel_yrotation,
                tickvals: [...Array(tickText.length).keys()],
                ticktext: tickText,
                showticklabels: settings.axis_labels,
                // @ts-ignore
                ticks: settings.axis_labels ? false : "",
                autorange: "reversed",
                showline: false,
                zeroline: false,
                showgrid: false,
                // add switch for forced square cells scaleanchor: 'x',
              },
            }}
            style={{ width: "100%", height: "100%" }}
          />
        ) : null}
      </div>
    </>
  );
};
