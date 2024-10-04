import React from "react";
import Plotly from "plotly.js-dist-min";
import createPlotlyComponent from "react-plotly.js/factory";
import { Colorscale, HeatmapData, HeatmapSettings } from "../plotTypes";
import { NumberInput } from "./NumberInput";
import { colorScales as defaultColorScales } from "../colorScales";
import tinycolor from "tinycolor2";

const Plot = createPlotlyComponent(Plotly);

export const Heatmap = ({
  settings,
  updateSettings,
  data,
  tickText,
  footer,
}: {
  settings: HeatmapSettings;
  updateSettings: (_: Partial<HeatmapSettings>) => void;
  data: HeatmapData;
  tickText: string[];
  footer?: React.ReactNode;
}) => {
  const colorScales = {
    ...defaultColorScales,
    Discrete: [
      [0, "rgb(0, 0, 255)"],
      [settings.cutoff_2 / 100 - 0.01, "rgb(0, 0, 255)"],
      [settings.cutoff_2 / 100, "rgb(0, 255, 0)"],
      [settings.cutoff_1 / 100, "rgb(0, 255, 0)"],
      [settings.cutoff_1 / 100 + 0.01, "rgb(200, 0, 0)"],
      [1, "rgb(255, 0, 0)"],
    ],
  };

  console.log(data);
  console.log(colorScales["Discrete"]);

  const annotations = React.useMemo(() => {
    const x: number[] = [];
    const y: number[] = [];
    const text: string[] = [];
    const textColors: (string | null)[] = [];

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
          text.push("");
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
  }, [tickText, settings]);

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
              {settings.colorscale === "Discrete" ? (
                <>
                  <div className="col-2">
                    <div className="field">
                      <div></div>
                      <NumberInput
                        label="Cutoff 1"
                        field="cutoff_1"
                        value={settings.cutoff_1}
                        updateValue={updateSettings}
                        min={settings.cutoff_2 + 1}
                        max={100}
                        step={1}
                      />
                    </div>
                    <div className="field">
                      <NumberInput
                        label="Cutoff 2"
                        field="cutoff_2"
                        value={settings.cutoff_2}
                        updateValue={updateSettings}
                        min={0}
                        max={settings.cutoff_1 - 1}
                        step={1}
                      />
                    </div>
                  </div>
                </>
              ) : null}
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
                  <label htmlFor="cbar-shrink">Height</label>
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
                  <label htmlFor="cbar-aspect">Width</label>
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
                  max={settings.vmax - 1}
                  step={1}
                />
                <NumberInput
                  label="Max"
                  field="vmax"
                  value={settings.vmax}
                  updateValue={updateSettings}
                  min={settings.vmin + 1}
                  max={100}
                  step={1}
                />
              </div>
            </div>
          </div>
        </div>
        <div className="app-sidebar-footer">{footer}</div>
      </div>
      <div className="app-main">
        {data ? (
          <Plot
            id="plot"
            data={[
              {
                z: data,
                colorscale: colorScales[settings.colorscale],
                reversescale: settings.reverse,
                type: "heatmap",
                hovertemplate:
                  "Seq 1: %{y}<br>Seq 2: %{x}<br>Percent Pairwise Identity: %{z}<extra></extra>",
                // @ts-ignore
                hoverongaps: false,
                transpose: false,
                zsmooth: false,
                autosize: true,
                set_aspect: 3,
                showscale: settings.showscale,
                zmin: settings.colorscale === "Discrete" ? 0 : settings.vmin,
                zmax: settings.colorscale === "Discrete" ? 100 : settings.vmax,
                // Gap between columns of cells. Actual 0 values cause blurriness on macOS
                xgap: settings.cellspace || 0.001,
                ygap: settings.cellspace || 0.001,
                colorbar: {
                  len: settings.cbar_shrink,
                  thickness: parseInt(settings.cbar_aspect, 10),
                  xpad: parseInt(settings.cbar_pad, 10),
                },
              },
              {
                type: "scattergl",
                mode: "text",
                // @ts-ignore
                scrollZoom: true,
                textfont: {
                  family: "Arial",
                  size: settings.annotation_font_size,
                  color: annotations?.textColors ?? "white",
                },
                hoverinfo: "skip",
                ...(settings.annotation
                  ? {
                      x: annotations.x,
                      y: annotations.y,
                      text: annotations.text,
                    }
                  : {}),
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
                l: 25,
                r: 25,
                b: 25,
                t: 25,
                pad: 5,
              },
              autosize: true,
              scaleanchor: true,
              dragmode: "pan",
              hovermode: "closest",
              xaxis: {
                minallowed: -2,
                maxallowed: tickText.length + 2,
                automargin: true,
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
                minallowed: -2,
                maxallowed: tickText.length + 2,
                automargin: true,
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
              },
            }}
            style={{ width: "100%", height: "100%" }}
          />
        ) : null}
      </div>
    </>
  );
};
