import Plotly from "plotly.js-dist-min";
import React from "react";
import { ToggleButton } from "react-aria-components";
import createPlotlyComponent from "react-plotly.js/factory";
import tinycolor from "tinycolor2";
import { colorScales as defaultColorScales } from "../colorScales";
import type { Colorscale, HeatmapData, HeatmapSettings } from "../plotTypes";
import { NumberInput } from "./NumberInput";
import { Select, SelectItem } from "./Select";
import { Slider } from "./Slider";
import { Switch } from "./Switch";
import { Tooltip } from "./Tooltip";

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
  const discreteColorScale: Array<[number, string]> = React.useMemo(() => {
    const scales = [
      [0, "#CDF0FF"],
      [Math.max(0, settings.cutoff_2 / 100 - 0.01), "#20B9FF"],
      [settings.cutoff_2 / 100, "#C3E8D3"],
      [settings.cutoff_1 / 100, "#009942"],
      [Math.min(1, settings.cutoff_1 / 100 + 0.01), "#FFDCDD"],
    ];
    if (settings.cutoff_1 < 100) {
      scales.push([1, "#FF6167"]);
    }
    return scales as Array<[number, string]>;
  }, [settings.cutoff_1, settings.cutoff_2]);

  const colorScales = React.useMemo(
    () => ({
      ...defaultColorScales,
      Discrete: discreteColorScale,
    }),
    [discreteColorScale],
  );

  const annotations = React.useMemo(() => {
    const x: number[] = [];
    const y: number[] = [];
    const text: string[] = [];
    const textColors: (string | null)[] = [];

    const normalizedData = data
      .flat()
      .map(Number.parseFloat)
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
          text.push(
            Number.parseFloat(datum).toFixed(settings.annotation_rounding),
          );
        }

        let parsedDatum = Number.parseFloat(datum);
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
  }, [
    settings.reverse,
    settings.colorscale,
    settings.annotation_rounding,
    colorScales,
    data,
  ]);

  return (
    <>
      <div className="app-sidebar heatmap-sidebar">
        <div className="app-sidebar-toolbar">
          <div className="form">
            <div className="group">
              <div className="field padded flush">
                <div className="col-2 aligned colorscale-setting">
                  <label className="header" htmlFor="colorscale">
                    Colorscale
                  </label>
                  <div className="controls">
                    <Tooltip tooltip="Reverse colorscale" delay={600}>
                      <ToggleButton
                        aria-label="Toggle reverse colorscale"
                        id="reverse"
                        isSelected={settings.reverse}
                        onChange={() =>
                          updateSettings({
                            reverse: !settings.reverse,
                          })
                        }
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
                            <path d="m15 1 5 5-5 5M9 23l-5-5 5-5" />
                            <path d="M4 18h13a6 6 0 0 0 6-6v-1M20 6H7a6 6 0 0 0-6 6v1" />
                          </g>
                        </svg>
                      </ToggleButton>
                    </Tooltip>
                    <Select
                      id="colorscale"
                      items={Object.keys(colorScales).map((name) => ({
                        id: name,
                        name,
                      }))}
                      selectedKey={settings.colorscale}
                      onSelectionChange={(value) => {
                        updateSettings({
                          colorscale: value as Colorscale,
                        });
                      }}
                    >
                      {(item) => (
                        <SelectItem key={item.name} textValue={item.name}>
                          {item.name}
                        </SelectItem>
                      )}
                    </Select>
                  </div>
                </div>
              </div>

              <div className="drawer">
                <Slider
                  label="Cell Spacing"
                  labelClassName="sublabel"
                  id="cellspace"
                  onChange={(value) => updateSettings({ cellspace: value })}
                  minValue={0}
                  maxValue={20}
                  value={settings.cellspace}
                />

                {settings.colorscale === "Discrete" ? (
                  <div className="col-2">
                    <div className="field">
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
                ) : null}
              </div>
            </div>
            <div className="group">
              <Switch
                isSelected={settings.annotation}
                onChange={(value) =>
                  updateSettings({
                    annotation: value,
                  })
                }
              >
                Percent Identities
              </Switch>
              <div
                className="drawer"
                data-hidden={!settings.annotation}
                aria-hidden={!settings.annotation}
              >
                <div className="col-2">
                  <div className="field">
                    <label htmlFor="round-vals">Precision</label>
                    <select
                      id="round-vals"
                      value={settings.annotation_rounding}
                      onChange={(event) =>
                        updateSettings({
                          annotation_rounding: Number.parseInt(
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
            </div>
            <div className="group">
              <Switch
                isSelected={settings.axis_labels}
                onChange={(value) =>
                  updateSettings({
                    axis_labels: value,
                  })
                }
              >
                Axis Labels
              </Switch>
              <div
                className="drawer"
                data-hidden={!settings.axis_labels}
                aria-hidden={!settings.axis_labels}
              >
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
            </div>
            <div className="group">
              <Switch
                isSelected={settings.showscale}
                onChange={(value) =>
                  updateSettings({
                    showscale: value,
                  })
                }
              >
                Scale Bar
              </Switch>
              <div
                className="drawer"
                data-hidden={!settings.showscale}
                aria-hidden={!settings.showscale}
              >
                <div className="range-group">
                  <Slider
                    label="Height"
                    minValue={0.1}
                    maxValue={1}
                    step={0.1}
                    value={settings.cbar_shrink}
                    onChange={(value) => updateSettings({ cbar_shrink: value })}
                  />
                  <Slider
                    label="Width"
                    id="cbar-aspect"
                    minValue={1}
                    maxValue={100}
                    step={1}
                    value={settings.cbar_aspect}
                    onChange={(value) => updateSettings({ cbar_aspect: value })}
                  />
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
                    isDisabled={settings.colorscale === "Discrete"}
                  />
                  <NumberInput
                    label="Max"
                    field="vmax"
                    value={settings.vmax}
                    updateValue={updateSettings}
                    min={settings.vmin + 1}
                    max={100}
                    step={1}
                    isDisabled={settings.colorscale === "Discrete"}
                  />
                </div>
              </div>
            </div>
          </div>
        </div>
        {footer ? <div className="app-sidebar-footer">{footer}</div> : null}
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
                  thickness: settings.cbar_aspect,
                  xpad: settings.cbar_pad,
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
              plot_bgcolor: "rgba(0,0,0,0)",
              paper_bgcolor: "rgba(0,0,0,0)",
              uirevision: "true",
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
