import Plotly from "plotly.js-cartesian-dist-min";
import React from "react";
import { Input, Label, TextField, ToggleButton } from "react-aria-components";
import createPlotlyComponent from "react-plotly.js/factory";
import useAppState, { type AppState } from "../appState";
import {
  type ColorScaleArray,
  colorScales as defaultColorScales,
} from "../colorScales";
import { plotFontMonospace, plotFontSansSerif } from "../constants";
import { formatTitle } from "../helpers";
import { useAnnotations } from "../hooks/useAnnotations";
import {
  useRelayoutHideSubtitle,
  useRelayoutUpdateTitles,
} from "../hooks/useRelayoutUpdateTitles";
import type { ColorScaleKey, HeatmapData, HeatmapSettings } from "../plotTypes";
import { NumberInput } from "./NumberInput";
import { Select, SelectItem } from "./Select";
import { Slider } from "./Slider";
import { Switch } from "./Switch";
import { Tooltip } from "./Tooltip";

const Plot = createPlotlyComponent(Plotly);

export const Heatmap = ({
  data,
  tickText,
  footer,
}: {
  data: HeatmapData;
  tickText: string[];
  footer?: React.ReactNode;
}) => {
  const {
    appState: {
      client: { heatmap: settings },
      sequences_count,
    },
    setAppState,
  } = useAppState();

  const updateSettings = React.useCallback(
    (values: Partial<AppState["client"]["heatmap"]>) =>
      setAppState((prev) => ({
        ...prev,
        client: {
          ...prev.client,
          heatmap: {
            ...prev.client.heatmap,
            ...values,
          },
        },
      })),
    [setAppState],
  );

  const discreteColorScale: ColorScaleArray = React.useMemo(() => {
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
    return scales as ColorScaleArray;
  }, [settings.cutoff_1, settings.cutoff_2]);

  const colorScales = React.useMemo(
    () => ({
      Discrete: discreteColorScale,
      ...defaultColorScales,
    }),
    [discreteColorScale],
  );

  const annotations = useAnnotations(
    settings.annotation,
    data,
    settings.vmin,
    settings.vmax,
    colorScales[settings.colorScaleKey],
    settings.reverse,
    settings.annotation_rounding,
  );

  const [textScale, setTextScale] = React.useState(1);

  const updateTextScale = React.useCallback(
    (event: Plotly.PlotRelayoutEvent) => {
      if (event["xaxis.autorange"]) {
        setTextScale(1);
        return;
      }
      if (!event["xaxis.range[0]"]) {
        return;
      }

      const dataLength = data.length;
      const earlyRamp = dataLength < 10 ? 0 : dataLength < 50 ? 10 : 0;
      const initialRange = -2;
      const x0 = event["xaxis.range[0]"] || initialRange;
      const x1 = event["xaxis.range[1]"] || initialRange;

      const length = dataLength + 2 + earlyRamp;
      const scale = Math.max(
        1,
        Math.min(4, Number((length / ((x1 - x0) * 1.5)).toFixed(1))),
      );
      setTextScale(scale);
    },
    [data.length],
  );

  const maybeWarnPerformance = React.useCallback(
    (enabled: boolean, fn: () => void) => {
      if (
        enabled &&
        sequences_count > 99 &&
        !confirm(
          "Warning: Enabling this setting may significantly impact render performance.",
        )
      ) {
        return;
      }
      fn();
    },
    [sequences_count],
  );

  const updateTitles = useRelayoutUpdateTitles(updateSettings);
  useRelayoutHideSubtitle(!settings.showTitles);

  return (
    <>
      <div className="app-sidebar heatmap-sidebar">
        <div className="app-sidebar-toolbar">
          <div className="form">
            <div className="group">
              <div className="field padded flush">
                <label className="header" htmlFor="colorscale">
                  Colorscale
                </label>
                <div className="col-2 onefr-atuo colorscale-setting">
                  <Select
                    id="colorscale"
                    items={Object.keys(colorScales).map((name) => ({
                      id: name,
                      name,
                    }))}
                    selectedKey={settings.colorScaleKey}
                    onSelectionChange={(value) => {
                      updateSettings({
                        colorScaleKey: value as ColorScaleKey,
                      });
                    }}
                  >
                    {(item) => (
                      <SelectItem key={item.name} textValue={item.name}>
                        <div className="color-scale-list-item">
                          <span
                            className="preview"
                            aria-hidden="true"
                            style={{
                              background: `linear-gradient(to right, ${colorScales[item.name as ColorScaleKey].map((v) => v[1]).join(", ")})`,
                            }}
                          />
                          <span>{formatTitle(item.name)}</span>
                        </div>
                      </SelectItem>
                    )}
                  </Select>
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
                </div>
              </div>

              <div className="drawer">
                <div className="field ">
                  <Slider
                    label="Cell Spacing"
                    labelClassName="sublabel"
                    id="cellspace"
                    onChange={(value) => updateSettings({ cellspace: value })}
                    minValue={0}
                    maxValue={20}
                    value={settings.cellspace}
                  />
                </div>

                {settings.colorScaleKey === "Discrete" ? (
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
                  maybeWarnPerformance(value, () =>
                    updateSettings({
                      annotation: value,
                    }),
                  )
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
                isSelected={settings.showTitles}
                onChange={(value) => {
                  updateSettings({
                    showTitles: value,
                  });
                }}
              >
                Plot Titles
              </Switch>
              <div
                className="drawer"
                data-hidden={!settings.showTitles}
                aria-hidden={!settings.showTitles}
              >
                <div className="col-2 auto-onefr align-items-center">
                  <Label htmlFor="font">Font Type</Label>
                  <Select
                    id="font"
                    data-compact
                    selectedKey={settings.titleFont}
                    onSelectionChange={(value) =>
                      updateSettings({
                        titleFont: value as typeof settings.titleFont,
                      })
                    }
                    items={["Sans Serif", "Monospace"].map((name) => ({
                      id: name,
                      name,
                    }))}
                  >
                    {(item) => (
                      <SelectItem textValue={item.name}>{item.name}</SelectItem>
                    )}
                  </Select>
                </div>

                <div className="field">
                  <TextField
                    onChange={(value) => updateSettings({ title: value })}
                    value={settings.title}
                  >
                    <Label>Title</Label>
                    <Input />
                  </TextField>
                </div>
                <div className="field">
                  <TextField
                    onChange={(value) => updateSettings({ subtitle: value })}
                    value={settings.subtitle}
                  >
                    <Label>Subtitle</Label>
                    <Input />
                  </TextField>
                </div>
                <div className="field">
                  <TextField
                    onChange={(value) => updateSettings({ xtitle: value })}
                    value={settings.xtitle}
                  >
                    <Label>X Axis Title</Label>
                    <Input />
                  </TextField>
                </div>
                <div className="field">
                  <TextField
                    onChange={(value) => updateSettings({ ytitle: value })}
                    value={settings.ytitle}
                  >
                    <Label>Y Axis Title</Label>
                    <Input />
                  </TextField>
                </div>
              </div>
            </div>
            <div className="group">
              <Switch
                isSelected={settings.axis_labels}
                onChange={(value) =>
                  maybeWarnPerformance(value, () =>
                    updateSettings({
                      axis_labels: value,
                    }),
                  )
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
                    isDisabled={settings.colorScaleKey === "Discrete"}
                  />
                  <NumberInput
                    label="Max"
                    field="vmax"
                    value={settings.vmax}
                    updateValue={updateSettings}
                    min={settings.vmin + 1}
                    max={100}
                    step={1}
                    isDisabled={settings.colorScaleKey === "Discrete"}
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
            onRelayout={(e) => {
              updateTitles(e);
              updateTextScale(e);
            }}
            data={[
              {
                z: data,
                colorscale: colorScales[settings.colorScaleKey],
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
                zmin: settings.colorScaleKey === "Discrete" ? 0 : settings.vmin,
                zmax:
                  settings.colorScaleKey === "Discrete" ? 100 : settings.vmax,
                // Gap between columns of cells. Actual 0 values cause blurriness on macOS
                xgap: settings.cellspace || 0.001,
                ygap: settings.cellspace || 0.001,
                colorbar: {
                  len: settings.cbar_shrink,
                  thickness: settings.cbar_aspect,
                  xpad: settings.cbar_pad,
                  title: "",
                },
              },
              {
                type: "scatter",
                mode: "text",
                // @ts-ignore
                scrollZoom: true,
                textfont: {
                  ...plotFontMonospace,
                  size: settings.annotation_font_size * textScale,
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
              editable: false,
              // showLink: true,
              // plotlyServerURL: "https://chart-studio.plotly.com",
            }}
            layout={{
              ...(settings.showTitles
                ? {
                    title: {
                      text: settings.title,
                      pad: {
                        t: 100,
                        r: 0,
                        b: 50,
                        l: 0,
                      },
                      subtitle: {
                        text: settings.subtitle,
                      },
                    },
                  }
                : {}),
              font: {
                ...(settings.titleFont === "Monospace"
                  ? plotFontMonospace
                  : plotFontSansSerif),
                // @ts-ignore
                weight: "bold",
              },
              plot_bgcolor: "rgba(0,0,0,0)",
              paper_bgcolor: "rgba(0,0,0,0)",
              uirevision: "true",
              autosize: true,
              dragmode: "pan",
              hovermode: "closest",
              xaxis: {
                ...(settings.showTitles
                  ? {
                      title: {
                        text: settings.xtitle,
                        font: {
                          ...(settings.titleFont === "Monospace"
                            ? plotFontMonospace
                            : plotFontSansSerif),
                          weight: "bold",
                        },
                      },
                    }
                  : {}),
                minallowed: -2,
                maxallowed: tickText.length + 2,
                automargin: true,
                tickfont: {
                  ...plotFontMonospace,
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
                ...(settings.showTitles
                  ? {
                      title: {
                        text: settings.ytitle,
                        font: {
                          ...(settings.titleFont === "Monospace"
                            ? plotFontMonospace
                            : plotFontSansSerif),
                          weight: "bold",
                        },
                        pad: {
                          r: 15,
                        },
                      },
                    }
                  : {}),
                minallowed: -2,
                maxallowed: tickText.length + 2,
                automargin: true,
                scaleanchor: "x",
                range: [tickText.length - 1, 0],
                autorange: false,
                tickfont: {
                  ...plotFontMonospace,
                  size: settings.axlabel_yfontsize,
                  color: "black",
                },
                tickangle: settings.axlabel_yrotation,
                tickvals: [...Array(tickText.length).keys()],
                ticktext: tickText,
                showticklabels: settings.axis_labels,
                // @ts-ignore
                ticks: settings.axis_labels ? false : "",
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
