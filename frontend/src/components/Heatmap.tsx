import Plotly from "plotly.js-cartesian-dist-min";
import React from "react";
import createPlotlyComponent from "react-plotly.js/factory";
import useAppState, {
  type DocState,
  type SetDocState,
  type UpdateDocState,
} from "../appState";
import {
  type ColorScaleArray,
  colorScales as defaultColorScales,
} from "../colorScales";
import { plotFontMonospace, plotFontSansSerif } from "../constants";
import { useAnnotations } from "../hooks/useAnnotations";
import {
  useRelayoutHideSubtitle,
  useRelayoutUpdateTitles,
} from "../hooks/useRelayoutUpdateTitles";
import type { HeatmapData } from "../plotTypes";
import { D3CanvasHeatmap } from "./D3CanvasHeatmap";
import { D3Heatmap } from "./D3Heatmap";
import { HeatmapSidebar } from "./HeatmapSidebar";

const Plot = createPlotlyComponent(Plotly);

export const Heatmap = ({
  data,
  tickText,
  docState,
  setDocState,
}: {
  docState: DocState;
  setDocState: SetDocState;
  updateDocState: UpdateDocState;
  data: HeatmapData;
  tickText: string[];
}) => {
  const { heatmap: settings, sequences_count } = docState;
  const { appState } = useAppState();
  const updateSettings = React.useCallback(
    (values: Partial<DocState["heatmap"]>) =>
      setDocState((prev) => ({
        ...prev,
        heatmap: {
          ...prev.heatmap,
          ...values,
        },
      })),
    [setDocState],
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

  const d3HeatmapData = React.useMemo(
    () =>
      data.flatMap((row, y) =>
        row.map((value, x) => ({
          x,
          y,
          value: Number(value),
        })),
      ),
    [data],
  );

  const [tempHeatmapComponent, tempSetHeatmapComponent] =
    React.useState("canvas");

  React.useEffect(() => {
    const handleKeydown = (event: KeyboardEvent) => {
      if (
        (event.ctrlKey || event.metaKey) &&
        [1, 2, 3].includes(Number(event.key))
      ) {
        event.preventDefault();
        const views = ["canvas", "svg", "plotly"];
        tempSetHeatmapComponent(() => views[Number(event.key) - 1] || "canvas");
      }
    };

    document.addEventListener("keydown", handleKeydown);

    return () => document.removeEventListener("keydown", handleKeydown);
  }, []);

  const elementRef = React.useRef<HTMLDivElement | null>(null);
  const [size, setSize] = React.useState<{ width: number; height: number }>({
    width: 0,
    height: 0,
  });

  const updateSize = React.useCallback(() => {
    if (elementRef.current) {
      const { offsetWidth, offsetHeight } = elementRef.current;
      setSize({ width: offsetWidth, height: offsetHeight });
    }
  }, []);

  React.useEffect(() => {
    updateSize();

    const handleResize = () => {
      updateSize();
    };

    window.addEventListener("resize", handleResize);
    return () => {
      window.removeEventListener("resize", handleResize);
    };
  }, [updateSize]);

  const updateTitles = useRelayoutUpdateTitles(updateSettings);
  useRelayoutHideSubtitle(!settings.showTitles);

  return (
    <>
      {data ? (
        <>
          <div
            className="app-main"
            ref={elementRef}
            style={{ background: "#fff" }}
          >
            <div
              style={{ position: "absolute", left: 10, top: 5, zIndex: "1000" }}
            >
              {tempHeatmapComponent}
            </div>
            {(
              appState.showExportModal
                ? appState.saveFormat !== "svg"
                : tempHeatmapComponent === "canvas"
            ) ? (
              <D3CanvasHeatmap
                data={d3HeatmapData}
                tickText={tickText}
                colorScale={colorScales[settings.colorScaleKey]}
                minVal={settings.vmin}
                maxVal={settings.vmax}
                width={size.width}
                height={size.height}
                cellSpace={settings.cellspace}
                roundTo={settings.annotation_rounding}
                showscale={settings.showscale}
                cbarHeight={settings.cbar_shrink}
                cbarWidth={settings.cbar_aspect}
                axlabel_xrotation={settings.axlabel_xrotation}
                axlabel_xfontsize={settings.axlabel_xfontsize}
                axlabel_yrotation={settings.axlabel_yrotation}
                axlabel_yfontsize={settings.axlabel_yfontsize}
                titleFont={settings.titleFont}
                showPercentIdentities={settings.annotation}
                showTitles={settings.showTitles}
                title={settings.title}
                subtitle={settings.subtitle}
              />
            ) : (
                appState.showExportModal
                  ? appState.saveFormat === "svg"
                  : tempHeatmapComponent === "svg"
              ) ? (
              <D3Heatmap
                data={d3HeatmapData}
                tickText={tickText}
                colorScale={colorScales[settings.colorScaleKey]}
                minVal={settings.vmin}
                maxVal={settings.vmax}
                width={size.width}
                height={size.height}
                cellSpace={settings.cellspace}
                showPercentIdentities={settings.annotation}
                roundTo={settings.annotation_rounding}
                showscale={settings.showscale}
                cbarHeight={settings.cbar_shrink}
                cbarWidth={settings.cbar_aspect}
                tempHeatmapComponent={"svg"}
                axlabel_xrotation={settings.axlabel_xrotation}
                axlabel_xfontsize={settings.axlabel_xfontsize}
                axlabel_yrotation={settings.axlabel_yrotation}
                axlabel_yfontsize={settings.axlabel_yfontsize}
                titleFont={settings.titleFont}
                showTitles={settings.showTitles}
                title={settings.title}
                subtitle={settings.subtitle}
              />
            ) : (
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
                    zmin:
                      settings.colorScaleKey === "Discrete" ? 0 : settings.vmin,
                    zmax:
                      settings.colorScaleKey === "Discrete"
                        ? 100
                        : settings.vmax,
                    // Gap between columns of cells. Actual 0 values cause blurriness on macOS
                    xgap: settings.cellspace || 0.001,
                    ygap: settings.cellspace || 0.001,
                    colorbar: {
                      len: settings.cbar_shrink,
                      thickness: settings.cbar_aspect,
                      xpad: settings.cbar_pad,
                      title: "",
                      // @ts-ignore
                      tickfont: {
                        ...plotFontMonospace,
                      },
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
                  displaylogo: false,
                  editable: false,
                }}
                layout={{
                  ...(settings.showTitles
                    ? {
                        title: {
                          text: settings.title,
                          pad: {
                            t: 0,
                            r: 0,
                            b: 0,
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
                  aspectratio: { x: 1, y: 1 },
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
                          scaleratio: 1,
                        }
                      : {}),
                    minallowed: -2,
                    maxallowed: tickText.length + 2,

                    tickfont: {
                      ...plotFontMonospace,
                      size: settings.axlabel_xfontsize,
                      color: "black",
                    },
                    tickangle: settings.axlabel_xrotation,
                    tickvals: [...Array(tickText.length).keys()],
                    ticktext: tickText,
                    showticklabels: settings.axis_labels,
                    ticklabelpadding: 0,
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
                    tickson: "all",
                    showticklabels: settings.axis_labels,
                    ticklabelpadding: 0,
                    // @ts-ignore
                    ticks: settings.axis_labels ? false : "",
                    showline: false,
                    zeroline: false,
                    showgrid: false,
                    scaleanchor: "x",
                    scaleratio: 1,
                  },
                }}
                style={{ width: "100%", height: "100%" }}
              />
            )}
          </div>
          <HeatmapSidebar
            settings={settings}
            updateSettings={updateSettings}
            sequences_count={sequences_count}
            colorScales={colorScales}
          />
        </>
      ) : null}
    </>
  );
};
