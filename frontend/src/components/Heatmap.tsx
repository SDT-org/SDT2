import React from "react";
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
import type { HeatmapData, HeatmapSettings } from "../plotTypes";
import { D3CanvasHeatmap } from "./D3CanvasHeatmap";
import { D3Heatmap } from "./D3Heatmap";
import { HeatmapSidebar } from "./HeatmapSidebar";

export type HeatmapRenderProps = {
  data: { x: number; y: number; value: number }[];
  tickText: string[];
  colorScale: ColorScaleArray;
  minVal: number;
  maxVal: number;
  width: number;
  height: number;
  cellSpace: number;
  roundTo: number;
  cbarWidth: number;
  cbarHeight: number;
  axlabel_xfontsize: number;
  axlabel_yfontsize: number;
  axlabel_xrotation: number;
  axlabel_yrotation: number;
  showPercentIdentities: boolean;
  showTitles: boolean;
  title: string;
  subtitle: string;
  showscale: boolean;
  axis_labels: boolean;
  titleFont: typeof plotFontMonospace | typeof plotFontSansSerif;
  margin: { top: number; bottom: number; left: number; right: number };
} & Pick<HeatmapSettings, "annotation_font_size" | "axis_labels">;

export const Heatmap = ({
  data,
  tickText,
  docState,
  setDocState,
  leftSidebarCollapsed,
}: {
  docState: DocState;
  setDocState: SetDocState;
  updateDocState: UpdateDocState;
  data: HeatmapData;
  tickText: string[];
  leftSidebarCollapsed: boolean;
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

  // biome-ignore lint/correctness/useExhaustiveDependencies(leftSidebarCollapsed): trigger updateSize
  React.useEffect(() => {
    updateSize();

    const handleResize = () => {
      updateSize();
    };

    window.addEventListener("resize", handleResize);

    return () => {
      window.removeEventListener("resize", handleResize);
    };
  }, [updateSize, leftSidebarCollapsed]);

  const cbar_shrink = settings.cbar_shrink * 60;
  const cbar_aspect = settings.cbar_aspect * 10;

  let colorScale = colorScales[settings.colorScaleKey];

  if (settings.reverse) {
    colorScale = [...colorScale]
      .reverse()
      .map((data, i) => [
        (colorScale[i] ?? colorScale[0])[0],
        data[1],
      ]) as ColorScaleArray;
  }

  const titleFont =
    settings.titleFont === "Monospace" ? plotFontMonospace : plotFontSansSerif;

  const longestTickWidth =
    Math.max(...tickText.map((tick) => tick.length)) *
    (settings.axlabel_yfontsize * 1.25);

  const margin = {
    top: 60,
    right: 60,
    bottom: longestTickWidth,
    left: longestTickWidth,
  };

  return (
    <>
      {data ? (
        <>
          <div
            className="app-main"
            ref={elementRef}
            style={{ background: "#fff" }}
          >
            {appState.showExportModal && appState.saveFormat === "svg" ? (
              <D3Heatmap
                data={d3HeatmapData}
                tickText={tickText}
                colorScale={colorScale}
                minVal={settings.vmin}
                maxVal={settings.vmax}
                width={size.width}
                height={size.height}
                cellSpace={settings.cellspace}
                showPercentIdentities={settings.annotation}
                roundTo={settings.annotation_rounding}
                cbarHeight={cbar_shrink}
                cbarWidth={cbar_aspect}
                annotation_font_size={settings.annotation_font_size}
                axlabel_xrotation={settings.axlabel_xrotation}
                axlabel_xfontsize={settings.axlabel_xfontsize}
                axlabel_yrotation={settings.axlabel_yrotation}
                axlabel_yfontsize={settings.axlabel_yfontsize}
                titleFont={titleFont}
                showTitles={settings.showTitles}
                title={settings.title}
                subtitle={settings.subtitle}
                axis_labels={settings.axis_labels}
                showscale={settings.showscale}
                margin={margin}
              />
            ) : (
              <D3CanvasHeatmap
                data={d3HeatmapData}
                tickText={tickText}
                colorScale={colorScale}
                minVal={settings.vmin}
                maxVal={settings.vmax}
                width={size.width}
                height={size.height}
                cellSpace={settings.cellspace}
                roundTo={settings.annotation_rounding}
                showscale={settings.showscale}
                cbarHeight={cbar_shrink}
                cbarWidth={cbar_aspect}
                annotation_font_size={settings.annotation_font_size}
                axlabel_xrotation={settings.axlabel_xrotation}
                axlabel_xfontsize={settings.axlabel_xfontsize}
                axlabel_yrotation={settings.axlabel_yrotation}
                axlabel_yfontsize={settings.axlabel_yfontsize}
                titleFont={titleFont}
                showPercentIdentities={settings.annotation}
                showTitles={settings.showTitles}
                title={settings.title}
                subtitle={settings.subtitle}
                axis_labels={settings.axis_labels}
                margin={margin}
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
