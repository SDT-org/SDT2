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
import { useHeatmapData, useMetrics, useSize } from "../hooks/map";
import type { HeatmapData, HeatmapSettings, MetaData } from "../plotTypes";
import { D3CanvasHeatmap } from "./D3CanvasHeatmap";
import { D3Heatmap } from "./D3Heatmap";
import { HeatmapSidebar } from "./HeatmapSidebar";

export type HeatmapRenderProps = {
  // TODO: just use settings
  data: { x: number; y: number; value: number }[];
  settings: HeatmapSettings;
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
  metaData,
  tickText,
  docState,
  setDocState,
  leftSidebarCollapsed,
}: {
  docState: DocState;
  setDocState: SetDocState;
  updateDocState: UpdateDocState;
  data: HeatmapData;
  metaData: MetaData;
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

  console.log(data);
  console.log(tickText);

  const discreteColorScale: ColorScaleArray = React.useMemo(() => {
    const scales = [
      [metaData.minVal, "#CDF0FF"],
      [settings.cutoff_2 - 1, "#20B9FF"],
      [settings.cutoff_2, "#C3E8D3"],
      [settings.cutoff_1, "#009942"],
      [Math.min(100, settings.cutoff_1 + 1), "#FFDCDD"],
    ];

    if (settings.cutoff_1 < 100) {
      scales.push([100, "#FF6167"]);
    }
    return scales as ColorScaleArray;
  }, [settings.cutoff_1, settings.cutoff_2, metaData.minVal]);

  const colorScales = React.useMemo(
    () => ({
      Discrete: discreteColorScale,
      ...defaultColorScales,
    }),
    [discreteColorScale],
  );

  const d3HeatmapData = useHeatmapData(data);

  const { cbar_shrink, cbar_aspect, margin } = useMetrics(settings, tickText);

  const elementRef = React.useRef<HTMLDivElement | null>(null);
  const size = useSize(elementRef, leftSidebarCollapsed);

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
                settings={settings}
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
                settings={settings}
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
