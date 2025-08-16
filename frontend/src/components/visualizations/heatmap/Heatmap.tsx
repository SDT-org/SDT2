import React from "react";
import useAppState, {
  type DocState,
  type SetDocState,
  type UpdateDocState,
} from "../../../appState";
import {
  type ColorScaleArray,
  colorScales as defaultColorScales,
} from "../../../colorScales";
import { plotFontMonospace, plotFontSansSerif } from "../../../constants";
import { useMetrics, useSize } from "../../../hooks/heatmap";
import { useHeatmapRenderToggle } from "../../../hooks/useHeatmapRenderToggle";
import type {
  ClusterDataItem,
  HeatmapData,
  HeatmapSettings,
  MetaData,
} from "../../../plotTypes";
import { AlignmentStats } from "../../stats/AlignmentStats";
import { D3CanvasHeatmap } from "./D3CanvasHeatmap";
import { D3SvgHeatmap } from "./D3SvgHeatmap";
import { HeatmapSidebar } from "./HeatmapSidebar";
import { type formatClustermapData, formatHeatmapData } from "./heatmapUtils";

export type HeatmapRenderProps = {
  // TODO: just use settings
  data:
    | ReturnType<typeof formatHeatmapData>
    | ReturnType<typeof formatClustermapData>;
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
  axlabel_fontsize: number;
  axlabel_xrotation: number;
  axlabel_yrotation: number;
  showPercentIdentities: boolean;
  showTitles: boolean;
  title: string;
  showscale: boolean;
  axis_labels: boolean;
  titleFont: typeof plotFontMonospace | typeof plotFontSansSerif;
  margin: { top: number; bottom: number; left: number; right: number };
  clusterData?: ClusterDataItem[];
  showLegend?: boolean;
  showClusterCounts?: boolean;
  clusterCounts?: { [key: number]: number };
  onRenderComplete?: () => void;
  metaData?: MetaData;
} & Pick<HeatmapSettings, "axis_labels">;

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

  const discreteColorScale: ColorScaleArray = React.useMemo(() => {
    const scales = [
      [settings.vmin, "#CDF0FF"],
      [settings.cutoff_2, "#20B9FF"],
      [settings.cutoff_2, "#C3E8D3"],
      [settings.cutoff_1, "#009942"],
      [Math.min(100, settings.cutoff_1), "#FFDCDD"],
    ];

    if (settings.cutoff_1 < 100) {
      scales.push([100, "#FF6167"]);
    }
    return scales as ColorScaleArray;
  }, [settings.cutoff_1, settings.cutoff_2, settings.vmin]);

  const colorScales = React.useMemo(
    () => ({
      Discrete: discreteColorScale,
      ...defaultColorScales,
    }),
    [discreteColorScale],
  );

  const colorScale = React.useMemo(() => {
    let scale = colorScales[settings.colorScaleKey];
    if (settings.reverse) {
      scale = [...scale]
        .reverse()
        .map((data, i) => [
          (scale[i] ?? scale[0])[0],
          data[1],
        ]) as ColorScaleArray;
    }
    return scale;
  }, [colorScales, settings.colorScaleKey, settings.reverse]);

  const heatmapData = React.useMemo(() => {
    const startTime = appState.debug ? performance.now() : 0;
    const result = formatHeatmapData(data, settings, colorScale);
    if (appState.debug) {
      const endTime = performance.now();
      console.debug(
        `[PERF] formatHeatmapData completed in ${(endTime - startTime).toFixed(2)}ms`,
      );
    }
    return result;
  }, [data, settings, colorScale, appState.debug]);

  const { cbar_shrink, cbar_aspect, margin } = useMetrics(settings, tickText);

  const elementRef = React.useRef<HTMLDivElement | null>(null);
  const size = useSize(elementRef, leftSidebarCollapsed);

  const titleFont =
    settings.titleFont === "Monospace" ? plotFontMonospace : plotFontSansSerif;

  const forceSvgRender = useHeatmapRenderToggle();

  const loaderRef = React.useRef<HTMLDivElement | null>(null);
  const onRenderComplete = React.useCallback(() => {
    loaderRef.current?.setAttribute("data-hidden", "true");
  }, []);

  const renderSvg =
    forceSvgRender ||
    (appState.exportStatus === "exporting" && appState.saveFormat === "svg");

  return (
    <>
      {data ? (
        <>
          <div
            className="app-main"
            ref={elementRef}
            style={{ position: "relative" }}
          >
            {!renderSvg && (
              <AlignmentStats
                metaData={metaData}
                dataLength={metaData.distribution_stats?.count || 0}
                activeDataSet="scores"
                heatmapSettings={settings}
                rawData={data}
              />
            )}
            <div
              className="app-loader app-main-loader delay"
              aria-hidden="true"
              data-hidden="false"
              ref={loaderRef}
            />
            {renderSvg ? (
              <>
                <div className="debug-toast">SVG</div>
                <D3SvgHeatmap
                  data={heatmapData}
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
                  cbarHeight={cbar_shrink ?? settings.cbar_shrink}
                  cbarWidth={cbar_aspect ?? settings.cbar_aspect}
                  axlabel_xrotation={settings.axlabel_xrotation}
                  axlabel_fontsize={settings.axlabel_fontsize}
                  axlabel_yrotation={settings.axlabel_yrotation}
                  titleFont={titleFont}
                  showTitles={settings.showTitles}
                  title={settings.title}
                  axis_labels={settings.axis_labels}
                  showscale={settings.showscale}
                  margin={margin}
                  metaData={metaData}
                />
              </>
            ) : (
              <D3CanvasHeatmap
                data={heatmapData}
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
                cbarHeight={cbar_shrink ?? settings.cbar_shrink}
                cbarWidth={cbar_aspect ?? settings.cbar_aspect}
                axlabel_xrotation={settings.axlabel_xrotation}
                axlabel_fontsize={settings.axlabel_fontsize}
                axlabel_yrotation={settings.axlabel_yrotation}
                titleFont={titleFont}
                showPercentIdentities={settings.annotation}
                showTitles={settings.showTitles}
                title={settings.title}
                axis_labels={settings.axis_labels}
                margin={margin}
                onRenderComplete={onRenderComplete}
                metaData={metaData}
              />
            )}
          </div>
          <HeatmapSidebar
            settings={settings}
            updateSettings={updateSettings}
            sequences_count={sequences_count}
            continuousColorScales={defaultColorScales}
          />
        </>
      ) : null}
    </>
  );
};
