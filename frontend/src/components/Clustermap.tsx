import React from "react";
import { type DocState, type SetDocState, useAppState } from "../appState";
import type { ColorScaleArray } from "../colorScales";
import { plotFontMonospace, plotFontSansSerif } from "../constants";
import { formatClustermapData } from "../heatmapUtils";
import { useMetrics, useSize } from "../hooks/heatmap";
import { useHeatmapRenderToggle } from "../hooks/useHeatmapRenderToggle";
import type {
  ClusterDataItem,
  ClusterStats,
  GetClustermapDataResponse,
  HeatmapData,
  MetaData,
} from "../plotTypes";
import { ClusterStatsDisplay } from "./ClusterStats";
import { ClustermapSidebar } from "./ClustermapSidebar";
import { D3CanvasHeatmap } from "./D3CanvasHeatmap";
import { D3SvgHeatmap } from "./D3SvgHeatmap";

export const Clustermap = ({
  data,
  metaData,
  docState,
  setDocState,
  tickText,
  leftSidebarCollapsed,
}: {
  data: HeatmapData;
  metaData: MetaData;
  docState: DocState;
  setDocState: SetDocState;
  tickText: string[];
  leftSidebarCollapsed: boolean;
}) => {
  const [clusterData, setClusterData] = React.useState<
    Array<ClusterDataItem> | undefined
  >(undefined);
  const [clusterStats, setClusterStats] = React.useState<
    ClusterStats | undefined
  >(undefined);
  const [orderedMatrix, setOrderedMatrix] = React.useState<HeatmapData>(data);
  const [orderedIds, setOrderedIds] = React.useState<string[]>(tickText);

  React.useEffect(() => {
    const start = performance.now();
    loaderRef.current?.setAttribute("data-hidden", "false");
    window.pywebview.api
      .get_clustermap_data(
        docState.id,
        docState.clustermap.threshold,
        docState.clustermap.method,
      )
      .then((response: GetClustermapDataResponse) => {
        setClusterData(response.clusterData);
        setOrderedIds(response.tickText);
        setOrderedMatrix(response.matrix);
        setClusterStats(response.cluster_stats);
        const end = performance.now();
        console.log(`Clustermap fetch time: ${end - start}ms`);
      });
  }, [docState.id, docState.clustermap.threshold, docState.clustermap.method]);

  const { appState } = useAppState();
  const elementRef = React.useRef<HTMLDivElement | null>(null);
  const size = useSize(elementRef, leftSidebarCollapsed);
  const { clustermap: settings } = docState;
  const { margin } = useMetrics(settings, orderedIds);
  const updateSettings = React.useCallback(
    (values: Partial<DocState["clustermap"]>) =>
      setDocState((prev) => ({
        ...prev,
        clustermap: {
          ...prev.clustermap,
          ...values,
        },
      })),
    [setDocState],
  );

  const colorScale: ColorScaleArray = [
    [0, "rgb(245,245,245)"],
    [1, "rgb(245,245,245)"],
  ];

  const clustermapData = React.useMemo(
    () => formatClustermapData(orderedMatrix, orderedIds, clusterData),
    [orderedMatrix, orderedIds, clusterData],
  );

  const forceSvgRender = useHeatmapRenderToggle();
  const renderSvg =
    forceSvgRender ||
    (appState.exportStatus === "exporting" && appState.saveFormat === "svg");

  const titleFont =
    settings.titleFont === "Monospace" ? plotFontMonospace : plotFontSansSerif;

  const loaderRef = React.useRef<HTMLDivElement | null>(null);
  const onRenderComplete = React.useCallback(() => {
    loaderRef.current?.setAttribute("data-hidden", "true");
  }, []);

  return (
    <>
      <div
        className="app-main"
        ref={elementRef}
        style={{ position: "relative" }}
      >
        {clusterStats && (
          <ClusterStatsDisplay
            metaData={metaData}
            clusterStats={clusterStats}
          />
        )}
        <div
          className="app-overlay app-loader app-main-loader delay"
          aria-hidden="true"
          data-hidden="false"
          ref={loaderRef}
        />
        {clusterData && forceSvgRender ? (
          <div className="debug-toast">SVG</div>
        ) : null}
        {clusterData ? (
          renderSvg ? (
            <>
              <div className="debug-toast">SVG</div>
              <D3SvgHeatmap
                data={clustermapData}
                clusterData={clusterData}
                settings={{ ...docState.heatmap, ...settings }}
                tickText={orderedIds}
                colorScale={colorScale}
                minVal={0}
                maxVal={100}
                width={size.width}
                height={size.height}
                cellSpace={settings.cellspace}
                showPercentIdentities={settings.annotation}
                roundTo={2}
                cbarHeight={0}
                cbarWidth={0}
                axlabel_xrotation={settings.axlabel_xrotation}
                axlabel_fontsize={settings.axlabel_fontsize}
                axlabel_yrotation={settings.axlabel_yrotation}
                titleFont={titleFont}
                showTitles={settings.showTitles}
                title={settings.title}
                axis_labels={settings.axis_labels}
                showscale={false}
                margin={margin}
                showLegend={settings.showLegend}
              />
            </>
          ) : (
            <D3CanvasHeatmap
              data={clustermapData}
              clusterData={clusterData}
              settings={{ ...docState.heatmap, ...settings }}
              tickText={orderedIds}
              colorScale={colorScale}
              cbarHeight={0}
              cbarWidth={0}
              minVal={0}
              maxVal={100}
              width={size.width}
              height={size.height}
              roundTo={2}
              showscale={false}
              axlabel_xrotation={settings.axlabel_xrotation}
              axlabel_fontsize={settings.axlabel_fontsize}
              axlabel_yrotation={settings.axlabel_yrotation}
              titleFont={titleFont}
              showPercentIdentities={settings.annotation}
              showTitles={settings.showTitles}
              title={settings.title}
              axis_labels={settings.axis_labels}
              margin={margin}
              cellSpace={settings.cellspace}
              showLegend={settings.showLegend}
              onRenderComplete={onRenderComplete}
            />
          )
        ) : null}
      </div>
      <ClustermapSidebar
        settings={settings}
        updateSettings={updateSettings}
        sequences_count={docState.sequences_count}
      />
    </>
  );
};
