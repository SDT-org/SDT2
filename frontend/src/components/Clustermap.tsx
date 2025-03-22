import React from "react";
import { type DocState, type SetDocState, useAppState } from "../appState";
import type { ColorScaleArray } from "../colorScales";
import { plotFontMonospace, plotFontSansSerif } from "../constants";
import { formatClustermapData } from "../heatmapUtils";
import { useMetrics, useSize } from "../hooks/heatmap";
import { useHeatmapRenderToggle } from "../hooks/useHeatmapRenderToggle";
import type { HeatmapData } from "../plotTypes";
import { ClustermapSidebar } from "./ClustermapSidebar";
import { D3CanvasHeatmap } from "./D3CanvasHeatmap";
import { D3SvgHeatmap } from "./D3SvgHeatmap";

type ClusterDataItem = {
  id: string;
  group: number;
};

export const Clustermap = ({
  data,
  docState,
  setDocState,
  tickText,
  leftSidebarCollapsed,
}: {
  data: HeatmapData;
  docState: DocState;
  setDocState: SetDocState;
  tickText: string[];
  leftSidebarCollapsed: boolean;
}) => {
  const [clusterData, setClusterData] = React.useState<
    Array<ClusterDataItem> | undefined
  >(undefined);
  const [orderedMatrix, setOrderedMatrix] = React.useState<HeatmapData>(data);
  const [orderedTickText, setOrderedTickText] =
    React.useState<string[]>(tickText);

  React.useEffect(() => {
    // @ts-ignore
    window.pywebview.api
      .get_cluster_ordered_data(
        docState.id,
        docState.clustermap.threshold,
        docState.clustermap.method,
      )
      // @ts-ignore
      .then((response) => {
        console.log("API response:", {
          matrix: response.matrix
            ? `${response.matrix.length}x${response.matrix[0]?.length || 0}`
            : "empty",
          tickText: response.tickText?.length || 0,
          clusterData: response.clusterData?.length || 0,
        });

        setClusterData(response.clusterData);
        setOrderedTickText(response.tickText);
        setOrderedMatrix(response.matrix);
      })
      .catch((error: Error) => {
        console.error("Error getting cluster ordered data:", error);
      });
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [docState.id, docState.clustermap.threshold, docState.clustermap.method]);

  const { appState } = useAppState();
  const elementRef = React.useRef<HTMLDivElement | null>(null);
  const size = useSize(elementRef, leftSidebarCollapsed);
  const { clustermap: settings } = docState;
  const { margin } = useMetrics(settings, orderedTickText);
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

  const clustermapData = React.useMemo(() => {
    console.log("Creating clustermapData with:", {
      matrixDimensions: `${orderedMatrix.length}x${orderedMatrix[0]?.length || 0}`,
      tickTextLength: orderedTickText.length,
      clusterDataLength: clusterData?.length || 0,
    });

    const result = formatClustermapData(
      orderedMatrix,
      orderedTickText,
      clusterData,
    );
    console.log(`Formatted data contains ${result.length} items`);
    return result;
  }, [orderedMatrix, clusterData, orderedTickText]);

  const forceSvgRender = useHeatmapRenderToggle();

  const titleFont =
    settings.titleFont === "Monospace" ? plotFontMonospace : plotFontSansSerif;

  console.log("Rendering D3CanvasHeatmap with:", {
    dataLength: clustermapData?.length || 0,
    tickTextLength: orderedTickText?.length || 0,
    clusterDataLength: clusterData?.length || 0,
  });

  return (
    <>
      <div className="app-main" ref={elementRef}>
        {clusterData && forceSvgRender ? (
          <div className="debug-toast">SVG</div>
        ) : null}
        {clusterData ? (
          forceSvgRender ||
          (appState.showExportModal && appState.saveFormat === "svg") ? (
            <D3SvgHeatmap
              data={clustermapData}
              clusterData={clusterData}
              settings={{ ...docState.heatmap, ...settings }}
              tickText={orderedTickText}
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
          ) : (
            <D3CanvasHeatmap
              data={clustermapData}
              clusterData={clusterData}
              settings={{ ...docState.heatmap, ...settings }}
              tickText={orderedTickText}
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
