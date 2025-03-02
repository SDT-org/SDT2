import React from "react";
import type { DocState } from "../appState";
import type { ColorScaleArray } from "../colorScales";
import { plotFontMonospace } from "../constants";
import { useHeatmapData, useMetrics, useSize } from "../hooks/map";
import type { HeatmapData } from "../plotTypes";
import { D3CanvasHeatmap } from "./D3CanvasHeatmap";

export const Clustermap = ({
  data,
  docState,
  tickText,
  leftSidebarCollapsed,
}: {
  data: HeatmapData;
  docState: DocState;
  tickText: string[];
  leftSidebarCollapsed: boolean;
}) => {
  const [clusterData, setClusterData] =
    React.useState<
      {
        id: string;
        group: number;
      }[]
    >();

  React.useEffect(() => {
    window.pywebview.api
      .generate_cluster_data(docState.id, 75, 0)
      .then((clusterData) => {
        setClusterData(clusterData);
      })
      .catch((e) => {
        alert(e);
        throw e;
      });
  }, [docState.id]);

  const elementRef = React.useRef<HTMLDivElement | null>(null);
  const size = useSize(elementRef, leftSidebarCollapsed);
  // TODO: make this cluster settings
  const { heatmap: settings } = docState;
  const { cbar_aspect, cbar_shrink, margin } = useMetrics(settings, tickText);

  const colorScale: ColorScaleArray = [
    [0, "rgb(0,0,0)"],
    [1, "rgb(235,235,235)"],
  ];

  const d3HeatmapData = useHeatmapData(data);

  return data && clusterData ? (
    <div className="app-main" ref={elementRef} style={{ background: "#fff" }}>
      <D3CanvasHeatmap
        data={d3HeatmapData}
        clusterData={clusterData}
        settings={settings}
        tickText={tickText}
        colorScale={colorScale}
        minVal={settings.vmin}
        maxVal={settings.vmax}
        width={size.width}
        height={size.height}
        cellSpace={settings.cellspace}
        roundTo={settings.annotation_rounding}
        showscale={false}
        cbarHeight={cbar_shrink}
        cbarWidth={cbar_aspect}
        annotation_font_size={settings.annotation_font_size}
        axlabel_xrotation={settings.axlabel_xrotation}
        axlabel_xfontsize={settings.axlabel_xfontsize}
        axlabel_yrotation={settings.axlabel_yrotation}
        axlabel_yfontsize={settings.axlabel_yfontsize}
        titleFont={plotFontMonospace}
        showPercentIdentities={settings.annotation}
        showTitles={settings.showTitles}
        title={settings.title}
        subtitle={settings.subtitle}
        axis_labels={settings.axis_labels}
        margin={margin}
      />
    </div>
  ) : null;
};
