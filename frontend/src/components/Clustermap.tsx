import React from "react";
import type { DocState, SetDocState } from "../appState";
import type { ColorScaleArray } from "../colorScales";
import { plotFontMonospace, plotFontSansSerif } from "../constants";
import { useHeatmapData, useMetrics, useSize } from "../hooks/map";
import type { ClustermapSettings, HeatmapData } from "../plotTypes";
import { ClustermapSidebar } from "./ClustermapSidebar";
import { D3CanvasHeatmap } from "./D3CanvasHeatmap";

export type ClustermapRenderProps = {
  data: { x: number; y: number; value: number }[];
  settings: ClustermapSettings;
  tickText: string[];
  width: number;
  height: number;
  roundTo: number;
  axlabel_xfontsize: number;
  axlabel_xrotation: number;
  axlabel_yrotation: number;
  showPercentIdentities: boolean;
  showTitles: boolean;
  title: string;
  showscale: boolean;
  axis_labels: boolean;
  titleFont: typeof plotFontMonospace | typeof plotFontSansSerif;
  margin: { top: number; bottom: number; left: number; right: number };
} & Pick<ClustermapSettings, "annotation_font_size" | "axis_labels"> & {
    clusterData?: {
      id: string;
      group: number;
    }[];
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
  const [clusterData, setClusterData] =
    React.useState<
      {
        id: string;
        group: number;
      }[]
    >();

  React.useEffect(() => {
    window.pywebview.api
      .generate_cluster_data(
        docState.id,
        docState.clustermap.threshold_one,
        docState.clustermap.threshold_two,
      )
      .then((clusterData) => {
        setClusterData(clusterData);
      });
  }, [docState.id, docState.clustermap]);

  const elementRef = React.useRef<HTMLDivElement | null>(null);
  const size = useSize(elementRef, leftSidebarCollapsed);
  const { clustermap: settings } = docState;
  const { margin } = useMetrics(settings, tickText);
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

  const d3HeatmapData = useHeatmapData(data);

  return (
    <>
      <div className="app-main" ref={elementRef} style={{ background: "#fff" }}>
        {data && clusterData ? (
          <D3CanvasHeatmap
            data={d3HeatmapData}
            clusterData={clusterData}
            settings={{ ...docState.heatmap, ...settings }}
            tickText={tickText}
            colorScale={colorScale}
            cbarHeight={0}
            cbarWidth={0}
            minVal={0}
            maxVal={100}
            width={size.width}
            height={size.height}
            roundTo={settings.annotation_rounding}
            showscale={false}
            annotation_font_size={settings.annotation_font_size}
            axlabel_xrotation={settings.axlabel_xrotation}
            axlabel_xfontsize={settings.axlabel_xfontsize}
            axlabel_yrotation={settings.axlabel_yrotation}
            titleFont={
              settings.titleFont === "Monospace"
                ? plotFontMonospace
                : plotFontSansSerif
            }
            showPercentIdentities={settings.annotation}
            showTitles={settings.showTitles}
            title={settings.title}
            axis_labels={settings.axis_labels}
            margin={margin}
            cellSpace={settings.cellspace}
          />
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
