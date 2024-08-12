import React from "react";
import { AppState, SetAppState } from "../appState";
import { initialAppState } from "../appState";
import {
  HeatmapData,
  HeatmapSettings,
  HistogramData,
  HistogramSettings,
} from "../plotTypes";
import { Heatmap } from "./Heatmap";
import { Histogram } from "./Histogram";
import { ExportData } from "./ExportData";
import { Tab, TabList, TabPanel, Tabs } from "react-aria-components";

export const Viewer = ({
  appState,
  setAppState,
}: {
  appState: AppState;
  setAppState: SetAppState;
}) => {
  const [loading, setLoading] = React.useState(false);
  const [heatmapTickText, setHeatmapTickText] = React.useState<string[]>([""]);
  const [heatmapData, setHeatmapData] = React.useState<HeatmapData>();
  const [histogramData, setHisogramData] = React.useState<HistogramData>();

  const [heatmapSettings, setHeatmapSettings] = React.useState<HeatmapSettings>(
    {
      colorscale: "Portland",
      reverse: false,
      vmax: 100,
      vmin: 65,
      cellspace: 0,
      annotation: false,
      annotation_font_size: 10,
      annotation_rounding: 0,
      annotation_alpha: "0",
      color: "white",
      showscale: true,
      cbar_shrink: 1,
      cbar_aspect: "25",
      cbar_pad: "10",
      axis_labels: false,
      axlabel_xrotation: 360  ,
      axlabel_xfontsize: 12,
      axlabel_yrotation:270,
      axlabel_yfontsize: 12,
    },
  );

  const [histogramSettings, setHistogramSettings] =
    React.useState<HistogramSettings>({
      lineColor: "tomato",
      lineWidth: 3,
      lineShape: "linear",
      barlineColor: "tomato",
      barOutlineWidth: 1,
      barColor: "lightblue",
      showHistogram: true,
      showLinePlot: false,
      showScatterPlot: false,
      markerSymbol: "square",
      markerColor: "tomato",
      markerSize: 7,
      showGrid: true,
      showLine: true,
      showZeroLine: true,
      plotTitle: "Distribution of Percent Identities",
      showTickLabels: true,
      showAxisLabels: true,
    });

  const getData = () => {
    setLoading(true);
    Promise.all([
      window.pywebview.api.get_heatmap_data().then((rawData) => {
        const {
          metadata,
          data,
        }: {
          metadata: {
            minVal: number;
            maxVal: number;
          };
          data: string[][];
        } = JSON.parse(rawData.replace(/\bNaN\b/g, "null"));
        const [tickText, ...parsedData] = data;
        setHeatmapTickText(tickText as string[]);
        setHeatmapData(parsedData);
        updateHeatmapState({
          vmin: metadata.minVal,
          ...(appState.sequences_count > 99 && {
            axis_labels: false,
            cellspace: 0,
          }),
        });
      }),
      window.pywebview.api.get_line_histo_data().then((data) => {
        const histogramData = JSON.parse(data.replace(/\bNaN\b/g, "null"));
        setHisogramData(histogramData);
      }),
    ]).then(() => {
      setLoading(false);
    });
  };

  React.useEffect(() => {
    getData();
  }, []);

  const updateHeatmapState = (newState: Partial<HeatmapSettings>) => {
    setHeatmapSettings((previous) => {
      return {
        ...previous,
        ...newState,
      };
    });
  };
  const updateHistogramSettings = (newState: Partial<HistogramSettings>) => {
    setHistogramSettings((previous) => {
      return {
        ...previous,
        ...newState,
      };
    });
  };

  return (
    <Tabs className="app-wrapper with-header">
      <div className="app-header">
        <div className="left">
          <button
            type="button"
            aria-label="Start a new sequence run"
            onClick={() => {
              window.pywebview.api.reset_state().then(() => {
                // Reset client state as well
                setAppState((previous) => {
                  return {
                    ...initialAppState,
                    client: { ...previous.client },
                  };
                });
              });
            }}
          >
            New
          </button>
          <ExportData />
          <div className="run-info">
            <strong>
              {appState.sequences_count > 0 ? (
                <>
                  {appState.sequences_count} Sequence
                  {appState.sequences_count === 1 ? "" : "s"}
                </>
              ) : null}
            </strong>
            <span className="filename">{appState.basename}</span>
          </div>
        </div>
        <TabList className="data-view">
          <Tab id="heatmap">Heatmap</Tab>
          <Tab id="plot">Distribution</Tab>
        </TabList>
      </div>

      <TabPanel id="heatmap" className="app-panel">
        <Heatmap
          data={heatmapData}
          settings={heatmapSettings}
          updateSettings={updateHeatmapState}
          tickText={heatmapTickText}
        />
      </TabPanel>
      <TabPanel id="plot" className="app-panel">
        <Histogram
          data={histogramData}
          settings={histogramSettings}
          updateSettings={updateHistogramSettings}
        />
      </TabPanel>
      {loading ? <div className="api-loader"></div> : null}
    </Tabs>
  );
};
