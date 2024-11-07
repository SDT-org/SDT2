import React from "react";
import { AppState, SetAppState } from "../appState";
import { HeatmapData, HeatmapSettings, DistributionData } from "../plotTypes";
import { Heatmap } from "./Heatmap";
import { Distribution } from "./Distribution";
import {
  Button,
  Key,
  Tab,
  TabList,
  TabPanel,
  Tabs,
} from "react-aria-components";

export const Viewer = ({
  appState,
  setAppState,
  mainMenu,
}: {
  appState: AppState;
  setAppState: SetAppState;
  mainMenu: React.ReactNode;
}) => {
  const [loading, setLoading] = React.useState(false);
  const [heatmapTickText, setHeatmapTickText] = React.useState<string[]>([""]);
  const [heatmapData, setHeatmapData] = React.useState<HeatmapData>();
  const [DistributionData, setHisogramData] =
    React.useState<DistributionData>();

  const [heatmapSettings, setHeatmapSettings] = React.useState<HeatmapSettings>(
    {
      colorscale: "Portland",
      reverse: false,
      vmax: 100,
      vmin: 65,
      cellspace: 1,
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
      axlabel_xrotation: 270,
      axlabel_xfontsize: 12,
      axlabel_yrotation: 360,
      axlabel_yfontsize: 12,
      cutoff_1: 95,
      cutoff_2: 75,
    },
  );

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
      window.pywebview.api.get_distribution_data().then((data) => {
        const DistributionData = JSON.parse(data.replace(/\bNaN\b/g, "null"));
        setHisogramData(DistributionData);
      }),
    ])
      .catch(() => {
        setLoading(false);
        alert(
          "An error occured while processing this file. Please ensure it is a valid, SDT-compatible file.",
        );
        window.pywebview.api.reset_state();
      })
      .finally(() => {
        setLoading(false);
      });
  };

  React.useEffect(() => {
    getData();
  }, [appState.basename]);

  const updateHeatmapState = (newState: Partial<HeatmapSettings>) => {
    setHeatmapSettings((previous) => {
      return {
        ...previous,
        ...newState,
      };
    });
  };

  const setDataView = (newValue: Key) =>
    setAppState((previous) => {
      return {
        ...previous,
        client: {
          ...previous.client,
          dataView: newValue as AppState["client"]["dataView"],
        },
      };
    });

  return (
    <Tabs
      className="app-wrapper with-header"
      selectedKey={appState.client.dataView}
      onSelectionChange={setDataView}
    >
      <div className="app-header">
        <div className="left">
          {mainMenu}
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
        <div>
          <TabList className="data-view">
            <Tab id="heatmap">Heatmap</Tab>
            <Tab id="plot">Distribution</Tab>
          </TabList>
        </div>
        <div className="right">
          <Button
            onPress={() =>
              setAppState((previous) => {
                return {
                  ...previous,
                  client: { ...previous.client, showExportModal: true },
                };
              })
            }
          >
            Export
          </Button>
        </div>
      </div>

      <TabPanel id="heatmap" className="app-panel">
        {heatmapData ? (
          <Heatmap
            data={heatmapData}
            settings={heatmapSettings}
            updateSettings={updateHeatmapState}
            tickText={heatmapTickText}
          />
        ) : null}
      </TabPanel>
      <TabPanel id="plot" className="app-panel">
        {DistributionData ? <Distribution data={DistributionData} /> : null}
      </TabPanel>
      {loading ? <div className="api-loader"></div> : null}
    </Tabs>
  );
};
