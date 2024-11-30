import React from "react";
import {
  Button,
  type Key,
  Tab,
  TabList,
  TabPanel,
  Tabs,
} from "react-aria-components";
import type { AppState, SetAppState } from "../appState";
import { useDistributionState } from "../distributionState";
import type {
  DistributionData,
  GetDataResponse,
  HeatmapData,
  HeatmapSettings,
} from "../plotTypes";
import { Distribution } from "./Distribution";
import { Heatmap } from "./Heatmap";

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
  const [tickText, setTickText] = React.useState<string[]>([""]);
  const [heatmapData, setHeatmapData] = React.useState<HeatmapData>();
  const [distributionData, setDistributionData] =
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
      cbar_aspect: 25,
      cbar_pad: 10,
      axis_labels: false,
      axlabel_xrotation: 270,
      axlabel_xfontsize: 12,
      axlabel_yrotation: 360,
      axlabel_yfontsize: 12,
      cutoff_1: 95,
      cutoff_2: 75,
    },
  );

  const {
    state: distributionState,
    setState: setDistributionState,
    updateHistogram,
    updateRaincloud,
    updateViolin,
  } = useDistributionState();

  const getData = () => {
    setLoading(true);

    window.pywebview.api
      .get_data()
      .then((rawData) => {
        const {
          data,
          metadata,
          identity_scores,
          gc_stats,
          length_stats,
        }: GetDataResponse = JSON.parse(rawData.replace(/\bNaN\b/g, "null"));

        const [tickText, ...parsedData] = data;

        setTickText(tickText as string[]);
        setHeatmapData(parsedData);
        setDistributionData({
          metadata,
          gc_stats,
          length_stats,
          raw_mat: identity_scores.map((i) => i[2]),
          identity_combos: identity_scores.map((i) => [i[0], i[1]]),
        });

        updateHeatmapState({
          vmin: metadata.minVal,
          ...(appState.sequences_count > 99 && {
            axis_labels: false,
            cellspace: 0,
          }),
        });
      })
      .catch((e) => {
        setLoading(false);
        console.error(e);
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
            {appState.sequences_count > 0 ? (
              <>
                {appState.sequences_count} Sequence
                {appState.sequences_count === 1 ? "" : "s"}
              </>
            ) : null}
            <span className="filename">{appState.basename}</span>
          </div>
        </div>
        <TabList>
          <Tab id="heatmap">Heatmap</Tab>
          <Tab id="plot">Distribution</Tab>
        </TabList>
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
            tickText={tickText}
          />
        ) : null}
      </TabPanel>
      <TabPanel id="plot" className="app-panel">
        {distributionData ? (
          <Distribution
            data={distributionData}
            state={distributionState}
            setState={setDistributionState}
            updateHistogram={updateHistogram}
            updateRaincloud={updateRaincloud}
            updateViolin={updateViolin}
          />
        ) : null}
      </TabPanel>
      {loading ? <div className="api-loader"></div> : null}
    </Tabs>
  );
};
