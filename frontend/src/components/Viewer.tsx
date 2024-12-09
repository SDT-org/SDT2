import React from "react";
import {
  Button,
  type Key,
  Tab,
  TabList,
  TabPanel,
  Tabs,
} from "react-aria-components";
import { type AppState, type SetAppState, initialAppState } from "../appState";
import type {
  DistributionData,
  GetDataResponse,
  HeatmapData,
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
  const [metaData, setMetaData] = React.useState<
    GetDataResponse["metadata"] | undefined
  >();

  const getScaledFontSize = React.useCallback(
    (base: number, count: number) =>
      count < 20 ? 16 : Math.min(16, Math.max(5, base / (1 + count / 90))),
    [],
  );

  const getData = React.useCallback(() => {
    setLoading(true);

    window.pywebview.api
      .get_data()
      .then((rawData) => {
        const parsedResponse: GetDataResponse = JSON.parse(
          rawData.replace(/\bNaN\b/g, "null"),
        );
        const { data, metadata, identity_scores, gc_stats, length_stats } =
          parsedResponse;

        const [tickText, ...parsedData] = data;

        setTickText(tickText as string[]);
        setHeatmapData(parsedData);
        setDistributionData({
          gc_stats,
          length_stats,
          raw_mat: identity_scores.map((i) => i[2]),
          identity_combos: identity_scores.map((i) => [i[0], i[1]]),
        });
        setMetaData(metadata);

        setAppState((prev) => ({
          ...prev,
          client: {
            ...prev.client,
            heatmap: {
              ...prev.client.heatmap,
              vmin: metadata.minVal,
              ...(appState.sequences_count > 99 && {
                annotation: false,
                axis_labels: false,
                cellspace: 0,
              }),
              annotation_font_size: getScaledFontSize(
                initialAppState.client.heatmap.annotation_font_size,
                appState.sequences_count,
              ),
            },
          },
        }));
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
  }, [appState.sequences_count, setAppState, getScaledFontSize]);

  React.useEffect(() => {
    getData();
  }, [getData]);

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
    <div className="app-wrapper with-header">
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
        <Tabs
          selectedKey={appState.client.dataView}
          onSelectionChange={setDataView}
        >
          <TabList>
            <Tab id="heatmap">Heatmap</Tab>
            <Tab id="distribution">Distribution</Tab>
          </TabList>
        </Tabs>
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
      <Tabs
        className={"app-panels"}
        selectedKey={appState.client.dataView}
        onSelectionChange={setDataView}
      >
        <TabPanel id="heatmap" className="app-panel">
          {appState.client.dataView === "heatmap" &&
          !loading &&
          heatmapData &&
          metaData ? (
            <Heatmap data={heatmapData} tickText={tickText} />
          ) : null}
        </TabPanel>
        <TabPanel id="distribution" className="app-panel">
          {distributionData && metaData ? (
            <Distribution data={distributionData} metaData={metaData} />
          ) : null}
        </TabPanel>
      </Tabs>
      {loading ? <div className="app-overlay app-loader" /> : null}
    </div>
  );
};
