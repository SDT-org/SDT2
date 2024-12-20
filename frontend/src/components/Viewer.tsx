import type React from "react";
import {
  Button,
  type Key,
  Tab,
  TabList,
  TabPanel,
  Tabs,
} from "react-aria-components";
import type { AppState, SetAppState } from "../appState";
import { useGetData } from "../hooks/useGetData";
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
  const { loading, tickText, heatmapData, distributionData, metaData } =
    useGetData();

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
