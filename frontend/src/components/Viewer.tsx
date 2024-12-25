import React from "react";
import {
  Button,
  type Key,
  Tab,
  TabList,
  TabPanel,
  Tabs,
} from "react-aria-components";
import {
  type DocState,
  type SetDocState,
  type UpdateDocState,
  useAppState,
} from "../appState";
import { useGetData } from "../hooks/useGetData";
import { Distribution } from "./Distribution";
import { Heatmap } from "./Heatmap";
import { SelectDocumentMenu } from "./SelectDocumentMenu";

export const Viewer = ({
  docState,
  setDocState,
  updateDocState,
  mainMenu,
  tabView,
}: {
  docState: DocState;
  setDocState: SetDocState;
  updateDocState: UpdateDocState;
  mainMenu: React.ReactNode;
  tabView: "tabs" | "select";
}) => {
  const { setAppState } = useAppState();
  const { loading, tickText, heatmapData, distributionData, metaData } =
    useGetData(docState, setDocState);

  const setDataView = React.useCallback(
    (newValue: Key) => {
      updateDocState({
        dataView: newValue as DocState["dataView"],
      });
    },
    [updateDocState],
  );

  const header = (
    <Tabs
      selectedKey={docState.dataView}
      onSelectionChange={setDataView}
      orientation="vertical"
    >
      <TabList>
        <Tab id="heatmap">
          <div style={{ display: "flex", alignItems: "center", gap: "0.8rem" }}>
            <svg
              height={12}
              xmlns="http://www.w3.org/2000/svg"
              viewBox="0 0 30 30"
              aria-hidden={true}
              color="currentcolor"
            >
              <rect id="one" x="0" y="0" width={8} height={8} />
              <rect x="0" y="10" width={8} height={8} />
              <rect x="10" y="10" width={8} height={8} />
              <rect x="0" y="20" width={8} height={8} />
              <rect x="10" y="20" width={8} height={8} />
              <rect x="20" y="20" width={8} height={8} />
            </svg>
            Heatmap
          </div>
        </Tab>
        <Tab id="distribution">
          <div style={{ display: "flex", alignItems: "center", gap: "0.8rem" }}>
            <svg
              height={12}
              xmlns="http://www.w3.org/2000/svg"
              viewBox="0 0 30 30"
              aria-hidden={true}
              color="currentcolor"
            >
              <rect x="0" y="18" width={6} height={12} fill="black" />
              <rect x="10" y="0" width={6} height={30} fill="black" />
              <rect x="20" y="10" width={6} height={20} fill="black" />
            </svg>
            Distribution
          </div>
        </Tab>
      </TabList>
    </Tabs>
  );

  const footer = (
    <Button
      onPress={() =>
        setAppState((prev) => ({
          ...prev,
          showExportModal: true,
        }))
      }
    >
      Export
    </Button>
  );

  return (
    <Tabs selectedKey={docState.dataView} onSelectionChange={setDataView}>
      <div className="app-wrapper">
        <TabPanel id="heatmap" className="app-panel">
          {docState.dataView === "heatmap" &&
          !loading &&
          heatmapData &&
          metaData ? (
            <Heatmap
              data={heatmapData}
              tickText={tickText}
              docState={docState}
              setDocState={setDocState}
              updateDocState={updateDocState}
              header={header}
              footer={footer}
            />
          ) : null}
        </TabPanel>
        <TabPanel id="distribution" className="app-panel">
          {distributionData && metaData ? (
            <Distribution
              data={distributionData}
              metaData={metaData}
              docState={docState}
              setDocState={setDocState}
              header={header}
              footer={footer}
            />
          ) : null}
        </TabPanel>
        {loading ? <div className="app-overlay app-loader" /> : null}
      </div>
    </Tabs>
  );
};
