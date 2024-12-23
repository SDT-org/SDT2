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
}: {
  docState: DocState;
  setDocState: SetDocState;
  updateDocState: UpdateDocState;
  mainMenu: React.ReactNode;
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

  return (
    <Tabs selectedKey={docState.dataView} onSelectionChange={setDataView}>
      <div className="app-wrapper with-header">
        <div className="app-header">
          <div className="left">
            {mainMenu}
            <SelectDocumentMenu />
            <div className="run-info">
              {docState.sequences_count > 0 ? (
                <>
                  {docState.sequences_count} Sequence
                  {docState.sequences_count === 1 ? "" : "s"}
                </>
              ) : null}
              <span className="filename">{docState.basename}</span>
            </div>
          </div>

          <TabList>
            <Tab id="heatmap">Heatmap</Tab>
            <Tab id="distribution">Distribution</Tab>
          </TabList>

          <div className="right">
            <Button
              onPress={() =>
                setAppState((prev) => ({ ...prev, showExportModal: true }))
              }
            >
              Export
            </Button>
          </div>
        </div>

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
            />
          ) : null}
        </TabPanel>
        {loading ? <div className="app-overlay app-loader" /> : null}
      </div>
    </Tabs>
  );
};
