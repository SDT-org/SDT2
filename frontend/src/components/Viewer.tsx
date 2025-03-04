import React from "react";
import { type Key, Tab, TabList, TabPanel, Tabs } from "react-aria-components";
import type { DocState, SetDocState, UpdateDocState } from "../appState";
import { useGetData } from "../hooks/useGetData";
import { Clustermap } from "./Clustermap";
import { DistributionPanels } from "./DistributionPanels";
import { Heatmap } from "./Heatmap";

export const Viewer = ({
  docState,
  setDocState,
  updateDocState,
  leftSidebarCollapsed,
}: {
  docState: DocState;
  setDocState: SetDocState;
  updateDocState: UpdateDocState;
  leftSidebarCollapsed: boolean;
  tabView: "tabs" | "select";
}) => {
  const { loading, tickText, heatmapData, distributionData, metaData } =
    useGetData(docState, setDocState);

  const setDataView = React.useCallback(
    (newValue: Key) => {
      updateDocState(
        {
          dataView: newValue as DocState["dataView"],
        },
        false,
      );
    },
    [updateDocState],
  );

  return (
    <TabPanel id={docState.id} key={docState.id}>
      <Tabs
        selectedKey={docState.dataView}
        onSelectionChange={setDataView}
        className="react-aria-Tabs app-panels with-sidebar"
        orientation="vertical"
        data-left-sidebar-collapsed={leftSidebarCollapsed}
      >
        <div className="app-sidebar app-sidebar-left heatmap-sidebar">
          <TabList data-collapsed={leftSidebarCollapsed}>
            <Tab
              id="pairwise"
              isDisabled={true}
              className={"react-aria-Tab header"}
              data-hidden={leftSidebarCollapsed}
            >
              Pairwise
            </Tab>
            <Tab id="heatmap">
              <div
                style={{ display: "flex", alignItems: "center", gap: "0.8rem" }}
              >
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
                <span>Heatmap</span>
              </div>
            </Tab>
            <Tab id="cluster">
              <div
                style={{ display: "flex", alignItems: "center", gap: "0.8rem" }}
              >
                <svg
                  height={12}
                  xmlns="http://www.w3.org/2000/svg"
                  viewBox="0 0 30 30"
                  aria-hidden={true}
                  color="currentcolor"
                >
                  <circle cx="6" cy="6" r="5" />
                  <circle cx="25" cy="6" r="5" />
                  <circle cx="6" cy="24" r="5" />
                  <circle cx="25" cy="24" r="5" />
                  <circle cx="16" cy="15" r="4" />
                </svg>
                <span>Clustermap</span>
              </div>
            </Tab>
            <Tab
              id="distribution"
              isDisabled={true}
              className={"react-aria-Tab header"}
            >
              {leftSidebarCollapsed ? <hr /> : "Distribution"}
            </Tab>
            <Tab id="distribution_histogram">
              <div
                style={{
                  display: "flex",
                  alignItems: "center",
                  gap: "0.8rem",
                }}
              >
                <svg
                  height={12}
                  xmlns="http://www.w3.org/2000/svg"
                  viewBox="0 0 30 30"
                  aria-hidden={true}
                  color="currentcolor"
                >
                  <rect x="0" y="18" width={6} height={10} fill="black" />
                  <rect x="10" y="0" width={6} height={28} fill="black" />
                  <rect x="20" y="10" width={6} height={18} fill="black" />
                </svg>
                <span>Histogram</span>
              </div>
            </Tab>
            <Tab id="distribution_violin">
              <div
                style={{
                  display: "flex",
                  alignItems: "center",
                  gap: "0.8rem",
                }}
              >
                <svg
                  height={12}
                  xmlns="http://www.w3.org/2000/svg"
                  viewBox="0 0 30 30"
                  aria-hidden={true}
                  color="currentcolor"
                >
                  <path
                    d="M15 1.5c1 0 5 5 5 9s-4 9-5 9-5-5-5-9 4-9 5-9z"
                    fill="black"
                  />
                  <rect x="12" y="12" width="6" height="9" fill="black" />
                  <rect x="7" y="11.25" width="16" height="7.5" fill="black" />
                  <circle cx="15" cy="23.25" r="1.75" fill="black" />
                  <circle cx="15" cy="26.75" r="1.75" fill="black" />
                </svg>
                <span>Violin</span>
              </div>
            </Tab>
          </TabList>
          <div className="app-sidebar-body" />
          <div className="app-sidebar-footer" />
        </div>
        <TabPanel id="heatmap" className="app-panel">
          {!loading && heatmapData && metaData ? (
            <Heatmap
              data={heatmapData}
              metaData={metaData}
              tickText={tickText}
              docState={docState}
              setDocState={setDocState}
              updateDocState={updateDocState}
              leftSidebarCollapsed={leftSidebarCollapsed}
            />
          ) : null}
        </TabPanel>
        <TabPanel id="cluster" className="app-panel">
          {!loading && heatmapData ? (
            <Clustermap
              data={heatmapData}
              docState={docState}
              setDocState={setDocState}
              tickText={tickText}
              leftSidebarCollapsed={leftSidebarCollapsed}
            />
          ) : null}
        </TabPanel>
        {distributionData && metaData ? (
          <DistributionPanels
            data={distributionData}
            metaData={metaData}
            docState={docState}
            setDocState={setDocState}
          />
        ) : null}
        {loading ? <div className="app-overlay app-loader" /> : null}
      </Tabs>
    </TabPanel>
  );
};
