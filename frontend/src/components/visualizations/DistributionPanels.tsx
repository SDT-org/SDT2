import React from "react";
import { TabPanel } from "react-aria-components";
import type { DocState, SetDocState } from "../../appState";
import { type DataSets, useDistributionState } from "../../distributionState";
import type { DistributionData, MetaData } from "../../plotTypes";
import { Select, SelectItem } from "../primitives/Select";
import { Histogram } from "./histogram/Histogram";
import { Violin } from "./violin/Violin";

const VisualizationSwitcher = ({
  activeDataSet,
  setActiveDataSet,
  gcData,
}: {
  activeDataSet: keyof DataSets;
  setActiveDataSet: React.Dispatch<keyof DataSets>;
  gcData: number[];
}) => {
  React.useEffect(() => {
    if (activeDataSet === "gc" && !gcData.some((value) => value > 0)) {
      setActiveDataSet("scores");
    }
  }, [gcData, activeDataSet, setActiveDataSet]);

  return (
    <>
      <div className="group padded">
        <div className="row">
          <div className="field">
            <label className="header" htmlFor="data-set">
              Data Set
            </label>
            <Select
              id="data-set"
              wide
              selectedKey={activeDataSet}
              onSelectionChange={(value) =>
                setActiveDataSet(value as keyof DataSets)
              }
              items={Object.entries({
                scores: "Scores",
                ...(gcData.some((value) => value > 0) && { gc: "GC" }),
                length: "Length",
              }).map(([id, name]) => ({
                id,
                name,
              }))}
            >
              {(item) => (
                <SelectItem textValue={item.name}>{item.name}</SelectItem>
              )}
            </Select>
          </div>
        </div>
      </div>
    </>
  );
};

export const DistributionPanels = ({
  docState,
  setDocState,
  data,
  metaData,
  header,
}: {
  docState: DocState;
  setDocState: SetDocState;
  data: DistributionData | undefined;
  metaData: MetaData;
  header?: React.ReactNode;
  footer?: React.ReactNode;
}) => {
  if (!data) {
    return (
      <div className="app-main centered">
        <h1>No data available</h1>
      </div>
    );
  }
  const dataSets: DataSets = {
    scores: data.raw_mat,
    gc: data.gc_stats,
    length: data.length_stats,
  };

  const {
    distributionState,
    updateDistributionState,
    updateHistogram,
    updateViolin,
  } = useDistributionState(docState, setDocState);

  const sidebarComponent = (
    <>
      <VisualizationSwitcher
        activeDataSet={distributionState.dataSet}
        setActiveDataSet={(value: keyof DataSets) =>
          updateDistributionState({ dataSet: value })
        }
        gcData={dataSets.gc}
      />
    </>
  );

  const commonProps = {
    data,
    metaData,
    dataSets,
    dataSetKey: distributionState.dataSet,
    sidebarComponent,
    header,
  };

  return (
    <>
      <div className="app-sidebar app-sidebar-left">{header}</div>
      <TabPanel id="distribution_histogram" className="app-panel">
        <Histogram
          {...commonProps}
          settings={distributionState.histogram}
          updateSettings={updateHistogram}
        />
      </TabPanel>

      <TabPanel id="distribution_violin" className="app-panel">
        <Violin
          {...commonProps}
          settings={distributionState.violin}
          updateSettings={updateViolin}
        />
      </TabPanel>
    </>
  );
};
