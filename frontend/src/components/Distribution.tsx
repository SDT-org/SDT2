import React from "react";
import { Histogram } from "./Histogram";
import { Violin } from "./Violin";
import { Raincloud } from "./Raincloud";
import { DistributionData } from "../plotTypes";

export type Visualization = "histogram" | "violin" | "raincloud";
type DataSet = number[];
export type DataSets = {
  scores: DataSet;
  gc: DataSet;
  length: DataSet;
};

const VisualizationSwitcher = ({
  activeDataSet,
  setActiveDataSet,
  visualization,
  setVisualization,
}: {
  activeDataSet: keyof DataSets;
  setActiveDataSet: React.Dispatch<keyof DataSets>;
  visualization: Visualization;
  setVisualization: React.Dispatch<Visualization>;
}) => (
  <>
    <div className="group">
      <div className="row">
        <div className="field">
          <label htmlFor="visualization">Visualization</label>
          <select
            id="visualization"
            value={visualization}
            onChange={(e) => setVisualization(e.target.value as Visualization)}
          >
            <option value="histogram">Histogram</option>
            <option value="violin">Violin</option>
            <option value="raincloud">Raincloud</option>
          </select>
        </div>
      </div>
      <div className="row">
        <div className="field">
          <label htmlFor="data-source">Data Set</label>
          <select
            id="data-set"
            value={activeDataSet}
            onChange={(e) => setActiveDataSet(e.target.value as keyof DataSets)}
          >
            <option value="scores">Scores</option>
            <option value="gc">GC</option>
            <option value="length">Length</option>
          </select>
        </div>
      </div>
    </div>
  </>
);

export const Distribution = ({
  data,
  tickText,
}: {
  data: DistributionData | undefined;
  tickText: string[];
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

  const [visualization, setVisualization] =
    React.useState<Visualization>("histogram");
  const [dataSetKey, setDataSetKey] = React.useState<keyof DataSets>("scores");

  const sidebarComponent = (
    <VisualizationSwitcher
      activeDataSet={dataSetKey}
      setActiveDataSet={setDataSetKey}
      visualization={visualization}
      setVisualization={setVisualization}
    />
  );

  return (
    <>
      {visualization === "histogram" && (
        <Histogram
          data={data}
          dataSets={dataSets}
          dataSetKey={dataSetKey}
          sidebarComponent={sidebarComponent}
        />
      )}
      {visualization === "violin" && (
        <Violin
          data={data}
          dataSets={dataSets}
          dataSetKey={dataSetKey}
          tickText={tickText}
          sidebarComponent={sidebarComponent}
        />
      )}
      {visualization === "raincloud" && (
        <Raincloud
          data={data}
          dataSets={dataSets}
          dataSetKey={dataSetKey}
          sidebarComponent={sidebarComponent}
        />
      )}
    </>
  );
};
