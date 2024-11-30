import type React from "react";
import type {
  DataSets,
  Visualization,
  useDistributionState,
} from "../distributionState";
import type { DistributionData } from "../plotTypes";
import { Histogram } from "./Histogram";
import { Raincloud } from "./Raincloud";
import { Violin } from "./Violin";

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
  state,
  setState,
  updateHistogram,
  updateRaincloud,
  updateViolin,
}: {
  data: DistributionData | undefined;
  footer?: React.ReactNode;
} & ReturnType<typeof useDistributionState>) => {
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

  const sidebarComponent = (
    <VisualizationSwitcher
      activeDataSet={state.dataSet}
      setActiveDataSet={(value: keyof DataSets) =>
        setState((prev) => ({ ...prev, dataSet: value }))
      }
      visualization={state.visualization}
      setVisualization={(value: Visualization) =>
        setState((prev) => ({ ...prev, visualization: value }))
      }
    />
  );

  const commonProps = {
    data,
    dataSets,
    dataSetKey: state.dataSet,
    sidebarComponent,
  };

  const components = {
    histogram: (
      <Histogram
        {...commonProps}
        settings={state.histogram}
        updateSettings={updateHistogram}
      />
    ),
    violin: (
      <Violin
        {...commonProps}
        settings={state.violin}
        updateSettings={updateViolin}
      />
    ),
    raincloud: (
      <Raincloud
        {...commonProps}
        settings={state.raincloud}
        updateSettings={updateRaincloud}
      />
    ),
  };

  return components[state.visualization];
};
