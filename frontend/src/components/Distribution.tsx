import useAppState from "../appState";
import {
  type DataSets,
  type Visualization,
  useDistributionState,
} from "../distributionState";
import type { DistributionData, MetaData } from "../plotTypes";
import { Histogram } from "./Histogram";
import { Raincloud } from "./Raincloud";
import { Select, SelectItem } from "./Select";
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
    <div className="group padded">
      <div className="row">
        <div className="field">
          <label className="header" htmlFor="visualization">
            Visualization
          </label>
          <Select
            id="visualization"
            wide
            selectedKey={visualization}
            onSelectionChange={(value) =>
              setVisualization(value as Visualization)
            }
            items={Object.entries({
              histogram: "Histogram",
              violin: "Violin",
              raincloud: "Raincloud",
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
              gc: "GC",
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

export const Distribution = ({
  data,
  metaData,
}: {
  data: DistributionData | undefined;
  metaData: MetaData;
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

  const { appState, setAppState } = useAppState();
  const {
    distributionState,
    updateDistributionState,
    updateHistogram,
    updateRaincloud,
    updateViolin,
  } = useDistributionState(appState, setAppState);

  const sidebarComponent = (
    <VisualizationSwitcher
      activeDataSet={distributionState.dataSet}
      setActiveDataSet={(value: keyof DataSets) =>
        updateDistributionState({ dataSet: value })
      }
      visualization={distributionState.visualization}
      setVisualization={(value: Visualization) =>
        updateDistributionState({ visualization: value })
      }
    />
  );

  const commonProps = {
    data,
    metaData,
    dataSets,
    dataSetKey: distributionState.dataSet,
    sidebarComponent,
  };

  const components = {
    histogram: (
      <Histogram
        {...commonProps}
        settings={distributionState.histogram}
        updateSettings={updateHistogram}
      />
    ),
    violin: (
      <Violin
        {...commonProps}
        settings={distributionState.violin}
        updateSettings={updateViolin}
      />
    ),
    raincloud: (
      <Raincloud
        {...commonProps}
        settings={distributionState.raincloud}
        updateSettings={updateRaincloud}
      />
    ),
  };

  return components[distributionState.visualization];
};
