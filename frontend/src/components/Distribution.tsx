import type {
  DataSets,
  Visualization,
  useDistributionState,
} from "../distributionState";
import type { DistributionData } from "../plotTypes";
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
          <label htmlFor="visualization">Visualization</label>
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
          <label htmlFor="data-set">Data Set</label>
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
