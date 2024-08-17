import React from "react";
import { Modal, Button, Dialog, DialogTrigger } from "react-aria-components";
import Plotly from "plotly.js-dist-min";
import useAppState, { AppState, SetAppState } from "../appState";
import { NumberInput } from "./NumberInput";

export const ExportData = () => {
  const [exportState, setExportState] = React.useState<
    "idle" | "exporting" | "success"
  >("idle");
  const { appState, setAppState } = useAppState();
  const [outputCluster, setOutputCluster] = React.useState(false);
  const [thresholds, setThresholds] = React.useState({ one: 79, two: 0 });

  const swapDataView = (view: AppState["client"]["dataView"]) =>
    setAppState((previous) => {
      return {
        ...previous,
        client: {
          ...previous.client,
          dataView: view,
        },
      };
    });

  const getImages = async () => {
    const previousDataView = appState.client.dataView;
    swapDataView("heatmap");

    let element = document.getElementsByClassName("js-plotly-plot")[0];
    const config = {
      format: "svg",
      width: 1000,
      height: 800,
    };

    await Plotly.toImage(element, config);
    const heatmapImage = await Plotly.toImage(element, config);

    swapDataView("plot");
    await new Promise((r) => setTimeout(r, 100));
    element = document.getElementsByClassName("js-plotly-plot")[0];

    await Plotly.toImage(element, config);
    const distributionImage = await Plotly.toImage(element, config);
    swapDataView(previousDataView);
    return [heatmapImage, distributionImage];
  };

  React.useEffect(() => {
    if (exportState !== "exporting") {
      return;
    }

    getImages().then((images) => {
      window.pywebview.api
        .export_data({
          output_cluster: outputCluster,
          cluster_threshold_one: thresholds.one,
          cluster_threshold_two: thresholds.two,
          heatmap_image_data: images[0],
          distribution_image_data: images[1],
        })
        .then((result) =>
          result ? setExportState("success") : setExportState("idle"),
        );
    });
  }, [exportState]);

  return (
    <DialogTrigger>
      <Button>Export Data</Button>
      <Modal
        isDismissable={exportState != "exporting"}
        className={`react-aria-Modal export-modal ${exportState}`}
      >
        <Dialog>
          {({ close }) => (
            <>
              <form className="form export-form">
                <h1>Export Data</h1>
                <div className="group">
                  <div className="field">
                    <label className="header">Output Folder</label>
                    <input type="text" readOnly value={appState.export_path} />
                    <button
                      type="button"
                      disabled={exportState === "exporting"}
                      onClick={() => {
                        setExportState("idle");
                        window.pywebview.api.select_export_path();
                      }}
                    >
                      Select Folder...
                    </button>
                  </div>
                </div>
                <div className="group">
                  <div className="field">
                    <label className="header">
                      <input
                        type="checkbox"
                        checked={outputCluster}
                        disabled={exportState === "exporting"}
                        onChange={() => {
                          setExportState("idle");
                          setOutputCluster(!outputCluster);
                        }}
                      />
                      <span>Cluster by Percent Identity</span>
                    </label>
                  </div>
                  <div className="col-2">
                    <NumberInput
                      label="Threshold 1"
                      field="one"
                      value={thresholds.one}
                      isDisabled={!outputCluster || exportState === "exporting"}
                      updateValue={(newValue) =>
                        setThresholds({ ...thresholds, ...newValue })
                      }
                      type="int"
                      min={0}
                      max={100}
                      step={1}
                    />
                    <NumberInput
                      label="Threshold 2"
                      field="two"
                      value={thresholds.two}
                      isDisabled={!outputCluster || exportState === "exporting"}
                      updateValue={(newValue) =>
                        setThresholds({ ...thresholds, ...newValue })
                      }
                      type="int"
                      min={0}
                      max={100}
                      step={1}
                    />
                  </div>
                </div>
                {exportState === "success" ? (
                  <div className="success">✓ Export complete</div>
                ) : null}
                <div className="actions">
                  <Button
                    onPress={() => {
                      setExportState("idle");
                      close();
                    }}
                  >
                    Cancel
                  </Button>
                  <Button
                    isDisabled={
                      exportState === "exporting" || !appState.export_path
                    }
                    onPress={() => setExportState("exporting")}
                  >
                    {exportState === "exporting" ? "Exporting..." : "Export"}
                  </Button>
                </div>
              </form>
            </>
          )}
        </Dialog>
      </Modal>
    </DialogTrigger>
  );
};
