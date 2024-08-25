import React from "react";
import { Modal, Button, Dialog, DialogTrigger } from "react-aria-components";
import Plotly from "plotly.js-dist-min";
import useAppState, { AppState, SaveableImageFormat } from "../appState";
import { NumberInput } from "./NumberInput";
import { assertDefined } from "../helpers";

export const ExportModal = () => {
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

  const getPlotlyElement = () =>
    assertDefined(
      (
        document.getElementsByClassName(
          "js-plotly-plot",
        ) as HTMLCollectionOf<HTMLElement>
      )[0],
    );

  const getImages = async () => {
    const previousDataView = appState.client.dataView;
    swapDataView("heatmap");

    let element = getPlotlyElement();

    const config = {
      format: appState.client.saveFormat,
      width: 1000,
      height: 800,
    };

    await Plotly.toImage(element, config);
    const heatmapImage = await Plotly.toImage(element, config);

    swapDataView("plot");
    await new Promise((r) => setTimeout(r, 100));
    element = getPlotlyElement();

    await Plotly.toImage(element, config);
    const distributionImage = await Plotly.toImage(element, config);

    swapDataView(previousDataView);
    return { heatmapImage, distributionImage };
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
          heatmap_image_data: images.heatmapImage,
          distribution_image_data: images.distributionImage,
          image_format: appState.client.saveFormat,
        })
        .then((result) =>
          result ? setExportState("success") : setExportState("idle"),
        );
    });
  }, [exportState]);

  return (
    <Modal
      isOpen={appState.client.showExportModal}
      onOpenChange={() =>
        setAppState((previous) => {
          return {
            ...previous,
            client: { ...previous.client, showExportModal: false },
          };
        })
      }
      isDismissable={exportState !== "exporting"}
      className={`react-aria-Modal export-modal ${exportState}`}
    >
      <Dialog>
        {({ close }) => (
          <>
            <form className="form export-form">
              <h1>Export Images & Data</h1>
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
              <div className="field">
                <label htmlFor="save-format" className="header">
                  Image Format
                </label>
                <select
                  onChange={(e) =>
                    setAppState((previous) => ({
                      ...previous,
                      client: {
                        ...previous.client,
                        saveFormat: e.target.value as SaveableImageFormat,
                      },
                    }))
                  }
                  value={appState.client.saveFormat}
                >
                  <option value="png">PNG</option>
                  <option value="jpeg">JPEG</option>
                  <option value="svg">SVG</option>
                </select>
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
  );
};
