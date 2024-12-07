import Plotly from "plotly.js-dist-min";
import React from "react";
import {
  Button,
  Dialog,
  Heading,
  Modal,
  ModalOverlay,
} from "react-aria-components";
import useAppState, {
  saveableImageFormats,
  type AppState,
  type SaveableImageFormat,
} from "../appState";
import { useDistributionState } from "../distributionState";
import { assertDefined } from "../helpers";
import { Select, SelectItem } from "./Select";
import { Slider } from "./Slider";
import { Switch } from "./Switch";

const getPlotlyElement = () =>
  assertDefined(
    (
      document.getElementsByClassName(
        "js-plotly-plot",
      ) as HTMLCollectionOf<HTMLElement>
    )[0],
  );

const SuccessMessage = () => (
  <div className="app-message">
    <svg
      className="app-success"
      xmlns="http://www.w3.org/2000/svg"
      viewBox="0 0 24 24"
      aria-hidden="true"
    >
      <g
        style={{
          fill: "none",
          stroke: "currentcolor",
          strokeWidth: 2,
          strokeLinecap: "round",
          strokeLinejoin: "round",
          strokeMiterlimit: 10,
        }}
      >
        <path d="m8.5 13 3 2 4-6" />
        <path d="M18 22H6a4 4 0 0 1-4-4V6a4 4 0 0 1 4-4h12a4 4 0 0 1 4 4v12a4 4 0 0 1-4 4z" />
      </g>
    </svg>
    <Heading slot="title">Export complete</Heading>
  </div>
);

export const ExportModal = () => {
  const [exportState, setExportState] = React.useState<
    "idle" | "exporting" | "success"
  >("idle");
  const { appState, setAppState } = useAppState();
  const [outputCluster, setOutputCluster] = React.useState(false);
  const [thresholds, setThresholds] = React.useState({ one: 79, two: 0 });
  const { updateDistributionState } = useDistributionState(
    appState,
    setAppState,
  );

  const swapDataView = React.useCallback(
    (view: AppState["client"]["dataView"]) =>
      setAppState((previous) => {
        return {
          ...previous,
          client: {
            ...previous.client,
            dataView: view,
          },
        };
      }),
    [setAppState],
  );

  const getImages = React.useCallback(async () => {
    const previousDataView = appState.client.dataView;
    const renderTimeout = 400;
    await new Promise((r) => setTimeout(r, 200));
    swapDataView("heatmap");
    await new Promise((r) => setTimeout(r, renderTimeout));

    let element = getPlotlyElement();

    const config = {
      format: appState.client.saveFormat,
      width: 1000,
      height: 800,
    };

    await Plotly.toImage(element, config);
    const heatmapImage = await Plotly.toImage(element, config);

    swapDataView("distribution");

    updateDistributionState({ visualization: "histogram" });
    await new Promise((r) => setTimeout(r, renderTimeout));
    element = getPlotlyElement();
    await Plotly.toImage(element, config);
    const histogramImage = await Plotly.toImage(element, config);

    updateDistributionState({ visualization: "violin" });
    await new Promise((r) => setTimeout(r, renderTimeout));
    element = getPlotlyElement();
    await Plotly.toImage(element, config);
    const violinImage = await Plotly.toImage(element, config);

    updateDistributionState({ visualization: "raincloud" });
    await new Promise((r) => setTimeout(r, renderTimeout));
    element = getPlotlyElement();
    await Plotly.toImage(element, config);
    const raincloudImage = await Plotly.toImage(element, config);

    swapDataView(previousDataView);
    return { heatmapImage, histogramImage, violinImage, raincloudImage };
  }, [appState, swapDataView, updateDistributionState]);

  const doExport = React.useCallback(() => {
    setExportState("exporting");
    getImages().then((images) => {
      window.pywebview.api
        .export_data({
          export_path: appState.client.dataExportPath,
          output_cluster: outputCluster,
          cluster_threshold_one: thresholds.one,
          cluster_threshold_two: thresholds.two,
          heatmap_image_data: images.heatmapImage,
          histogram_image_data: images.histogramImage,
          violin_image_data: images.violinImage,
          raincloud_image_data: images.raincloudImage,
          image_format: appState.client.saveFormat,
        })
        .then((result) =>
          result ? setExportState("success") : setExportState("idle"),
        );
    });
  }, [getImages, appState, outputCluster, thresholds]);

  React.useEffect(() => {
    if (exportState !== "success") {
      return;
    }

    const handler = () => {
      setExportState("idle");
      setAppState((prev) => ({
        ...prev,
        client: { ...prev.client, showExportModal: false },
      }));
    };

    const id = setTimeout(handler, 2000);

    return () => {
      clearTimeout(id);
    };
  }, [exportState, setAppState]);

  return (
    <ModalOverlay
      className={`react-aria-ModalOverlay export-modal-overlay ${exportState}`}
      isOpen={appState.client.showExportModal}
      onOpenChange={() =>
        setAppState((previous) => {
          return {
            ...previous,
            client: { ...previous.client, showExportModal: false },
          };
        })
      }
    >
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
        {exportState === "exporting" ? (
          <div className="app-overlay app-loader" />
        ) : (
          <Dialog>
            {({ close }) =>
              exportState === "success" ? (
                <SuccessMessage />
              ) : (
                <>
                  <Heading slot="title">Export images and data</Heading>
                  <form className="form export-form">
                    <label htmlFor="export-path" className="header">
                      Output Folder
                    </label>

                    <div className="group  col-2 onefr-auto">
                      <div className="filename">
                        {appState.client.dataExportPath}
                      </div>
                      <Button
                        type="button"
                        onPress={() => {
                          setExportState("idle");
                          window.pywebview.api
                            .select_path_dialog(appState.client.dataExportPath)
                            .then((result) => {
                              if (!result) {
                                return;
                              }
                              setAppState((prev) => ({
                                ...prev,
                                client: {
                                  ...prev.client,
                                  dataExportPath: result,
                                },
                              }));
                            });
                        }}
                      >
                        Set...
                      </Button>
                    </div>
                    <div className="group col-2 onefr-auto force-flat">
                      <label htmlFor="save-format" className="header">
                        Image Format
                      </label>
                      <Select
                        data-flat
                        selectedKey={appState.client.saveFormat}
                        onSelectionChange={(value) =>
                          setAppState((previous) => ({
                            ...previous,
                            client: {
                              ...previous.client,
                              saveFormat: value as SaveableImageFormat,
                            },
                          }))
                        }
                        items={Object.entries(saveableImageFormats).map(
                          ([id, name]) => ({
                            id,
                            name,
                          }),
                        )}
                      >
                        {(item) => (
                          <SelectItem textValue={item.name}>
                            {item.name}
                          </SelectItem>
                        )}
                      </Select>
                    </div>
                    <div className="group">
                      <Switch
                        data-split
                        isSelected={outputCluster}
                        onChange={() => {
                          setExportState("idle");
                          setOutputCluster(!outputCluster);
                        }}
                      >
                        Cluster by Percent Identity
                      </Switch>
                      <div
                        className="field col-2 subfield"
                        data-hidden={!outputCluster}
                        aria-hidden={!outputCluster}
                      >
                        <Slider
                          label="Threshold 1"
                          value={thresholds.one}
                          isDisabled={!outputCluster}
                          onChange={(newValue) =>
                            setThresholds((prev) => ({
                              ...prev,
                              one: newValue,
                            }))
                          }
                          minValue={0}
                          maxValue={100}
                          step={1}
                          includeField
                        />
                        <Slider
                          label="Threshold 2"
                          value={thresholds.two}
                          isDisabled={!outputCluster}
                          onChange={(newValue) =>
                            setThresholds((prev) => ({
                              ...prev,
                              two: newValue,
                            }))
                          }
                          minValue={0}
                          maxValue={100}
                          step={1}
                          includeField
                        />
                      </div>
                    </div>
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
                        data-primary
                        isDisabled={!appState.client.dataExportPath}
                        onPress={doExport}
                      >
                        Export
                      </Button>
                    </div>
                  </form>
                </>
              )
            }
          </Dialog>
        )}
      </Modal>
    </ModalOverlay>
  );
};
