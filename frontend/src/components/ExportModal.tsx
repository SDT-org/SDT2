import Plotly from "plotly.js-cartesian-dist-min";
import React from "react";
import {
  Button,
  Dialog,
  Heading,
  Modal,
  ModalOverlay,
} from "react-aria-components";
import useAppState, {
  type DocState,
  type SaveableImageFormat,
  saveableImageFormats,
} from "../appState";
import { assertDefined } from "../helpers";
import { useDocState } from "../hooks/useDocState";
import { useHeatmapRef } from "../hooks/useHeatmapRef";
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
  const { appState, setAppState } = useAppState();
  if (
    !appState.activeDocumentId ||
    !appState.documents.find((doc) => doc.id === appState.activeDocumentId)
  ) {
    return null;
  }
  const [exportState, setExportState] = React.useState<
    "idle" | "exporting" | "success"
  >("idle");
  const [outputCluster, setOutputCluster] = React.useState(false);
  const [thresholds, setThresholds] = React.useState({ one: 79, two: 0 });
  const { docState, updateDocState } = useDocState(
    appState.activeDocumentId,
    appState,
    setAppState,
  );

  const swapDataView = React.useCallback(
    (view: DocState["dataView"]) => {
      updateDocState({ dataView: view });
    },
    [updateDocState],
  );

  const heatmapRef = useHeatmapRef();

  const getImages = React.useCallback(async () => {
    const config = {
      format: appState.saveFormat,
      width: 1000,
      height: 800,
    };
    const previousDataView = docState.dataView;
    const renderTimeout = 400;
    let heatmapImage = "";

    swapDataView("heatmap");
    await new Promise((r) => setTimeout(r, renderTimeout));

    if (config.format === "svg") {
      if (!heatmapRef.current) {
        throw new Error("Expected heatmapRef to have a current value");
      }
      const encoded64Svg = encodeURIComponent(heatmapRef.current.outerHTML);
      heatmapImage = `data:image/svg+xml;base64,${encoded64Svg}`;
    } else {
      heatmapImage = await new Promise((resolve) => {
        if (!heatmapRef.current) {
          throw new Error("Expected heatmapRef to have a current value");
        }

        (heatmapRef.current as HTMLCanvasElement).toBlob(async (blob) => {
          if (blob) {
            const arrayBuffer = await blob.arrayBuffer();
            const binary = Array.from(new Uint8Array(arrayBuffer))
              .map((byte) => String.fromCharCode(byte))
              .join("");
            resolve(`data:image/${config.format};base64,${btoa(binary)}`);
          } else {
            resolve("");
          }
        }, `image/${config.format}`);
      });
    }

    await new Promise((r) => setTimeout(r, renderTimeout));

    // Plotly exports

    swapDataView("distribution_histogram");

    await new Promise((r) => setTimeout(r, renderTimeout));

    let element = getPlotlyElement();

    element = getPlotlyElement();
    await Plotly.toImage(element, config);
    const histogramImage = await Plotly.toImage(element, config);

    swapDataView("distribution_violin");

    await new Promise((r) => setTimeout(r, renderTimeout));
    element = getPlotlyElement();
    await Plotly.toImage(element, config);
    const violinImage = await Plotly.toImage(element, config);

    // swapDataView("distribution_raincloud");

    // await new Promise((r) => setTimeout(r, renderTimeout));
    // element = getPlotlyElement();
    // await Plotly.toImage(element, config);
    // const raincloudImage = await Plotly.toImage(element, config);

    swapDataView(previousDataView);
    return { heatmapImage, histogramImage, violinImage };
  }, [heatmapRef, appState, docState.dataView, swapDataView]);

  const doExport = React.useCallback(() => {
    setExportState("exporting");
    getImages().then((images) => {
      window.pywebview.api
        .export_data({
          doc_id: appState.activeDocumentId,
          export_path: appState.dataExportPath,
          output_cluster: outputCluster,
          cluster_threshold_one: thresholds.one,
          cluster_threshold_two: thresholds.two,
          heatmap_image_data: images.heatmapImage,
          histogram_image_data: images.histogramImage,
          violin_image_data: images.violinImage,
          image_format: appState.saveFormat,
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
        showExportModal: false,
      }));
    };

    const id = setTimeout(handler, 2000);

    return () => {
      clearTimeout(id);
    };
  }, [exportState, setAppState]);

  const groupCss = "group col-2 onefr-auto";

  return (
    <ModalOverlay
      className={`react-aria-ModalOverlay export-modal-overlay ${exportState}`}
      isOpen={appState.showExportModal}
      onOpenChange={() =>
        setAppState((previous) => {
          return {
            ...previous,
            showExportModal: false,
          };
        })
      }
    >
      <Modal
        isOpen={appState.showExportModal}
        onOpenChange={() =>
          setAppState((previous) => {
            return {
              ...previous,
              showExportModal: false,
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
                    <div className={appState.dataExportPath ? "" : groupCss}>
                      <label htmlFor="export-path" className="header">
                        Output Folder
                      </label>

                      <div
                        className={`${appState.dataExportPath ? groupCss : ""}`}
                      >
                        {appState.dataExportPath ? (
                          <div className="filename">
                            {appState.dataExportPath}
                          </div>
                        ) : null}
                        <Button
                          type="button"
                          onPress={() => {
                            setExportState("idle");
                            window.pywebview.api
                              .select_path_dialog(appState.dataExportPath)
                              .then((result) => {
                                if (!result) {
                                  return;
                                }
                                setAppState((prev) => ({
                                  ...prev,
                                  dataExportPath: result,
                                }));
                              });
                          }}
                        >
                          {appState.dataExportPath ? "Set" : "Set folder"}
                          ...
                        </Button>
                      </div>
                    </div>
                    <div className="group col-2 onefr-auto force-flat">
                      <label htmlFor="save-format" className="header">
                        Image Format
                      </label>
                      <Select
                        data-flat
                        selectedKey={appState.saveFormat}
                        onSelectionChange={(value) =>
                          setAppState((previous) => ({
                            ...previous,
                            saveFormat: value as SaveableImageFormat,
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
                        isDisabled={!appState.dataExportPath}
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
