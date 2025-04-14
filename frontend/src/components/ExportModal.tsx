import React from "react";
import {
  Button,
  Dialog,
  Heading,
  Modal,
  ModalOverlay,
} from "react-aria-components";
import {
  type AppState,
  type SaveableImageFormat,
  type SetAppState,
  saveableImageFormats,
} from "../appState";
import { useDocState } from "../hooks/useDocState";
import { Select, SelectItem } from "./Select";
import { Switch } from "./Switch";

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

export const ExportModal = ({
  appState,
  setAppState,
}: {
  appState: AppState;
  setAppState: SetAppState;
}) => {
  if (
    !appState.activeDocumentId ||
    !appState.documents.find((doc) => doc.id === appState.activeDocumentId)
  ) {
    return null;
  }
  const { docState, setDocState } = useDocState(
    appState.activeDocumentId,
    appState,
    setAppState,
  );

  const setExportState = React.useCallback(
    (status: AppState["exportStatus"]) =>
      setAppState((prev) => ({ ...prev, exportStatus: status })),
    [setAppState],
  );

  const groupCss = "group col-2 onefr-auto";
  const defaultExportPrefix = docState.basename
    .split(".")
    .slice(0, -1)
    .join("");

  const [prefix, setPrefix] = React.useState<string>(
    docState.exportPrefix || defaultExportPrefix,
  );
  const [showFolder, setShowFolder] = React.useState(false);

  React.useEffect(() => {
    setPrefix(docState.exportPrefix || defaultExportPrefix);
  }, [docState.exportPrefix, defaultExportPrefix]);

  const startExport = React.useCallback(() => {
    setExportState("exporting");
    setDocState((previous) => ({
      ...previous,
      exportPrefix: prefix,
      openExportFolder: showFolder,
    }));
  }, [setDocState, prefix, showFolder, setExportState]);

  return (
    <ModalOverlay
      className={`react-aria-ModalOverlay export-modal-overlay ${appState.exportStatus}`}
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
        isDismissable={appState.exportStatus !== "exporting"}
        className={`react-aria-Modal export-modal ${appState.exportStatus}`}
      >
        {appState.exportStatus === "exporting" ? (
          <div className="app-overlay app-loader" />
        ) : (
          <Dialog>
            {({ close }) =>
              appState.exportStatus === "success" ? (
                <SuccessMessage />
              ) : (
                <>
                  <Heading slot="title">Export images and data</Heading>
                  <form className="form export-form">
                    <div className={appState.dataExportPath ? "" : groupCss}>
                      <label htmlFor="export-path" className="header">
                        Output folder
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
                    <div className="group col-2 onefr-auto">
                      <label htmlFor="prefix" className="header">
                        File prefix
                      </label>
                      <input
                        type="text"
                        id="prefix"
                        value={prefix}
                        onChange={(e) => {
                          setPrefix(e.target.value);
                        }}
                        placeholder={defaultExportPrefix}
                      />
                    </div>
                    <div className="group col-2 onefr-auto force-flat">
                      <label htmlFor="save-format" className="header">
                        Image format
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
                    <div className="group col-2 onefr-auto">
                      <label htmlFor="open-export-folder" className="header">
                        Open folder after export
                      </label>
                      <Switch
                        id="open-export-folder"
                        isSelected={showFolder}
                        onChange={setShowFolder}
                      />
                    </div>
                    <div className="actions">
                      <Button
                        onPress={() => {
                          setDocState((previous) => ({
                            ...previous,
                            exportPrefix: prefix,
                            openExportFolder: showFolder,
                          }));
                          setExportState("idle");
                          close();
                        }}
                      >
                        Cancel
                      </Button>
                      <Button
                        data-primary
                        isDisabled={
                          !appState.dataExportPath ||
                          appState.exportStatus === "exporting"
                        }
                        onPress={startExport}
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
