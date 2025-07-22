import React from "react";
import {
  Button,
  Dialog,
  Heading,
  Modal,
  ModalOverlay,
} from "react-aria-components";
import type { AppState, SaveableImageFormat, SetAppState } from "../appState";
import { saveableImageFormats } from "../constants";
import { useDocState } from "../hooks/useDocState";
import { Select, SelectItem } from "./Select";
import { Switch } from "./Switch";

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
  const [showFolder, setShowFolder] = React.useState(appState.openExportFolder);
  const [saveFormat, setSaveFormat] = React.useState<SaveableImageFormat>(
    appState.saveFormat,
  );

  React.useEffect(() => {
    setPrefix(docState.exportPrefix || defaultExportPrefix);
  }, [docState.exportPrefix, defaultExportPrefix]);

  const startExport = React.useCallback(() => {
    setDocState((previous) => ({
      ...previous,
      exportPrefix: prefix,
    }));
    setAppState((previous) => ({
      ...previous,
      exportStatus: "preparing",
      saveFormat,
      showExportModal: false,
      openExportFolder: showFolder,
    }));
  }, [setAppState, setDocState, prefix, showFolder, saveFormat]);

  return (
    <ModalOverlay
      className={`react-aria-ModalOverlay modal-overlay ${appState.exportStatus}`}
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
        <Dialog>
          {({ close }) => (
            <>
              <Heading slot="title">Export images and data</Heading>
              <div className="form export-form">
                <div className={appState.dataExportPath ? "" : groupCss}>
                  <label htmlFor="export-path" className="header">
                    Output folder
                  </label>

                  <div className={`${appState.dataExportPath ? groupCss : ""}`}>
                    {appState.dataExportPath ? (
                      <div className="filename">
                        {appState.dataExportPath.replace(
                          appState.config?.userPath || "",
                          "~",
                        )}
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
                    selectedKey={saveFormat}
                    onSelectionChange={(value) =>
                      setSaveFormat(value as SaveableImageFormat)
                    }
                    items={Object.entries(saveableImageFormats).map(
                      ([id, name]) => ({
                        id,
                        name,
                      }),
                    )}
                  >
                    {(item) => (
                      <SelectItem textValue={item.name}>{item.name}</SelectItem>
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
                      setAppState((previous) => ({
                        ...previous,
                        openExportModal: showFolder,
                      }));
                      setDocState((previous) => ({
                        ...previous,
                        exportPrefix: prefix,
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
              </div>
            </>
          )}
        </Dialog>
      </Modal>
    </ModalOverlay>
  );
};
