import React from "react";
import {
  Button,
  Input,
  Label,
  Radio,
  RadioGroup,
  TabPanel,
  Text,
} from "react-aria-components";
import { TbFile } from "react-icons/tb";
import useAppState, {
  type DocState,
  type SetDocState,
  type UpdateDocState,
} from "../appState";
import { reorderMethods } from "../constants";
import { splitFilePath } from "../helpers";
import useOpenFileDialog from "../hooks/useOpenFileDialog";
import { useRecentFiles } from "../hooks/useRecentFiles";
import { useStartRun } from "../hooks/useStartRun";
import messages from "../messages";
import { RunnerPerformance } from "./RunnerPerformance";
import { Select, SelectItem } from "./Select";
import { Switch } from "./Switch";

export type RunSettings = Pick<DocState, "compute_cores" | "analysisMethod"> & {
  doc_id: string;
  cluster_method: DocState["cluster_method"] | "None";
};

const RunnerSettings = ({
  docState,
  setDocState,
  updateDocState,
}: {
  docState: DocState;
  setDocState: SetDocState;
  updateDocState: UpdateDocState;
}) => {
  const { appState, setAppState } = useAppState();
  const startRun = useStartRun(docState);
  const [startingRun, setStartingRun] = React.useState(false);
  const initialized = React.useRef(false);
  const openFileDialog = useOpenFileDialog(appState, setAppState);
  const openRecentFile = useRecentFiles(appState, setAppState);
  const fileName =
    docState.filename?.length && docState.filename
      ? docState.filename.split("/").pop()
      : undefined;
  const isFastaType = docState.filetype === "text/fasta";

  React.useEffect(() => {
    const handleEnter = (event: KeyboardEvent) => {
      if (
        event.key === "Enter" &&
        (event.ctrlKey || event.metaKey || event.altKey) &&
        Boolean(docState.filename) &&
        !docState.validation_error_id
      ) {
        event.preventDefault();
        startRun();
      }
    };
    document.addEventListener("keydown", handleEnter);

    return () => {
      document.removeEventListener("keydown", handleEnter);
    };
  }, [docState.filename, docState.validation_error_id, startRun]);

  return (
    <div className="form-wrapper runner-wrapper">
      <div className="form runner-form">
        {isFastaType && !docState.validation_error_id ? (
          <>
            <div className="field">
              <label htmlFor="data-file" className="header">
                Data file{" "}
                <span className="text-normal-weight">FASTA or SDT matrix</span>
              </label>
              <div className="setting input-with-button">
                <Input
                  id="data-file"
                  type="text"
                  readOnly
                  value={docState.validation_error_id ? "" : (fileName ?? "")}
                />
                <Button
                  type="button"
                  onPress={() => openFileDialog(docState.id)}
                >
                  Select file&#8230;
                </Button>
              </div>
            </div>
            <div className="field">
              <RadioGroup
                data-card
                onChange={(value) => {
                  value &&
                    updateDocState({
                      analysisMethod: value as DocState["analysisMethod"],
                    });
                }}
                value={docState.analysisMethod}
              >
                <Label data-header>Analysis Method</Label>
                <div className="cards">
                  <Radio value="parasail">
                    Parasail
                    <p className="text-deemphasis">
                      Best for small datasets. Supports amino acid calculations.
                    </p>
                  </Radio>
                  <Radio value="lzani">
                    LZ-ANI
                    <div style={{}}>
                      <p className="text-deemphasis">
                        Best for large datasets. Only supports nucleotide ANI
                        calculations.
                      </p>
                    </div>
                  </Radio>
                </div>
              </RadioGroup>
            </div>
            <div className="field clustering inline-toggle">
              <Switch
                isSelected={docState.enableClustering}
                onChange={(value) => {
                  updateDocState({ enableClustering: value });
                }}
              >
                Reorder data by linkage clustering method
              </Switch>
              {docState.enableClustering ? (
                <div className="setting clustering-method">
                  <Select
                    data-flat
                    selectedKey={docState.cluster_method}
                    onSelectionChange={(value) => {
                      updateDocState({
                        cluster_method: value as typeof docState.cluster_method,
                      });
                    }}
                    items={Object.entries(reorderMethods).map(
                      ([key, value]) => ({
                        id: key,
                        name: value.name,
                        description: value.description,
                      }),
                    )}
                  >
                    {(item) => (
                      <SelectItem textValue={item.name}>
                        <Text slot="label">{item.name}</Text>
                        <Text slot="description">{item.description}</Text>
                      </SelectItem>
                    )}
                  </Select>
                </div>
              ) : null}
            </div>

            {docState.analysisMethod === "parasail" ? (
              <RunnerPerformance
                appState={appState}
                updateDocState={updateDocState}
                setDocState={setDocState}
                docState={docState}
                initialized={initialized}
              />
            ) : null}

            <div className="actions">
              <Button
                data-primary
                type="button"
                onPress={() => {
                  setStartingRun(true);
                  startRun().finally(() => {
                    setStartingRun(false);
                  });
                }}
                isDisabled={Boolean(
                  startingRun ||
                    appState.active_run_document_id ||
                    !fileName ||
                    docState.validation_error_id,
                )}
              >
                {startingRun ? "Starting..." : "Start Analysis"}
              </Button>
            </div>
          </>
        ) : (
          <>
            <div className="field file-selector">
              <Button type="button" onPress={() => openFileDialog(docState.id)}>
                Select FASTA or SDT Matrix file&#8230;
              </Button>
            </div>
            {appState.recentFiles.some(Boolean) ? (
              <div className="recent-files">
                <h2>Recent Files</h2>
                <div className="grid">
                  {appState.recentFiles.map((file) => (
                    <Button
                      className={"react-aria-Button flat compact"}
                      key={file}
                      onPress={() => openRecentFile(file, docState)}
                    >
                      <TbFile size={16} />
                      <div className="file-info">
                        <h4>{splitFilePath(file).name}</h4>
                        <div className="dir">
                          {splitFilePath(file).dir.replace(
                            appState.config?.userPath || "",
                            "~",
                          )}
                        </div>
                      </div>
                    </Button>
                  ))}
                </div>
              </div>
            ) : null}
          </>
        )}
        {docState.validation_error_id ? (
          <div className="validation-error">
            {messages[docState.validation_error_id]}
          </div>
        ) : null}
      </div>
    </div>
  );
};

export const Runner = ({
  docState,
  setDocState,
  updateDocState,
}: {
  docState: DocState;
  setDocState: SetDocState;
  updateDocState: UpdateDocState;
  tabView: "tabs" | "select";
}) => {
  return (
    <TabPanel
      id={docState.id}
      key={docState.id}
      className="app-panel full-width"
    >
      <div className="app-main centered runner">
        <RunnerSettings
          docState={docState}
          updateDocState={updateDocState}
          setDocState={setDocState}
        />
      </div>
    </TabPanel>
  );
};
