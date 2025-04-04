import React from "react";
import {
  Button,
  Input,
  Label,
  Meter,
  Slider,
  SliderOutput,
  SliderThumb,
  SliderTrack,
  TabPanel,
  Text,
} from "react-aria-components";
import { TbAlertTriangleFilled, TbFile } from "react-icons/tb";
import useAppState, {
  type AppState,
  type DocState,
  type SetDocState,
  type UpdateDocState,
} from "../appState";
import { reorderMethods } from "../constants";
import { formatBytes, splitFilePath } from "../helpers";
import useOpenFileDialog from "../hooks/useOpenFileDialog";
import { useStartRun } from "../hooks/useStartRun";
import messages from "../messages";
import { openFile } from "../services/files";
import { Select, SelectItem } from "./Select";
import { Switch } from "./Switch";

export type RunProcessDataArgs = Pick<AppState, "compute_cores"> & {
  doc_id: string;
  cluster_method: AppState["cluster_method"] | "None";
  export_alignments: "True" | "False";
};

const RunnerSettings = ({
  docState,
  setDocState,
}: {
  docState: DocState;
  setDocState: SetDocState;
  updateDocState: UpdateDocState;
}) => {
  const { appState, setAppState } = useAppState();
  const startRun = useStartRun(docState, appState);
  const [startingRun, setStartingRun] = React.useState(false);
  const openFileDialog = useOpenFileDialog(appState, setAppState);
  const fileName =
    docState.filename?.length && docState.filename
      ? docState.filename.split("/").pop()
      : undefined;

  const updateAppState = React.useCallback(
    (value: Partial<AppState>) =>
      setAppState((previous) => {
        return {
          ...previous,
          ...value,
        };
      }),
    [setAppState],
  );

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

  const initialized = React.useRef(false);
  React.useEffect(() => {
    if (!docState.compute_stats || initialized.current) {
      return;
    }
    updateAppState({
      compute_cores: docState.compute_stats.recommended_cores,
    });
    initialized.current = true;
  }, [docState.compute_stats, updateAppState]);

  React.useEffect(() => {
    const id = setInterval(() => {
      if (!docState.filename) {
        return;
      }

      window.pywebview.api.get_available_memory().then((available_memory) =>
        setDocState((previous) => {
          return {
            ...previous,
            compute_stats: {
              ...(previous.compute_stats || {
                recommended_cores: 1,
                required_memory: 1,
              }),
              available_memory,
            },
          };
        }),
      );
    }, 3000);

    return () => clearInterval(id);
  }, [docState.filename, setDocState]);

  const [recentFiles, setRecentFiles] = React.useState<string[]>([]);
  React.useEffect(() => {
    window.pywebview.api.app_settings().then((data) => {
      setRecentFiles(data.recent_files);
    });
  }, []);

  const estimatedMemory =
    (docState.compute_stats?.required_memory || 1) *
    (appState.compute_cores || 1);

  const estimatedMemoryValue = docState.compute_stats
    ? ((docState.compute_stats.required_memory *
        (appState.compute_cores || 1)) /
        docState.compute_stats.available_memory) *
      100
    : 1;

  const impactName = (percentage: number) =>
    percentage > 99
      ? "extreme"
      : percentage > 90
        ? "high"
        : percentage > 80
          ? "medium"
          : "low";

  const coresImpact =
    docState.compute_stats && docState.compute_stats.recommended_cores === 0
      ? "extreme"
      : impactName((appState.compute_cores / appState.platform.cores) * 100);

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
            <div className="field clustering inline-toggle">
              <Switch
                isSelected={appState.enableClustering}
                onChange={(value) => {
                  updateAppState({ enableClustering: value });
                }}
              >
                Reorder sequences by cluster method:
              </Switch>
              <div
                className="setting clustering-method"
                data-hidden={!appState.enableClustering}
                aria-hidden={!appState.enableClustering}
              >
                <Select
                  selectedKey={appState.cluster_method}
                  onSelectionChange={(value) => {
                    updateAppState({
                      cluster_method: value as typeof appState.cluster_method,
                    });
                  }}
                  items={Object.entries(reorderMethods).map(([key, value]) => ({
                    id: key,
                    name: value.name,
                    description: value.description,
                  }))}
                >
                  {(item) => (
                    <SelectItem textValue={item.name}>
                      <Text slot="label">{item.name}</Text>
                      <Text slot="description">{item.description}</Text>
                    </SelectItem>
                  )}
                </Select>
              </div>
            </div>

            {/* <div className="field output inline-toggle">
              <Switch
                isSelected={appState.enableOutputAlignments}
                onChange={(value) => {
                  updateAppState({ enableOutputAlignments: value });
                }}
              >
                Save alignments
              </Switch>

              {appState.enableOutputAlignments ? (
                <div className="setting">
                  {appState.alignmentExportPath ? (
                    <div aria-live="polite" className="folder">
                      <span>
                        {appState.alignmentExportPath.replace("/", "")}
                      </span>
                    </div>
                  ) : null}
                  <Button
                    onPress={() => {
                      window.pywebview.api
                        .select_path_dialog(appState.alignmentExportPath)
                        .then((result) => {
                          if (!result) {
                            return;
                          }
                          setAppState((prev) => ({
                            ...prev,
                            alignmentExportPath: result,
                          }));
                        });
                    }}
                  >
                    Set folder&#8230;
                  </Button>
                </div>
              ) : null}
            </div> */}

            <div className="field performance">
              <label className="header" htmlFor="compute-cores">
                Performance
              </label>
              {docState.compute_stats ? (
                <div className="setting performance-settings">
                  <div className="cores-used inline">
                    <div>
                      <Label id="compute-cores-label">Cores</Label>
                      <small className="text-deemphasis">
                        {docState.compute_stats.recommended_cores} Recommended /{" "}
                        {appState.platform.cores} Total
                      </small>
                    </div>
                    <Slider
                      aria-labelledby="compute-cores-label"
                      onChange={(value) =>
                        updateAppState({ compute_cores: value })
                      }
                      minValue={1}
                      maxValue={appState.platform.cores}
                      value={appState.compute_cores}
                    >
                      <SliderOutput data-impact={coresImpact}>
                        {({ state }) => (
                          <>
                            {docState.compute_stats &&
                            (docState.compute_stats.recommended_cores === 0 ||
                              appState.compute_cores >
                                docState.compute_stats.recommended_cores) ? (
                              <TbAlertTriangleFilled />
                            ) : null}
                            {state.getThumbValueLabel(0)} /{" "}
                            {appState.platform.cores}
                          </>
                        )}
                      </SliderOutput>
                      <SliderTrack data-impact={coresImpact}>
                        {({ state }) => (
                          <>
                            <div className="track" />
                            <div
                              className="fill"
                              style={{
                                width: `${state.getThumbPercent(0) * 100}%`,
                              }}
                            />
                            <SliderThumb />
                          </>
                        )}
                      </SliderTrack>
                    </Slider>
                  </div>
                  <div className="memory-used inline">
                    <div>
                      <Label id="memory-used-label" htmlFor="meter">
                        Memory
                      </Label>
                      <div>
                        <small className="text-deemphasis">
                          {docState.compute_stats?.available_memory
                            ? `${formatBytes(
                                docState.compute_stats?.available_memory || 0,
                                0,
                              )} Available`
                            : null}{" "}
                          / {formatBytes(appState.platform.memory)} Total
                        </small>
                      </div>
                    </div>
                    <div className="estimated-memory">
                      <Meter
                        value={estimatedMemoryValue}
                        aria-labelledby="memory-used-label"
                      >
                        {({ percentage }) => (
                          <>
                            <span className="value">
                              {formatBytes(estimatedMemory, 1)} /{" "}
                              {docState.compute_stats?.available_memory
                                ? `${formatBytes(
                                    docState.compute_stats?.available_memory ||
                                      0,
                                    1,
                                  )}`
                                : null}
                            </span>
                            <div className="bar">
                              <div
                                className="fill"
                                data-impact={impactName(percentage)}
                                style={{
                                  width: `${percentage}%`,
                                  minWidth: "4px",
                                }}
                              />
                            </div>
                          </>
                        )}
                      </Meter>
                      <small className="swap">
                        {docState.compute_stats &&
                        estimatedMemoryValue > 100 ? (
                          <>
                            {formatBytes(
                              estimatedMemory -
                                docState.compute_stats.available_memory,
                              0,
                            )}{" "}
                            Swap
                          </>
                        ) : null}
                      </small>
                    </div>
                  </div>
                </div>
              ) : null}
            </div>

            {docState.compute_stats &&
            (docState.compute_stats?.recommended_cores === 0 ||
              estimatedMemory > appState.platform.memory ||
              appState.compute_cores >
                docState.compute_stats?.recommended_cores) ? (
              <div className="compute-forecast">
                <b>Warning:</b> System instability may occur when exceeding
                recommended cores or available memory. When exceeding available
                memory, it is recommended to limit cores used.
              </div>
            ) : null}
            <div className="actions">
              <Button
                data-primary
                type="button"
                onPress={() => {
                  setStartingRun(true);
                  startRun();
                }}
                isDisabled={Boolean(
                  startingRun ||
                    appState.active_run_document_id ||
                    !fileName ||
                    docState.validation_error_id ||
                    (appState.enableOutputAlignments &&
                      !appState.alignmentExportPath.length),
                )}
              >
                Start Analysis
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
            {recentFiles.some(Boolean) ? (
              <div className="recent-files">
                <h2>Recent Files</h2>
                <div className="grid">
                  {recentFiles.map((file) => (
                    <Button
                      className={"react-aria-Button flat compact"}
                      key={file}
                      onPress={() => openFile(file, docState.id)}
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
