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
} from "react-aria-components";
import useAppState, {
  type AppState,
  type DocState,
  type SetDocState,
  type UpdateDocState,
  clusterMethods,
} from "../appState";
import { formatBytes } from "../helpers";
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

const WarningIcon = () => (
  <svg
    aria-hidden="true"
    xmlns="http://www.w3.org/2000/svg"
    fillRule="evenodd"
    clipRule="evenodd"
    imageRendering="optimizeQuality"
    shapeRendering="geometricPrecision"
    textRendering="geometricPrecision"
    viewBox="0 0 512 463.43"
    fill="currentColor"
  >
    <path d="M189.46 44.02c34.26-58.66 99.16-58.77 133.24.12l.97 1.81 175.27 304.4c33.71 56.4-1.2 113.76-66.17 112.96v.12H73.53c-.9 0-1.78-.04-2.66-.11-58.34-.79-86.64-54.22-61.9-106.84.39-.85.82-1.67 1.28-2.46l-.04-.03 179.3-309.94-.05-.03zm50.32 302.4c4.26-4.13 9.35-6.19 14.45-6.56 3.4-.24 6.8.29 9.94 1.48 3.13 1.19 6.01 3.03 8.39 5.41 6.92 6.91 8.72 17.38 4.64 26.16-2.69 5.8-7.08 9.7-12.11 11.78-3.03 1.27-6.3 1.84-9.56 1.76-3.27-.08-6.49-.82-9.41-2.18-5.02-2.33-9.3-6.43-11.7-12.2-2.65-6.36-2.27-12.96.63-19.15 1.15-2.46 2.75-4.81 4.73-6.5zm33.86-47.07c-.8 19.91-34.51 19.93-35.28-.0-3.41-34.1-12.13-110.53-11.85-142.58.28-9.87 8.47-15.72 18.94-17.95 3.23-.69 6.78-1.03 10.35-1.02 3.6.01 7.16.36 10.39 1.05 10.82 2.3 19.31 8.39 19.31 18.45l-.05 1-11.81 141.06z" />
  </svg>
);

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
            <Button type="button" onPress={() => openFileDialog(docState.id)}>
              Select file&#8230;
            </Button>
          </div>
        </div>
        {isFastaType && !docState.validation_error_id ? (
          <>
            <div className="field clustering inline-toggle">
              <Switch
                isSelected={appState.enableClustering}
                onChange={(value) => {
                  updateAppState({ enableClustering: value });
                }}
              >
                Cluster sequences
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
                      cluster_method: value as (typeof clusterMethods)[number],
                    });
                  }}
                  items={clusterMethods.map((name) => ({
                    id: name,
                    name,
                  }))}
                >
                  {(item) => (
                    <SelectItem textValue={item.name}>
                      {item.name} method
                    </SelectItem>
                  )}
                </Select>
              </div>
            </div>

            <div className="field output inline-toggle">
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
            </div>

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
                              <WarningIcon />
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
                onPress={startRun}
                isDisabled={Boolean(
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
          <div className="recent-files">
            <h2>Recent Files</h2>
            {recentFiles.map((file) => (
              <Button
                className={"react-aria-Button flat compact"}
                key={file}
                onPress={() => openFile(file, docState.id)}
              >
                {file.split(/(\/|\\)/).pop()}
              </Button>
            ))}
          </div>
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
    <TabPanel id={docState.id} key={docState.id}>
      <div className="app-main full-height centered runner">
        <RunnerSettings
          docState={docState}
          updateDocState={updateDocState}
          setDocState={setDocState}
        />
      </div>
    </TabPanel>
  );
};
