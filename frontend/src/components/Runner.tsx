import React from "react";
import { AppState, SetAppState, clusterMethods } from "../appState";
import messages from "../messages";
import {
  Label,
  Meter,
  Slider,
  SliderOutput,
  SliderThumb,
  SliderTrack,
} from "react-aria-components";
import { formatBytes } from "../helpers";

export type RunProcessDataArgs = Pick<
  AppState["client"],
  "cluster_method" | "compute_cores"
>;

const WarningIcon = () => (
  <svg
    xmlns="http://www.w3.org/2000/svg"
    fillRule="evenodd"
    clipRule="evenodd"
    imageRendering="optimizeQuality"
    shapeRendering="geometricPrecision"
    textRendering="geometricPrecision"
    viewBox="0 0 512 463.43"
    fill="currentColor"
  >
    <path d="M189.46 44.02c34.26-58.66 99.16-58.77 133.24.12l.97 1.81 175.27 304.4c33.71 56.4-1.2 113.76-66.17 112.96v.12H73.53c-.9 0-1.78-.04-2.66-.11-58.34-.79-86.64-54.22-61.9-106.84.39-.85.82-1.67 1.28-2.46l-.04-.03 179.3-309.94-.05-.03zm50.32 302.4c4.26-4.13 9.35-6.19 14.45-6.56 3.4-.24 6.8.29 9.94 1.48 3.13 1.19 6.01 3.03 8.39 5.41 6.92 6.91 8.72 17.38 4.64 26.16-2.69 5.8-7.08 9.7-12.11 11.78-3.03 1.27-6.3 1.84-9.56 1.76-3.27-.08-6.49-.82-9.41-2.18-5.02-2.33-9.3-6.43-11.7-12.2-2.65-6.36-2.27-12.96.63-19.15 1.15-2.46 2.75-4.81 4.73-6.5zm33.86-47.07c-.8 19.91-34.51 19.93-35.28-.01-3.41-34.1-12.13-110.53-11.85-142.58.28-9.87 8.47-15.72 18.94-17.95 3.23-.69 6.78-1.03 10.35-1.02 3.6.01 7.16.36 10.39 1.05 10.82 2.3 19.31 8.39 19.31 18.45l-.05 1-11.81 141.06z" />
  </svg>
);

const RunnerSettings = ({
  appState,
  setAppState,
  startProcessData,
}: {
  appState: AppState;
  setAppState: SetAppState;
  startProcessData: () => void;
}) => {
  const [showOutputSelector, setShowOutputSelector] = React.useState(false);
  const clusterOptionRefs = React.useRef<(HTMLLabelElement | null)[]>([]);

  const handleClusterOptionFocus = (index: number) => {
    clusterOptionRefs.current[index]?.setAttribute("data-focused", "true");
  };

  const handleClusterOptionBlur = (index: number) => {
    clusterOptionRefs.current[index]?.removeAttribute("data-focused");
  };

  const handleChangeClusterMethod = (
    value: typeof appState.client.cluster_method,
  ) =>
    setAppState((previous) => {
      return {
        ...previous,
        client: {
          ...previous.client,
          cluster_method: value,
        },
      };
    });

  const handleChangeComputeCores = (
    value: typeof appState.client.compute_cores,
  ) =>
    setAppState((previous) => {
      return {
        ...previous,
        client: {
          ...previous.client,
          compute_cores: value,
        },
      };
    });

  const fileName =
    appState.filename?.length && appState.filename[0]
      ? appState.filename[0].split("/").pop()
      : undefined;

  const isFastaType = appState.filetype === "text/fasta";

  React.useEffect(() => {
    const handleEnter = (event: KeyboardEvent) => {
      if (
        event.key === "Enter" &&
        (event.ctrlKey || event.metaKey || event.altKey) &&
        Boolean(appState.filename) &&
        !appState.validation_error_id
      ) {
        event.preventDefault();
        startProcessData();
      }
    };
    document.addEventListener("keydown", handleEnter);

    return () => {
      document.removeEventListener("keydown", handleEnter);
    };
  }, []);

  React.useEffect(() => {
    if (!appState.compute_stats) {
      return;
    }
    const stats = appState.compute_stats;
    handleChangeComputeCores(stats.recommended_cores);
  }, [appState.compute_stats?.recommended_cores]);

  React.useEffect(() => {
    const id = setInterval(() => {
      if (!appState.compute_stats) {
        return;
      }

      window.pywebview.api.get_available_memory().then((available_memory) =>
        setAppState((previous) => {
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
  }, [appState.compute_stats]);

  const estimatedMemory =
    (appState.compute_stats?.required_memory || 1) *
    (appState.client.compute_cores || 1);

  const estimatedMemoryValue = appState.compute_stats
    ? ((appState.compute_stats.required_memory *
        (appState.client.compute_cores || 1)) /
        appState.compute_stats.available_memory) *
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
    appState.compute_stats && appState.compute_stats.recommended_cores === 0
      ? "extreme"
      : impactName(
          (appState.client.compute_cores / appState.platform.cores) * 100,
        );

  return (
    <div className="form-wrapper runner-wrapper">
      <div className="form">
        <div className="field">
          <label className="header">FASTA or SDT Matrix File</label>
          <div className="input-with-button">
            <input
              type="text"
              readOnly
              value={appState.validation_error_id ? "" : fileName}
            />
            <button
              type="button"
              onClick={() => window.pywebview.api.open_file_dialog()}
            >
              Select file...
            </button>
          </div>
        </div>
        {isFastaType && !appState.validation_error_id ? (
          <>
            <div className="field runner-settings">
              <label className="header">Clustering Method</label>
              <div className="clustering-method">
                {clusterMethods.map((value, index) => (
                  <label
                    className="radio"
                    key={value}
                    ref={(el) => (clusterOptionRefs.current[index] = el)}
                    data-selected={appState.client.cluster_method === value}
                  >
                    <input
                      key={value}
                      type="radio"
                      id={value}
                      name="cluster-method"
                      value={value}
                      checked={appState.client.cluster_method === value}
                      onChange={() => handleChangeClusterMethod(value)}
                      onFocus={() => handleClusterOptionFocus(index)}
                      onBlur={() => handleClusterOptionBlur(index)}
                    />
                    <span>{value}</span>
                  </label>
                ))}
              </div>
            </div>

            <div className="field runner-settings performance">
              <label className="header" htmlFor="compute-cores">
                Performance
              </label>
              {appState.compute_stats ? (
                <div className="performance-settings col-2">
                  <div className="cores-used">
                    <Slider
                      id="compute-cores"
                      onChange={handleChangeComputeCores}
                      minValue={1}
                      maxValue={appState.platform.cores}
                      value={appState.client.compute_cores}
                    >
                      <Label>Cores</Label>
                      <SliderOutput data-impact={coresImpact}>
                        {({ state }) => (
                          <>
                            {appState.compute_stats &&
                            (appState.compute_stats.recommended_cores === 0 ||
                              appState.client.compute_cores >
                                appState.compute_stats.recommended_cores) ? (
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
                                width: state.getThumbPercent(0) * 100 + "%",
                              }}
                            />
                            <SliderThumb />
                          </>
                        )}
                      </SliderTrack>
                    </Slider>
                    <small>
                      Recommended: {appState.compute_stats.recommended_cores}
                    </small>
                  </div>
                  <div className="estimated-memory">
                    <Meter value={estimatedMemoryValue}>
                      {({ percentage }) => (
                        <>
                          <Label>Memory</Label>
                          <span className="value">
                            {formatBytes(estimatedMemory, 1)} /{" "}
                            {appState.compute_stats?.available_memory
                              ? `${formatBytes(
                                  appState.compute_stats?.available_memory || 0,
                                  0,
                                )} Available`
                              : null}
                          </span>
                          <div className="bar">
                            <div
                              className="fill"
                              data-impact={impactName(percentage)}
                              style={{
                                width: percentage + "%",
                                minWidth: "4px",
                              }}
                            />
                          </div>
                        </>
                      )}
                    </Meter>
                    <small>
                      <span>
                        {appState.compute_stats &&
                        estimatedMemoryValue > 100 ? (
                          <>
                            Swap:{" "}
                            {formatBytes(
                              estimatedMemory -
                                appState.compute_stats.available_memory,
                            )}
                          </>
                        ) : null}
                      </span>
                      <span>
                        {formatBytes(appState.platform.memory, 0)} Total
                      </span>
                    </small>
                  </div>
                </div>
              ) : null}
            </div>

            <div className="advanced-settings">
              <label>
                <input
                  type="checkbox"
                  checked={showOutputSelector}
                  onChange={() => {
                    setShowOutputSelector(!showOutputSelector);
                    if (!showOutputSelector) {
                      setAppState((previous) => {
                        return { ...previous, alignment_output_path: "" };
                      });
                    }
                  }}
                />
                Output alignments to folder
              </label>

              {showOutputSelector ? (
                <div className="field">
                  <div className="input-with-button">
                    <input
                      type="text"
                      value={appState.alignment_output_path}
                      readOnly
                    />
                    <button
                      type="button"
                      onClick={() =>
                        window.pywebview.api.select_alignment_output_path()
                      }
                    >
                      Select folder...
                    </button>
                  </div>
                </div>
              ) : null}
            </div>

            {appState.compute_stats &&
            (appState.compute_stats?.recommended_cores === 0 ||
              estimatedMemory > appState.platform.memory ||
              appState.client.compute_cores >
                appState.compute_stats?.recommended_cores) ? (
              <div className="compute-forecast">
                <b>Warning:</b> System instability may occur when exceeding
                recommended cores or available memory. When exceeding available
                memory, it is recommended to limit cores used.
              </div>
            ) : null}

            <div className="actions">
              <button
                type="button"
                onClick={startProcessData}
                disabled={!fileName || appState.validation_error_id}
              >
                Run
              </button>
            </div>
          </>
        ) : null}
        {appState.validation_error_id ? (
          <div className="validation-error">
            {messages[appState.validation_error_id]}
          </div>
        ) : null}
      </div>
    </div>
  );
};

export const Runner = ({
  appState,
  setAppState,
  mainMenu,
  startProcessData,
}: {
  appState: AppState;
  setAppState: SetAppState;
  mainMenu: React.ReactNode;
  startProcessData: () => void;
}) => {
  const [appConfig, setAppConfig] = React.useState<{ appVersion: string }>();
  const fetchAppConfig = () => {
    window.pywebview.api
      .app_config()
      .then((result) => setAppConfig(JSON.parse(result)));
  };

  React.useEffect(() => {
    const waitForPywebview = () =>
      new Promise((resolve) => {
        if (!window.pywebview) {
          setTimeout(waitForPywebview, 10);
        } else {
          resolve(true);
        }
      });

    waitForPywebview().then(() => fetchAppConfig());
  }, []);

  return (
    <div className="app-wrapper with-header with-footer">
      <div className="app-header runner">{mainMenu}</div>
      <div className="app-main centered runner">
        <RunnerSettings
          appState={appState}
          setAppState={setAppState}
          startProcessData={startProcessData}
        />
      </div>
      <div className="app-footer centered">
        <div>{appConfig ? appConfig.appVersion : null}</div>
      </div>
    </div>
  );
};
