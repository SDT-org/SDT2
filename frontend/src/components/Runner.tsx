import React from "react";
import {
  AppState,
  SetAppState,
  clusterMethods,
  PerformanceProfile,
} from "../appState";
import messages from "../messages";
import {
  Label,
  Select,
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

const RunnerSettings = ({
  appState,
  setAppState,
}: {
  appState: AppState;
  setAppState: SetAppState;
}) => {
  const [computeModes, setComputeModes] = React.useState<
    Record<Exclude<PerformanceProfile, "custom">, number>
  >({
    recommended: 1,
    best: 1,
    balanced: 1,
    low: 1,
  });

  const handleRun = () => {
    window.pywebview.api.run_process_data({
      cluster_method: appState.client.cluster_method,
      compute_cores: appState.client.compute_cores,
    });
  };

  const handleChangePerformanceProfile = (value: PerformanceProfile) =>
    setAppState((previous) => ({
      ...previous,
      client: {
        ...previous.client,
        performanceProfile: value,
        ...(value !== "custom" && { compute_cores: computeModes[value] }),
      },
    }));

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
        Boolean(appState.filename) &&
        !appState.validation_error_id
      ) {
        event.preventDefault();
        handleRun();
      }
    };
    document.addEventListener("keydown", handleEnter);

    return () => {
      document.removeEventListener("keydown", handleEnter);
    };
  }, []);

  React.useEffect(() => {
    if (appState.client.performanceProfile !== "custom") {
      handleChangeComputeCores(
        computeModes[appState.client.performanceProfile],
      );
    }
  }, [
    computeModes,
    appState.filename,
    appState.compute_stats,
    appState.client.performanceProfile,
  ]);

  React.useEffect(() => {
    if (!appState.compute_stats) {
      return;
    }
    const stats = appState.compute_stats;
    const platform = appState.platform;

    setAppState((previous) => ({
      ...previous,
      client: {
        ...previous.client,
        compute_cores: stats.recommended_cores,
        performanceProfile: "recommended",
      },
    }));

    setComputeModes({
      recommended: stats.recommended_cores,
      best: platform.cores,
      balanced: Math.floor(Math.max(platform.cores / 2, 1)),
      low: 1,
    });
  }, [appState.compute_stats]);

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
            <details className="advanced-settings">
              <summary>
                <strong>Advanced</strong>
              </summary>

              <div className="group">
                <div className="field">
                  <label className="header">Alignment Output Folder</label>
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
                      Select...
                    </button>
                  </div>
                </div>
              </div>
            </details>

            <div className="col-2">
              <div className="field runner-settings">
                <label className="header">Clustering Method</label>
                {clusterMethods.map((value) => (
                  <label className="radio" key={value}>
                    <input
                      key={value}
                      type="radio"
                      id={value}
                      name="cluster-method"
                      value={value}
                      checked={appState.client.cluster_method === value}
                      onChange={() => handleChangeClusterMethod(value)}
                    />
                    <span>{value}</span>
                  </label>
                ))}
              </div>

              <div className="field runner-settings performance">
                <label className="header">Compute Mode</label>
                <select
                  onChange={(event) =>
                    handleChangePerformanceProfile(
                      event.target.value as PerformanceProfile,
                    )
                  }
                  value={appState.client.performanceProfile}
                >
                  {Object.keys(computeModes).map((key) => (
                    <option key={key} value={key}>
                      {key}
                    </option>
                  ))}
                  <option key="custom" value="custom">
                    Custom
                  </option>
                </select>
                <div className="performance-details">
                  {appState.client.performanceProfile === "custom" &&
                  appState.compute_stats ? (
                    <>
                      <Slider
                        onChange={handleChangeComputeCores}
                        minValue={1}
                        maxValue={appState.platform.cores}
                        defaultValue={appState.client.compute_cores}
                      >
                        <Label>Cores</Label>
                        <SliderOutput />
                        <SliderTrack>
                          <SliderThumb />
                        </SliderTrack>
                      </Slider>
                    </>
                  ) : (
                    <>
                      <span className="cores">
                        {appState.client.performanceProfile === "custom"
                          ? appState.client.compute_cores
                          : computeModes[appState.client.performanceProfile]}
                        <span>/</span>
                        {appState.platform.cores}
                      </span>
                      cores
                    </>
                  )}
                </div>
              </div>
            </div>

            {appState.compute_stats &&
            (appState.client.compute_cores *
              appState.compute_stats.required_memory >
              appState.platform.memory ||
              appState.client.compute_cores >
                appState.compute_stats?.recommended_cores) ? (
              <div className="compute-forecast">
                <p>
                  <b>Warning:</b> Analysing these sequences may cause system
                  instability if you select more than the recommended amout of
                  cores.
                </p>
                <p>
                  <b>Recommended cores:</b>{" "}
                  {appState.compute_stats.recommended_cores}
                </p>
                <p>
                  <b>Required memory:</b>
                  <br />
                  {formatBytes(
                    appState.compute_stats.required_memory *
                      appState.client.compute_cores,
                  )}{" "}
                  / {formatBytes(appState.platform.memory)}{" "}
                  <small>
                    ({formatBytes(appState.compute_stats.required_memory)} per
                    core)
                  </small>
                </p>
              </div>
            ) : null}

            <div className="actions">
              <button
                type="button"
                onClick={handleRun}
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
}: {
  appState: AppState;
  setAppState: SetAppState;
  mainMenu: React.ReactNode;
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
        <RunnerSettings appState={appState} setAppState={setAppState} />
      </div>
      <div className="app-footer centered">
        <div>{appConfig ? appConfig.appVersion : null}</div>
      </div>
    </div>
  );
};
