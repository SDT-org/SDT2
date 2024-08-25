import React from "react";
import {
  AppState,
  SetAppState,
  performanceProfiles,
  clusterMethods,
} from "../appState";
import { NumberInput } from "./NumberInput";
import messages from "../messages";

export type RunProcessDataArgs = Pick<
  AppState["client"],
  "performance_profile" | "cluster_method"
>;

const RunnerSettings = ({
  appState,
}: {
  appState: AppState;

  setAppState: SetAppState;
}) => {
  const [runnerSettings, setRunnerSettings] =
    React.useState<RunProcessDataArgs>(
      (({ client: { performance_profile, cluster_method } }) => ({
        performance_profile,
        cluster_method,
      }))(appState),
    );

  const handleRun = () => {
    window.pywebview.api.run_process_data({
      cluster_method: runnerSettings.cluster_method,
      performance_profile: runnerSettings.performance_profile,
    });
  };

  const handleChangePerformanceProfile = (
    value: typeof appState.client.performance_profile,
  ) =>
    setRunnerSettings({
      ...runnerSettings,
      performance_profile: value,
    });

  const handleChangeClusterMethod = (
    value: typeof appState.client.cluster_method,
  ) =>
    setRunnerSettings({
      ...runnerSettings,
      cluster_method: value,
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
                      checked={runnerSettings.cluster_method === value}
                      onChange={() => handleChangeClusterMethod(value)}
                    />
                    <span>{value}</span>
                  </label>
                ))}
              </div>
            </div>
            <div className="field runner-settings performance">
              <label className="header">Compute Mode</label>
              <div className="col-2">
                <select
                  className="="
                  onChange={(event) =>
                    handleChangePerformanceProfile(event.target.value as any)
                  }
                  value={runnerSettings.performance_profile as string}
                >
                  {Object.keys(performanceProfiles).map((key) => (
                    <option key={key} value={key}>
                      {
                        performanceProfiles[
                          key as keyof typeof performanceProfiles
                        ]
                      }
                    </option>
                  ))}
                </select>
                <div className="field runner-settings performance-details">
                  <span className="cores">
                    {
                      appState.performance_profiles[
                        runnerSettings.performance_profile
                      ]
                    }
                    <span>/</span>
                    {appState.performance_profiles["best"]}
                  </span>
                  cores
                </div>
              </div>
            </div>

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
}: {
  appState: AppState;
  setAppState: SetAppState;
}) => {
  const [appConfig, setAppConfig] = React.useState<{ appVersion: string }>();
  const fetchAppConfig = () => {
    window.pywebview.api
      .app_config()
      .then((result) => setAppConfig(JSON.parse(result)));
  };

  React.useEffect(() => {
    if (!window.pywebview) {
      setTimeout(fetchAppConfig, 500);
      return;
    }
    fetchAppConfig();
  }, []);

  return (
    <div className="app-wrapper with-footer">
      <div className="app-main centered runner">
        <RunnerSettings appState={appState} setAppState={setAppState} />
      </div>
      <div className="app-footer centered">
        <div>{appConfig ? appConfig.appVersion : null}</div>
      </div>
    </div>
  );
};
