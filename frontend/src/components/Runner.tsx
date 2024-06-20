import React from "react";
import {
  AppState,
  SetAppState,
  performanceProfiles,
  alignmentTypes,
  clusterMethods,
} from "../appState";
import { NumberInput } from "./NumberInput";

export type RunSDT2Args = Pick<
  AppState["client"],
  "performance_profile" | "alignment_type" | "cluster_method"
>;

const RunnerSettings = ({
  appState,
}: {
  appState: AppState;

  setAppState: SetAppState;
}) => {
  const [runnerSettings, setRunnerSettings] = React.useState<RunSDT2Args>(
    (({ client: { performance_profile, alignment_type, cluster_method } }) => ({
      performance_profile,
      alignment_type,
      cluster_method,
    }))(appState),
  );

  const [advancedSettingsState, setAdvancedSettingsState] = React.useState({
    match: 1.5,
    mismatch: -1,
    iog: -2,
    ieg: -2,
    log: -3,
    leg: -2,
    rog: -3,
    reg: -2,
  });

  const handleRun = () => {
    window.pywebview.api.run_sdt2({
      cluster_method: runnerSettings.cluster_method,
      performance_profile: runnerSettings.performance_profile,
      alignment_type: runnerSettings.alignment_type,
      ...advancedSettingsState,
    });
  };

  const handleChangeAdvancedSetting =
    (field: keyof typeof advancedSettingsState) => (value: number) =>
      setAdvancedSettingsState({
        ...advancedSettingsState,
        [field]: value,
      });

  const handleChangePerformanceProfile = (
    value: typeof appState.client.performance_profile,
  ) =>
    setRunnerSettings({
      ...runnerSettings,
      performance_profile: value,
    });

  const handleChangeAlignmentType = (
    value: typeof appState.client.alignment_type,
  ) =>
    setRunnerSettings({
      ...runnerSettings,
      alignment_type: value,
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

  return (
    <>
      <div className="app-main centered runner">
        <div className="form-wrapper runner-wrapper">
          <div className="form">
            <div className="field">
              <label className="header">FASTA or SDT Matrix File</label>
              <div className="input-with-button">
                <input type="text" readOnly value={fileName} />
                <button
                  type="button"
                  onClick={() => window.pywebview.api.open_file_dialog()}
                >
                  Select file...
                </button>
              </div>
            </div>
            {isFastaType ? (
              <>
                <details className="advanced-settings">
                  <summary>
                    <strong>Advanced</strong>
                  </summary>
                  <div className="group">
                    <h4>Aligner</h4>
                    <div className="col-2">
                      <NumberInput
                        label="Match Score"
                        type="float"
                        onChange={handleChangeAdvancedSetting("match")}
                        value={advancedSettingsState.match}
                        min={-50}
                        max={50}
                        step={0.5}
                      />
                      <NumberInput
                        label="Mismatch Score"
                        onChange={handleChangeAdvancedSetting("mismatch")}
                        value={advancedSettingsState.mismatch}
                        min={-50}
                        max={50}
                        step={1}
                      />
                    </div>
                  </div>

                  <div className="group">
                    <h4>Internal</h4>
                    <div className="col-2">
                      <NumberInput
                        label="Open Gap Score"
                        onChange={handleChangeAdvancedSetting("iog")}
                        value={advancedSettingsState.iog}
                        min={-50}
                        max={50}
                        step={1}
                      />
                      <NumberInput
                        type="int"
                        label="Extend Gap Score"
                        onChange={handleChangeAdvancedSetting("ieg")}
                        min={-50}
                        max={50}
                        step={1}
                        value={advancedSettingsState.ieg}
                      />
                    </div>
                  </div>

                  <div className="group">
                    <h4>Left</h4>
                    <div className="col-2">
                      <NumberInput
                        type="int"
                        label="Open Gap Score"
                        onChange={handleChangeAdvancedSetting("log")}
                        min={-50}
                        max={50}
                        step={1}
                        value={advancedSettingsState.log}
                      />
                      <NumberInput
                        type="int"
                        label="Extend Gap Score"
                        onChange={handleChangeAdvancedSetting("leg")}
                        min={-50}
                        max={50}
                        step={1}
                        value={advancedSettingsState.leg}
                      />
                    </div>
                  </div>

                  <div className="group">
                    <h4>Right</h4>
                    <div className="col-2">
                      <NumberInput
                        type="int"
                        label="Open Gap Score"
                        onChange={handleChangeAdvancedSetting("rog")}
                        min={-50}
                        max={50}
                        step={1}
                        value={advancedSettingsState.rog}
                      />
                      <NumberInput
                        type="int"
                        label="Extend Gap Score"
                        onChange={handleChangeAdvancedSetting("reg")}
                        min={-50}
                        max={50}
                        step={1}
                        value={advancedSettingsState.reg}
                      />
                    </div>
                  </div>

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
                  <div className="field runner-settings">
                    <label className="header">Alignment Type</label>
                    {alignmentTypes.map((value) => (
                      <label className="radio" key={value}>
                        <input
                          type="radio"
                          id={value}
                          name="alignment-type"
                          value={value}
                          checked={runnerSettings.alignment_type === value}
                          onChange={() => handleChangeAlignmentType(value)}
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
                        handleChangePerformanceProfile(
                          event.target.value as any,
                        )
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
                    disabled={!Boolean(fileName) ?? false}
                  >
                    Run
                  </button>
                </div>
              </>
            ) : null}
            {fileName && !isFastaType ? (
              <div className="incompatible-filetype">
                This file is not compatible with SDT2. Please select a file with
                an extension of .fasta or .csv.
              </div>
            ) : null}
          </div>
        </div>
      </div>
    </>
  );
};

export const Runner = ({
  appState,
  setAppState,
}: {
  appState: AppState;
  setAppState: SetAppState;
}) => (
  <div className="app-wrapper">
    <RunnerSettings appState={appState} setAppState={setAppState} />
  </div>
);
