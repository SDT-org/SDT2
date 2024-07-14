import React from "react";
import {
  AppState,
  SetAppState,
  performanceProfiles,
  alignmentTypes,
} from "../appState";

export type RunSDT2Args = Pick<
  AppState["client"],
  "neighbor_joining" | "performance_profile" | "alignment_type"
>;

const RunnerSettings = ({
  appState,
}: {
  appState: AppState;

  setAppState: SetAppState;
}) => {
  const [runnerSettings, setRunnerSettings] = React.useState<RunSDT2Args>(
    (({
      client: { neighbor_joining, performance_profile, alignment_type },
    }) => ({
      neighbor_joining,
      performance_profile,
      alignment_type,
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
      neighbor_joining: runnerSettings.neighbor_joining,
      performance_profile: runnerSettings.performance_profile,
      alignment_type: runnerSettings.alignment_type,
      ...advancedSettingsState,
    });
  };

  const toggleNeighborJoining = () =>
    setRunnerSettings({
      ...runnerSettings,
      neighbor_joining: !runnerSettings.neighbor_joining,
    });

  const handleChangePerformanceProfile = (
    value: (typeof performanceProfiles)[keyof typeof performanceProfiles],
  ) =>
    setRunnerSettings({
      ...runnerSettings,
      performance_profile: value,
    });

  const handleChangeAlignmentType = (
    value: (typeof alignmentTypes)[keyof typeof alignmentTypes],
  ) =>
    setRunnerSettings({
      ...runnerSettings,
      alignment_type: value,
    });

  const fileName =
    appState.filename?.length && appState.filename[0]
      ? appState.filename[0].split("/").pop()
      : undefined;

  const isFastaType = appState.filetype === "text/fasta";

  return (
    <>
      <div className="app-main centered">
        <div className="form-wrapper runner-wrapper">
          <div className="form">
            <div className="field">
              <label className="header">FASTA or SDT Matrix File</label>
              <input type="text" readOnly value={fileName} />
              <button
                type="button"
                onClick={() => window.pywebview.api.open_file_dialog()}
              >
                Select file
              </button>
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
                      <div className="field">
                        <label>Match Score</label>
                        <input
                          type="number"
                          onChange={(v) =>
                            setAdvancedSettingsState({
                              ...advancedSettingsState,
                              match: parseFloat(v.target.value),
                            })
                          }
                          value={advancedSettingsState.match}
                          min={-50}
                          max={50}
                          step={0.5}
                        />
                      </div>
                      <div className="field">
                        <label>Mismatch Score</label>
                        <input
                          type="number"
                          onChange={(v) =>
                            setAdvancedSettingsState({
                              ...advancedSettingsState,
                              mismatch: parseFloat(v.target.value),
                            })
                          }
                          value={advancedSettingsState.mismatch}
                          min={-50}
                          max={50}
                          step={1}
                        />
                      </div>
                    </div>
                  </div>

                  <div className="group">
                    <h4>Internal</h4>
                    <div className="col-2">
                      <div className="field">
                        <label>Open Gap Score</label>
                        <input
                          type="number"
                          onChange={(v) =>
                            setAdvancedSettingsState({
                              ...advancedSettingsState,
                              iog: parseFloat(v.target.value),
                            })
                          }
                          value={advancedSettingsState.iog}
                          min={-50}
                          max={50}
                          step={1}
                        />
                      </div>
                      <div className="field">
                        <label>Extend Gap Score</label>
                        <input
                          type="number"
                          onChange={(v) =>
                            setAdvancedSettingsState({
                              ...advancedSettingsState,
                              ieg: parseFloat(v.target.value),
                            })
                          }
                          value={advancedSettingsState.ieg}
                          min={-50}
                          max={50}
                          step={1}
                        />
                      </div>
                    </div>
                  </div>

                  <div className="group">
                    <h4>Left</h4>
                    <div className="col-2">
                      <div className="field">
                        <label>Open Gap Score</label>
                        <input
                          type="number"
                          onChange={(v) =>
                            setAdvancedSettingsState({
                              ...advancedSettingsState,
                              log: parseFloat(v.target.value),
                            })
                          }
                          value={advancedSettingsState.log}
                          min={-50}
                          max={50}
                          step={1}
                        />
                      </div>
                      <div className="field">
                        <label>Extend Gap Score</label>
                        <input
                          type="number"
                          onChange={(v) =>
                            setAdvancedSettingsState({
                              ...advancedSettingsState,
                              leg: parseFloat(v.target.value),
                            })
                          }
                          value={advancedSettingsState.leg}
                          min={-50}
                          max={50}
                          step={1}
                        />
                      </div>
                    </div>
                  </div>

                  <div className="group">
                    <h4>Right</h4>
                    <div className="col-2">
                      <div className="field">
                        <label>Open Gap Score</label>
                        <input
                          type="number"
                          onChange={(v) =>
                            setAdvancedSettingsState({
                              ...advancedSettingsState,
                              rog: parseFloat(v.target.value),
                            })
                          }
                          value={advancedSettingsState.rog}
                          min={-50}
                          max={50}
                          step={1}
                        />
                      </div>
                      <div className="field">
                        <label>Extend Gap Score</label>
                        <input
                          type="number"
                          onChange={(v) =>
                            setAdvancedSettingsState({
                              ...advancedSettingsState,
                              reg: parseFloat(v.target.value),
                            })
                          }
                          value={advancedSettingsState.reg}
                          min={-50}
                          max={50}
                          step={1}
                        />
                      </div>
                    </div>
                  </div>

                  <div className="group">
                    <div className="field">
                      <label className="header">
                        Alignment Output Directory
                      </label>
                      <input
                        type="text"
                        value={appState.alignment_output_path}
                      />
                      <button
                        type="button"
                        onClick={() =>
                          window.pywebview.api.select_alignment_output_path()
                        }
                      >
                        Select Folder...
                      </button>
                    </div>
                  </div>
                </details>

                <div className="field">
                  <label className="header">
                    <input
                      type="checkbox"
                      checked={runnerSettings.neighbor_joining}
                      defaultChecked={runnerSettings.neighbor_joining}
                      onChange={toggleNeighborJoining}
                    />
                    <span>Cluster Sequences by Neighbor Joining</span>
                  </label>
                </div>

                <div className="field runner-settings">
                  <label className="header">Alignment Type</label>
                  {alignmentTypes.map((value) => (
                    <div key={value}>
                      <input
                        type="radio"
                        id={value}
                        name="alignment-type"
                        value={value}
                        checked={runnerSettings.alignment_type === value}
                        onChange={() =>
                          handleChangeAlignmentType(
                            value as (typeof alignmentTypes)[keyof typeof alignmentTypes],
                          )
                        }
                      />
                      <label htmlFor={value}>{value}</label>
                    </div>
                  ))}
                </div>

                <div className="field runner-settings">
                  <label className="header">Performance</label>
                  {performanceProfiles.map((value) => (
                    <div key={value}>
                      <input
                        type="radio"
                        id={value}
                        name="performance-profile"
                        value={value}
                        checked={runnerSettings.performance_profile === value}
                        onChange={() =>
                          handleChangePerformanceProfile(
                            value as (typeof performanceProfiles)[keyof typeof performanceProfiles],
                          )
                        }
                      />
                      <label htmlFor={value}>{value}</label>
                    </div>
                  ))}
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
