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
  type AppState,
  type SetAppState,
  type DocState,
  type SetDocState,
  type UpdateDocState,
} from "../appState";
import { lzaniPresets } from "../config/lzaniSettings";
import {
  getRecommendedMatrix,
  scoringMatrices,
} from "../config/scoringMatrices";
import { reorderMethods } from "../constants";
import { splitFilePath } from "../helpers";
import useOpenFileDialog from "../hooks/useOpenFileDialog";
import { useRecentFiles } from "../hooks/useRecentFiles";
import { useStartRun } from "../hooks/useStartRun";
import messages from "../messages";
import { NumberInput } from "./NumberInput";
import { RunnerPerformance } from "./RunnerPerformance";
import { Select, SelectItem } from "./Select";
import { Switch } from "./Switch";

export type ParasailRunSettings = {
  scoring_matrix?: string;
  open_penalty?: number;
  extend_penalty?: number;
};

export type LzaniRunSettings = {
  aw?: number;
  am?: number;
  mal?: number;
  msl?: number;
  mrd?: number;
  mqd?: number;
  reg?: number;
  ar?: number;
};

export type RunSettings = Pick<DocState, "compute_cores" | "analysisMethod"> & {
  doc_id: string;
  cluster_method: DocState["cluster_method"] | "None";
  lzani_score_type?: string;
  export_alignments: boolean;
  alignment_export_path: string;
} & Partial<ParasailRunSettings> &
  Partial<LzaniRunSettings>;

const DefaultForm = ({
  appState,
  setAppState,
  docState,
}: {
  appState: AppState;
  setAppState: SetAppState;
  docState: DocState;
}) => {
  const openFileDialog = useOpenFileDialog(appState, setAppState);
  const openRecentFile = useRecentFiles(appState, setAppState);

  return (
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
  );
};

const ParasailSettings = ({
  docState,
  updateDocState,
  setDocState,
}: {
  docState: DocState;
  updateDocState: UpdateDocState;
  setDocState: SetDocState;
}) => {
  const isAminoAcid = docState.result_metadata?.is_aa;
  const [showAdvanced, setShowAdvanced] = React.useState<boolean>(false);

  React.useEffect(() => {
    // Always enable override for Parasail to match LZ-ANI behavior
    if (!docState.overrideParasail) {
      updateDocState({ overrideParasail: true });
    }
  }, [docState.overrideParasail, updateDocState]);

  return (
    <div className="field">
      <div className="header">Parasail Settings</div>
      <div
        className="setting-group"
        style={{
          backgroundColor: "#f8f8f8",
          padding: "0.75rem",
          borderRadius: "6px",
          border: "1px solid #e0e0e0",
          marginTop: "0.5rem",
        }}
      >
        <div
          style={{
            display: "flex",
            gap: "0.75rem",
            alignItems: "flex-end",
            justifyContent: "space-between",
          }}
        >
          <div style={{ flex: "0 0 auto", minWidth: "200px" }}>
            <Label
              htmlFor="scoring-matrix"
              style={{
                fontSize: "0.75rem",
                fontWeight: "500",
                color: "#666",
                marginBottom: "0.25rem",
                display: "block",
              }}
            >
              Scoring Matrix
            </Label>
            <Select
              id="scoring-matrix"
              wide
              data-compact
              selectedKey={
                docState.parasail_settings?.scoring_matrix ||
                getRecommendedMatrix(!!isAminoAcid)
              }
              onSelectionChange={(value) => {
                setDocState((previous) => ({
                  ...previous,
                  parasail_settings: {
                    ...previous.parasail_settings,
                    scoring_matrix: value as string,
                  },
                }));
              }}
              items={scoringMatrices.map((matrix) => ({
                ...matrix,
                recommended: matrix.id === getRecommendedMatrix(!!isAminoAcid),
              }))}
            >
              {(item) => (
                <SelectItem textValue={item.name}>
                  <Text slot="label">{item.name}</Text>
                  {item.recommended ? (
                    <Text slot="description">
                      <small>
                        <em>
                          Recommended for
                          {!isAminoAcid
                            ? " Nucletide datasets"
                            : " Amino-acid datasets"}
                        </em>
                      </small>
                    </Text>
                  ) : null}
                </SelectItem>
              )}
            </Select>
          </div>
          <Button
            type="button"
            className="react-aria-Button"
            onPress={() => setShowAdvanced(!showAdvanced)}
            style={{
              padding: "0.375rem 0.875rem",
              fontSize: "0.8125rem",
              fontWeight: "500",
              minHeight: "2rem",
              backgroundColor: showAdvanced ? "#007acc" : "#f0f0f0",
              color: showAdvanced ? "white" : "#333",
              border: showAdvanced ? "1px solid #007acc" : "1px solid #d0d0d0",
              borderRadius: "4px",
              transition: "all 0.2s ease",
              cursor: "pointer",
            }}
          >
            {showAdvanced ? "← Hide Advanced" : "Advanced →"}
          </Button>
        </div>
      </div>
      {showAdvanced ? (
        <div className="setting-group form col-2" style={{ marginTop: "1rem" }}>
          <NumberInput
            id="open-penalty"
            label="Gap open penalty"
            value={
              docState.parasail_settings?.open_penalty || (isAminoAcid ? 10 : 8)
            }
            onChange={(value) => {
              setDocState((previous) => ({
                ...previous,
                parasail_settings: {
                  ...previous.parasail_settings,
                  open_penalty: value,
                },
              }));
            }}
            min={1}
            max={99}
          />
          <NumberInput
            id="extend-penalty"
            label="Gap extend penalty"
            value={docState.parasail_settings?.extend_penalty || 1}
            onChange={(value) => {
              setDocState((previous) => ({
                ...previous,
                parasail_settings: {
                  ...previous.parasail_settings,
                  extend_penalty: value,
                },
              }));
            }}
            min={1}
            max={99}
          />
        </div>
      ) : null}
    </div>
  );
};

const LzaniSettings = ({
  docState,
  updateDocState,
  setDocState,
}: {
  docState: DocState;
  updateDocState: UpdateDocState;
  setDocState: SetDocState;
}) => {
  const [selectedPreset, setSelectedPreset] =
    React.useState<string>("balanced");
  const [showAdvanced, setShowAdvanced] = React.useState<boolean>(false);

  React.useEffect(() => {
    // Always enable override for LZ-ANI
    if (!docState.overrideLzani) {
      updateDocState({ overrideLzani: true });
    }
  }, [docState.overrideLzani, updateDocState]);

  const handlePresetChange = (presetId: string) => {
    setSelectedPreset(presetId);
    const preset = lzaniPresets.find((p) => p.id === presetId);
    if (preset) {
      setDocState((previous) => ({
        ...previous,
        lzani_settings: preset.settings,
      }));
    }
  };

  return (
    <div className="field">
      <div className="header">LZ-ANI Settings</div>
      <div
        className="setting-group"
        style={{
          backgroundColor: "#f8f8f8",
          padding: "0.75rem",
          borderRadius: "6px",
          border: "1px solid #e0e0e0",
          marginTop: "0.5rem",
        }}
      >
        <div
          style={{
            display: "flex",
            gap: "0.75rem",
            alignItems: "flex-end",
            flexWrap: "wrap",
            justifyContent: "space-between",
          }}
        >
          <div
            style={{ display: "flex", gap: "0.75rem", alignItems: "flex-end" }}
          >
            <div style={{ flex: "0 0 auto", minWidth: "140px" }}>
              <Label
                style={{
                  fontSize: "0.75rem",
                  fontWeight: "500",
                  color: "#666",
                  marginBottom: "0.25rem",
                  display: "block",
                }}
              >
                Settings
              </Label>
              <Select
                id="lzani-preset"
                wide
                data-compact
                isDisabled={showAdvanced}
                selectedKey={selectedPreset}
                onSelectionChange={(value) => {
                  handlePresetChange(value as string);
                }}
                items={lzaniPresets}
                aria-label="LZ-ANI Preset"
              >
                {(item) => (
                  <SelectItem textValue={`${item.name} - ${item.description}`}>
                    <Text slot="label">{item.name}</Text>
                    <Text slot="description">
                      <small style={{ opacity: 0.7 }}>{item.description}</small>
                    </Text>
                  </SelectItem>
                )}
              </Select>
            </div>
            <div
              style={{ flex: "0 0 auto", minWidth: "100px", maxWidth: "150px" }}
            >
              <Label
                style={{
                  fontSize: "0.75rem",
                  fontWeight: "500",
                  color: "#666",
                  marginBottom: "0.25rem",
                  display: "block",
                }}
              >
                Score Type
              </Label>
              <Select
                id="lzani-score-type"
                wide
                data-compact
                selectedKey={docState.lzaniScoreType}
                onSelectionChange={(value) => {
                  updateDocState({
                    lzaniScoreType: value as DocState["lzaniScoreType"],
                  });
                }}
                items={[
                  {
                    id: "ani",
                    name: "ANI",
                  },
                  {
                    id: "tani",
                    name: "TANI",
                  },
                ]}
                aria-label="Score Type"
              >
                {(item) => (
                  <SelectItem textValue={item.name}>{item.name}</SelectItem>
                )}
              </Select>
            </div>
          </div>
          <Button
            type="button"
            className="react-aria-Button"
            onPress={() => setShowAdvanced(!showAdvanced)}
            style={{
              padding: "0.375rem 0.875rem",
              fontSize: "0.8125rem",
              fontWeight: "500",
              minHeight: "2rem",
              backgroundColor: showAdvanced ? "#007acc" : "#f0f0f0",
              color: showAdvanced ? "white" : "#333",
              border: showAdvanced ? "1px solid #007acc" : "1px solid #d0d0d0",
              borderRadius: "4px",
              transition: "all 0.2s ease",
              cursor: "pointer",
            }}
          >
            {showAdvanced ? "← Presets" : "Advanced →"}
          </Button>
        </div>
      </div>
      {showAdvanced ? (
        <div className="setting-group form col-3" style={{ marginTop: "1rem" }}>
          <NumberInput
            id="lzani-aw"
            label="Anchor width"
            value={docState.lzani_settings?.aw || 15}
            onChange={(value) => {
              setDocState((previous) => ({
                ...previous,
                lzani_settings: {
                  ...previous.lzani_settings,
                  aw: value,
                },
              }));
            }}
            min={1}
            max={30}
          />
          <NumberInput
            id="lzani-am"
            label="Anchor mismatch"
            value={docState.lzani_settings?.am || 7}
            onChange={(value) => {
              setDocState((previous) => ({
                ...previous,
                lzani_settings: {
                  ...previous.lzani_settings,
                  am: value,
                },
              }));
            }}
            min={0}
            max={15}
          />
          <NumberInput
            id="lzani-mal"
            label="Min anchor length"
            value={docState.lzani_settings?.mal || 11}
            onChange={(value) => {
              setDocState((previous) => ({
                ...previous,
                lzani_settings: {
                  ...previous.lzani_settings,
                  mal: value,
                },
              }));
            }}
            min={1}
            max={25}
          />
          <NumberInput
            id="lzani-msl"
            label="Min seed length"
            value={docState.lzani_settings?.msl || 7}
            onChange={(value) => {
              setDocState((previous) => ({
                ...previous,
                lzani_settings: {
                  ...previous.lzani_settings,
                  msl: value,
                },
              }));
            }}
            min={1}
            max={15}
          />
          <NumberInput
            id="lzani-mrd"
            label="Max ref. distance"
            value={docState.lzani_settings?.mrd || 40}
            onChange={(value) => {
              setDocState((previous) => ({
                ...previous,
                lzani_settings: {
                  ...previous.lzani_settings,
                  mrd: value,
                },
              }));
            }}
            min={1}
            max={120}
          />
          <NumberInput
            id="lzani-mqd"
            label="Max query distance"
            value={docState.lzani_settings?.mqd || 40}
            onChange={(value) => {
              setDocState((previous) => ({
                ...previous,
                lzani_settings: {
                  ...previous.lzani_settings,
                  mqd: value,
                },
              }));
            }}
            min={1}
            max={120}
          />
          <NumberInput
            id="lzani-reg"
            label="Min region length"
            value={docState.lzani_settings?.reg || 35}
            onChange={(value) => {
              setDocState((previous) => ({
                ...previous,
                lzani_settings: {
                  ...previous.lzani_settings,
                  reg: value,
                },
              }));
            }}
            min={1}
            max={80}
          />
          <NumberInput
            id="lzani-ar"
            label="Min end length"
            value={docState.lzani_settings?.ar || 3}
            onChange={(value) => {
              setDocState((previous) => ({
                ...previous,
                lzani_settings: {
                  ...previous.lzani_settings,
                  ar: value,
                },
              }));
            }}
            min={0.1}
            max={10}
            step={0.1}
          />
        </div>
      ) : null}
    </div>
  );
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
  const fileName =
    docState.filename?.length && docState.filename
      ? docState.filename.split("/").pop()
      : undefined;
  const isFastaType = docState.filetype === "text/fasta";
  const openFileDialog = useOpenFileDialog(appState, setAppState);

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
                <Label data-header>
                  Analysis Method
                  {docState.result_metadata?.is_aa ? (
                    <span className="text-deemphasis is-aa-detected">
                      Amino acid detected
                    </span>
                  ) : null}
                </Label>
                <div className="cards">
                  <Radio value="parasail">
                    <div className="col-2 analysis-method-body">
                      <div>
                        Parasail
                        <p className="text-deemphasis">
                          Best for small datasets. Supports amino acid
                          calculations.
                        </p>
                      </div>
                    </div>
                  </Radio>
                  <Radio
                    value="lzani"
                    isDisabled={docState.result_metadata?.is_aa || false}
                  >
                    <div className="col-2 analysis-method-body">
                      <div>
                        LZ-ANI
                        <p className="text-deemphasis">
                          Best for large datasets. Only supports nucleotide ANI
                          calculations.
                        </p>
                      </div>
                    </div>
                  </Radio>
                </div>
              </RadioGroup>
            </div>

            {docState.analysisMethod === "parasail" ? (
              <ParasailSettings
                docState={docState}
                updateDocState={updateDocState}
                setDocState={setDocState}
              />
            ) : null}

            {docState.analysisMethod === "lzani" ? (
              <LzaniSettings
                docState={docState}
                updateDocState={updateDocState}
                setDocState={setDocState}
              />
            ) : null}

            <div className="field">
              <Switch
                isSelected={docState.exportAlignments}
                onChange={(value) => {
                  updateDocState({ exportAlignments: value });
                }}
              >
                Export alignment files
              </Switch>
              {docState.exportAlignments ? (
                <div
                  className="setting input-with-button"
                  style={{ marginTop: "0.75rem", width: "100%" }}
                >
                  <Input
                    id="alignment-export-path"
                    type="text"
                    readOnly
                    value={docState.alignmentExportPath}
                    placeholder="Select output directory..."
                  />
                  <Button
                    type="button"
                    onPress={async () => {
                      const path =
                        await window.pywebview.api.files.select_path_dialog(
                          docState.alignmentExportPath || "",
                        );
                      if (path) {
                        updateDocState({ alignmentExportPath: path });
                      }
                    }}
                  >
                    Browse...
                  </Button>
                </div>
              ) : null}
            </div>

            <div className="field clustering inline-setting">
              <Switch
                isSelected={docState.enableClustering}
                onChange={(value) => {
                  updateDocState({ enableClustering: value });
                }}
              >
                Reorder data by linkage clustering method
              </Switch>
              <div className="setting clustering-method">
                <Select
                  isDisabled={!docState.enableClustering}
                  data-flat
                  selectedKey={docState.cluster_method}
                  onSelectionChange={(value) => {
                    updateDocState({
                      cluster_method: value as typeof docState.cluster_method,
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
          <DefaultForm
            appState={appState}
            setAppState={setAppState}
            docState={docState}
          />
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
      className="app-panel app-panel-full-width app-panel-full-height"
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
