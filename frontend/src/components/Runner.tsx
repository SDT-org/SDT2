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
  type ScoringMatrix,
} from "../appState";
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

export type RunSettings = Pick<DocState, "compute_cores" | "analysisMethod"> & {
  doc_id: string;
  cluster_method: DocState["cluster_method"] | "None";
};

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

  return (
    <div className="field">
      <Switch
        isSelected={docState.overrideParasail}
        onChange={(value) => {
          updateDocState({ overrideParasail: value });
        }}
      >
        Override default Parasail settings
      </Switch>
      {docState.overrideParasail ? (
        <div className="setting-group form col-3">
          <div className="field">
            <Label htmlFor="scoring-matrix">Scoring matrix</Label>
            <Select
              id="scoring-matrix"
              wide
              data-compact
              selectedKey={
                docState.parasail_settings?.scoring_matrix ||
                (isAminoAcid ? "blosum62" : "pam250")
              }
              onSelectionChange={(value) => {
                setDocState((previous) => ({
                  ...previous,
                  parasail_settings: {
                    ...previous.parasail_settings,
                    scoring_matrix: value as ScoringMatrix["id"],
                  },
                }));
              }}
              items={[
                {
                  id: "blosum62",
                  name: "BLOSUM62",
                  description: "Default scoring matrix for amino acids",
                },
                {
                  id: "pam250",
                  name: "PAM250",
                  description: "TODO",
                },
              ]}
            >
              {(item) => (
                <SelectItem textValue={item.name}>
                  <Text slot="label">{item.name}</Text>
                  <Text slot="description">{item.description}</Text>
                </SelectItem>
              )}
            </Select>
          </div>
          <NumberInput
            id="open-penalty"
            label="Open Penalty"
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
            label="Extend Penalty"
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
                <Label data-header>Analysis Method</Label>
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
                  <Radio value="lzani">
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
