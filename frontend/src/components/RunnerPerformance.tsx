import React from "react";
import {
  Label,
  Meter,
  Slider,
  SliderOutput,
  SliderThumb,
  SliderTrack,
} from "react-aria-components";
import { TbAlertTriangleFilled } from "react-icons/tb";
import type {
  AppState,
  DocState,
  SetDocState,
  UpdateDocState,
} from "../appState";
import { formatBytes } from "../helpers";

export const RunnerPerformance = ({
  appState,
  docState,
  setDocState,
  updateDocState,
  initialized,
}: {
  appState: AppState;
  docState: DocState;
  setDocState: SetDocState;
  updateDocState: UpdateDocState;
  initialized: React.MutableRefObject<boolean>;
}) => {
  React.useEffect(() => {
    const id = setInterval(() => {
      if (!docState.filename) {
        return;
      }

      window.pywebview.api.system
        .get_available_memory()
        .then((available_memory) =>
          setDocState((previous) => ({
            ...previous,
            compute_stats: {
              ...(previous.compute_stats || {
                recommended_cores: 1,
                required_memory: 1,
              }),
              available_memory,
            },
          })),
        );
    }, 3000);

    return () => clearInterval(id);
  }, [docState.filename, setDocState]);

  React.useEffect(() => {
    if (
      !docState.compute_stats ||
      initialized.current ||
      docState.compute_cores
    ) {
      return;
    }
    updateDocState({
      compute_cores: docState.compute_stats.recommended_cores,
    });
    initialized.current = true;
  }, [
    docState.compute_stats,
    docState.compute_cores,
    updateDocState,
    initialized,
  ]);

  const estimatedMemory =
    (docState.compute_stats?.required_memory || 1) *
    (docState.compute_cores || 1);

  const estimatedMemoryValue = docState.compute_stats
    ? ((docState.compute_stats.required_memory *
        (docState.compute_cores || 1)) /
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
      : impactName(
          ((docState.compute_cores || 1) / appState.platform.cores) * 100,
        );

  const showWarning =
    docState.compute_stats &&
    (docState.compute_stats?.recommended_cores === 0 ||
      estimatedMemory > appState.platform.memory ||
      (docState.compute_cores || 1) >
        docState.compute_stats?.recommended_cores);

  return (
    <>
      <div
        className="field performance"
        data-hidden={docState.analysisMethod === "lzani"}
      >
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
                onChange={(value) => updateDocState({ compute_cores: value })}
                minValue={1}
                maxValue={appState.platform.cores}
                value={
                  docState.compute_cores ||
                  docState.compute_stats.recommended_cores
                }
              >
                <SliderOutput data-impact={coresImpact}>
                  {({ state }) => (
                    <>
                      {docState.compute_stats &&
                      (docState.compute_stats.recommended_cores === 0 ||
                        (docState.compute_cores || 1) >
                          docState.compute_stats.recommended_cores) ? (
                        <TbAlertTriangleFilled />
                      ) : null}
                      {state.getThumbValueLabel(0)} / {appState.platform.cores}
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
                              docState.compute_stats?.available_memory || 0,
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
                  {docState.compute_stats && estimatedMemoryValue > 100 ? (
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
      {showWarning ? (
        <div className="compute-forecast">
          <b>Warning:</b> System instability may occur when exceeding
          recommended cores or available memory. When exceeding available
          memory, it is recommended to limit cores used.
        </div>
      ) : null}
    </>
  );
};
