import React from "react";
import { Button, Label, ProgressBar } from "react-aria-components";
import type { AppState } from "../appState";
import { formatBytes } from "../helpers";
import { LoadingAnimation } from "./LoadingAnimation";

export const Loader = ({
  appState: { stage, progress, estimated_time, debug, compute_stats },
  mainMenu,
}: {
  appState: AppState;
  mainMenu: React.ReactNode;
}) => {
  const [canceling, setCanceling] = React.useState(false);
  const [estimatedDisplay, setEstimatedDisplay] = React.useState("");
  const startTime = React.useRef(Date.now());
  const updatedTime = React.useRef(Date.now());
  const [processInfo, setProcessInfo] = React.useState<
    [number, number, number, string][]
  >([]);

  React.useEffect(() => {
    if (!estimated_time || Date.now() - updatedTime.current < 3000) {
      return;
    }

    if (estimated_time < 60) {
      setEstimatedDisplay("Estimated: <1 minute");
      updatedTime.current = Date.now();
      return;
    }

    const estimatedMinute = Math.round(estimated_time / 60);
    setEstimatedDisplay(
      `Estimated: ${estimatedMinute} minute${estimatedMinute > 1 ? "s" : ""}`,
    );
    updatedTime.current = Date.now();
  }, [estimated_time]);

  React.useEffect(() => {
    if (!debug) {
      return;
    }

    const id = setInterval(
      () =>
        window.pywebview.api.processes_info().then((data) => {
          const parsedData = JSON.parse(data);
          if (!parsedData) {
            return;
          }

          setProcessInfo(parsedData.sort());
        }),
      500,
    );
    return () => clearInterval(id);
  }, [debug]);

  return (
    <div className="app-wrapper with-header loader">
      <div className="app-header loader">
        <div className="left">{mainMenu}</div>
        {Date.now() - startTime.current > 3000 && estimated_time ? (
          <LoadingAnimation className="loader-loading-animation" />
        ) : null}
        <div className="right" />
      </div>
      <div className="app-main centered loader">
        <div className="form-wrapper loader-wrapper">
          <ProgressBar value={progress}>
            {({ percentage, valueText }) => (
              <>
                <Label>{stage ? <>{stage}...</> : null}</Label>
                <span className="value">
                  {stage === "Analyzing" && valueText}
                </span>
                <div className="bar">
                  <div
                    className="fill"
                    style={{
                      width: `${percentage}%`,
                    }}
                    data-animation={estimated_time && estimated_time > 10}
                  />
                </div>
                <div className="estimate text-secondary">
                  {(stage === "Analyzing" && estimatedDisplay) || <>&nbsp;</>}
                </div>
              </>
            )}
          </ProgressBar>
          {debug ? (
            <details
              style={{
                position: "absolute",
                bottom: "1.6rem",
                left: "1.6rem",
                fontSize: "1rem",
              }}
            >
              <summary>🔬</summary>
              <pre>
                Required: {formatBytes(compute_stats?.required_memory || 0)}
                <br />
                {processInfo.map(
                  (i) =>
                    `[${i[0]}] ${i[1].toString().padStart(5, " ")} ${formatBytes(i[2]).padStart(5, " ")} ${i[3]}\n`,
                )}
              </pre>
            </details>
          ) : null}
        </div>
        <Button
          className="react-aria-Button cancel-run"
          onPress={() => {
            if (
              window.confirm(
                "Are you sure you want to cancel? All analysis progress will be lost.",
              )
            ) {
              setCanceling(true);
              window.pywebview.api.cancel_run();
            }
          }}
          isDisabled={canceling}
        >
          {canceling ? "Canceling..." : "Cancel Run"}
        </Button>
      </div>
    </div>
  );
};
