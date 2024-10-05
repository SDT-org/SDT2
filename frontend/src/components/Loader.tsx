import React from "react";
import { AppState } from "../appState";
import { Button, Label, ProgressBar } from "react-aria-components";
import { formatBytes } from "../helpers";

export const Loader = ({
  appState: { stage, progress, estimated_time, debug, compute_stats },
  mainMenu,
}: {
  appState: AppState;
  mainMenu: React.ReactNode;
}) => {
  const [canceling, setCanceling] = React.useState(false);
  const [estimatedDisplay, setEstimatedDisplay] = React.useState("");
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
  }, []);

  return (
    <div className="app-wrapper with-header loader">
      <div className="app-header loader">
        <div>{mainMenu}</div>
      </div>
      <div className="app-main centered loader">
        <ProgressBar value={progress}>
          {({ percentage, valueText }) => (
            <>
              <Label>{stage ? `${stage}...` : null}</Label>
              <span className="value">
                {stage === "Analyzing" && valueText}
              </span>
              <div className="bar">
                <div
                  className="fill"
                  style={{ width: percentage + "%" }}
                  data-animation={estimated_time && estimated_time > 10}
                />
              </div>
              <div className="estimate">
                {(stage === "Analyzing" && estimatedDisplay) || <>&nbsp;</>}
              </div>
            </>
          )}
        </ProgressBar>
        {debug ? (
          <pre>
            Required: {formatBytes(compute_stats?.required_memory || 0)}
            <br />
            {processInfo.map(
              (i) =>
                `[${i[0]}] ${i[1].toString().padStart(5, " ")} ${formatBytes(i[2]).padStart(5, " ")} ${i[3]}\n`,
            )}
          </pre>
        ) : null}
        <Button
          className="cancel-run"
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
