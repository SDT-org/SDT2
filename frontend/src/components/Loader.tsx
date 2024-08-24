import React from "react";
import { AppState } from "../appState";
import { Button, Label, ProgressBar } from "react-aria-components";

export const Loader = ({
  appState: { stage, progress, estimated_time },
}: {
  appState: AppState;
}) => {
  const [canceling, setCanceling] = React.useState(false);
  const [estimatedDisplay, setEstimatedDisplay] = React.useState("");
  const updatedTime = React.useRef(Date.now());

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

  return (
    <div className="app-wrapper">
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
