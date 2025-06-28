import React from "react";
import { Button, Label, ProgressBar, TabPanel } from "react-aria-components";
import type { DocState, SetDocState, UpdateDocState } from "../appState";
import { LoadingAnimation } from "./LoadingAnimation";

export const Loader = ({
  docState: { id: docId, estimated_time },
  updateDocState,
}: {
  docState: DocState;
  setDocState: SetDocState;
  updateDocState: UpdateDocState;
  tabView: "tabs" | "select";
}) => {
  const [canceling, setCanceling] = React.useState(false);
  const [estimatedDisplay, setEstimatedDisplay] = React.useState("");
  const startTime = React.useRef(Date.now());
  const updatedTime = React.useRef(Date.now());
  const [stage, setStage] = React.useState("Initializing");
  const [progress, setProgress] = React.useState<number | undefined>();

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
    const handler = () => {
      window.pywebview.api.get_workflow_run_status(docId).then((data) => {
        setStage(data.stage);
        setProgress(data.progress);

        // TODO: start passing down an active workflow state from app tabs so this isn't needed
        updateDocState({
          ...(progress !== undefined && { progress: data.progress }),
          stage: data.stage,
        });
      });
    };
    const id = setInterval(handler, 120);
    return () => clearInterval(id);
  }, [docId, updateDocState, progress]);

  return (
    <TabPanel id={docId} key={docId} className={"app-panel full-width"}>
      <div className="app-main centered loader">
        <div className="loader-wrapper">
          <ProgressBar
            isIndeterminate={Number.isNaN(progress)}
            value={progress || 0}
          >
            {({ percentage, valueText, isIndeterminate }) => (
              <>
                <Label>
                  {stage ? (
                    <>
                      {stage}
                      {Date.now() - startTime.current > 3000 ? (
                        <LoadingAnimation className="loader-loading-animation" />
                      ) : null}
                    </>
                  ) : null}
                </Label>
                <span className="value">
                  {stage === "Analyzing" && valueText}
                </span>
                <div className="bar">
                  <div
                    className="fill"
                    style={{
                      width: `${percentage}%`,
                    }}
                    data-animation
                    data-indeterminate={isIndeterminate}
                  />
                </div>
                <div className="estimate text-secondary">
                  {(stage === "Analyzing" && estimatedDisplay) || <>&nbsp;</>}
                </div>
              </>
            )}
          </ProgressBar>
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
              window.pywebview.api.cancel_run(docId, "preserve");
            }
          }}
          isDisabled={canceling}
        >
          {canceling ? "Canceling..." : "Cancel"}
        </Button>
      </div>
    </TabPanel>
  );
};
