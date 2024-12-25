import React from "react";
import { Button, Label, ProgressBar } from "react-aria-components";
import type { DocState, SetDocState, UpdateDocState } from "../appState";
import { services } from "../services";
import { LoadingAnimation } from "./LoadingAnimation";

export const Loader = ({
  docState: { id: docId, stage, progress, estimated_time },
  setDocState,
  mainMenu,
}: {
  docState: DocState;
  setDocState: SetDocState;
  updateDocState: UpdateDocState;
  mainMenu: React.ReactNode;
}) => {
  const [canceling, setCanceling] = React.useState(false);
  const [estimatedDisplay, setEstimatedDisplay] = React.useState("");
  const startTime = React.useRef(Date.now());
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

  React.useEffect(() => {
    const handler = () => {
      services.getDocument(docId).then((data) => {
        setDocState((prev) => ({ ...prev, ...data }));
      });
    };
    const id = setInterval(handler, 500);
    return () => clearInterval(id);
  }, [docId, setDocState]);

  return (
    <div className="app-wrapper with-header loader">
      <div className="app-header loader">
        <div className="left">{mainMenu}</div>
        <div className="right" />
      </div>
      <div className="app-main centered loader">
        <div className="loader-wrapper">
          <ProgressBar value={progress}>
            {({ percentage, valueText }) => (
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
          {canceling ? "Canceling..." : "Cancel Run"}
        </Button>
      </div>
    </div>
  );
};
