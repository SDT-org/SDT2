import React from "react";
import { AppState } from "../appState";
import { Button, Label, ProgressBar } from "react-aria-components";

export const Loader = ({ appState }: { appState: AppState }) => {
  const [canceling, setCanceling] = React.useState(false);

  return (
    <div className="app-wrapper">
      <div className="app-main centered loader">
        <ProgressBar value={appState.progress}>
          {({ percentage, valueText }) => (
            <>
              <Label>{appState.stage ? `${appState.stage}...` : null}</Label>
              <span className="value">{valueText}</span>
              <div className="bar">
                <div className="fill" style={{ width: percentage + "%" }} />
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
          {canceling ? "Canceling..." : "Cancel Run..."}
        </Button>
      </div>
    </div>
  );
};
