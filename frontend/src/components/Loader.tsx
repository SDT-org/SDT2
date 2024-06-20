import React from "react";
import { AppState } from "../appState";
import { Label, ProgressBar } from "react-aria-components";

export const Loader = ({ appState }: { appState: AppState }) => (
  <div className="app-wrapper">
    <div className="app-main centered loader">
      <ProgressBar value={appState.progress}>
        {({ percentage, valueText }) => (
          <>
            <Label>{appState.stage}...</Label>
            <span className="value">{valueText}</span>
            <div className="bar">
              <div className="fill" style={{ width: percentage + "%" }} />
            </div>
          </>
        )}
      </ProgressBar>
    </div>
  </div>
);
