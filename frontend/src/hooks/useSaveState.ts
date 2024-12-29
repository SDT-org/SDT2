import React from "react";
import { type AppState, clientStateKey } from "../appState";

export const useSaveState = (initialized: boolean, appState: AppState) =>
  React.useEffect(() => {
    if (!initialized) {
      return;
    }
    localStorage.setItem(clientStateKey, JSON.stringify(appState));
    if (appState.debug) {
      window.APP_STATE = appState;
    }
  }, [initialized, appState]);
