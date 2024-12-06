import React from "react";
import { type AppState, clientStateKey } from "../appState";

export const useSaveState = (initialized: boolean, appState: AppState) =>
  React.useEffect(() => {
    if (!initialized) {
      return;
    }
    localStorage.setItem(
      clientStateKey,
      JSON.stringify({
        ...appState.client,
        lastDataFilePath:
          typeof appState.filename === "string"
            ? appState.filename
            : appState.filename[0],
      }),
    );
    window.APP_STATE = appState;
  }, [initialized, appState]);
