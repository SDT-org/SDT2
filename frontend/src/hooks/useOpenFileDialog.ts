import React from "react";
// import useAppState from "../appState";
import type { AppState, SetAppState } from "../appState";
import { services } from "../services";

const useOpenFileDialog = (appState: AppState, setAppState: SetAppState) => {
  // const useOpenFileDialog = () => {
  // const { appState, setAppState } = useAppState();

  return React.useCallback(() => {
    services.openFileDialog(appState.client.lastDataFilePath).then((data) => {
      const [success, value] = data;
      if (!success) {
        return;
      }
      setAppState((prev) => ({
        ...prev,
        client: { ...prev.client, lastDataFilePath: value },
      }));
    });
  }, [appState.client.lastDataFilePath, setAppState]);
};

export default useOpenFileDialog;
