import {
  type AppState,
  clientStateKey,
  clientStateSchema,
  initialAppState,
} from "./appState";
import { partialSafeParse } from "./zodUtils";

export const restoreClientState = (baseClientState: AppState["client"]) => {
  try {
    const state = localStorage.getItem(clientStateKey);
    if (state) {
      const parsedState = JSON.parse(state);
      const parsedClient = partialSafeParse(
        clientStateSchema,
        parsedState.client,
      );
      const merged: AppState["client"] = {
        ...baseClientState,
        ...parsedClient.validData,
      };
      return merged;
    }
  } catch (e) {
    return initialAppState.client;
  }
  return initialAppState.client;
};
