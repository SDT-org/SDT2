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
      parsedState.error = null;
      parsedState.errorInfo = null;
      const parsedClient = partialSafeParse(clientStateSchema, parsedState);
      const validData = parsedClient.validData;
      console.debug(parsedClient.error);
      const merged: AppState["client"] = {
        ...baseClientState,
        ...parsedClient.validData,
        distribution: {
          ...baseClientState.distribution,
          ...validData.distribution,
          histogram: {
            ...baseClientState.distribution.histogram,
            ...validData.distribution?.histogram,
          },
          violin: {
            ...baseClientState.distribution.violin,
            ...validData.distribution?.violin,
          },
          raincloud: {
            ...baseClientState.distribution.raincloud,
            ...validData.distribution?.raincloud,
          },
        },
        heatmap: {
          ...baseClientState.heatmap,
          ...validData.heatmap,
        },
      };
      return merged;
    }
  } catch (e) {
    console.error(e);
    return initialAppState.client;
  }
  return initialAppState.client;
};
