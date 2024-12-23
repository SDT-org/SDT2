import {
  type AppState,
  // clientStateKey,
  // clientStateSchema,
  initialAppState,
} from "./appState";
// import { partialSafeParse } from "./zodUtils";

export const restoreClientState = (baseClientState: AppState) => {
  try {
    return baseClientState || initialAppState;
    // const state = localStorage.getItem(clientStateKey);
    // if (state) {
    // const parsedState = JSON.parse(state);
    //   parsedState.error = null;
    //   parsedState.errorInfo = null;
    //   const parsedClient = partialSafeParse(clientStateSchema, parsedState);
    //   const validData = parsedClient.validData;
    //   const merged: AppState = {
    //     ...baseClientState,
    //     ...parsedClient.validData,
    //     distribution: {
    //       ...baseClientState.distribution,
    //       ...validData.distribution,
    //       histogram: {
    //         ...baseClientState.distribution.histogram,
    //         ...validData.distribution?.histogram,
    //       },
    //       violin: {
    //         ...baseClientState.distribution.violin,
    //         ...validData.distribution?.violin,
    //       },
    //       raincloud: {
    //         ...baseClientState.distribution.raincloud,
    //         ...validData.distribution?.raincloud,
    //       },
    //     },
    //     heatmap: {
    //       ...baseClientState.heatmap,
    //       ...validData.heatmap,
    //     },
    //   };
    //   return merged;
    // }
  } catch (e) {
    return initialAppState;
  }
  // return initialAppState;s
};
