import { type DocState, docStateSchema, initialDocState } from "./appState";
import { partialSafeParse } from "./zodUtils";

export const parseDocState = (state: DocState) => {
  try {
    const parsedState = partialSafeParse(docStateSchema, state);
    const validData = parsedState.validData;
    const merged: DocState = {
      ...initialDocState,
      ...validData,
      ...("compute_stats" in state
        ? { compute_stats: state.compute_stats }
        : {}),
      distribution: {
        ...initialDocState.distribution,
        ...validData.distribution,
        histogram: {
          ...initialDocState.distribution.histogram,
          ...validData.distribution?.histogram,
        },
        violin: {
          ...initialDocState.distribution.violin,
          ...validData.distribution?.violin,
        },
        raincloud: {
          ...initialDocState.distribution.raincloud,
          ...validData.distribution?.raincloud,
        },
      },
      heatmap: {
        ...initialDocState.heatmap,
        ...validData.heatmap,
      },
      parsed: true,
    };

    if (parsedState.error) {
      console.warn(parsedState.error);
    }

    return merged;
  } catch (e) {
    return initialDocState;
  }
};
