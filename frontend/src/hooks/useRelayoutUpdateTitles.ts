import React from "react";

type ExtendedPlotRelayoutEvent = Plotly.PlotRelayoutEvent & {
  "title.text"?: string;
  "title.subtitle.text"?: string;
  "yaxis.title.text"?: string;
  "xaxis.title.text"?: string;
};

type UpdateSettingsPayload = Partial<{
  title: string;
  subtitle: string;
}>;

export const useRelayoutUpdateTitles = (
  updateSettings: (payload: UpdateSettingsPayload) => void,
): ((event: ExtendedPlotRelayoutEvent) => void) =>
  React.useCallback(
    (event: ExtendedPlotRelayoutEvent) => {
      const updates: UpdateSettingsPayload = {};

      if (event["title.text"]) updates.title = event["title.text"];
      if (event["title.subtitle.text"])
        updates.subtitle = event["title.subtitle.text"];

      if (Object.keys(updates).length > 0) {
        updateSettings(updates);
      }
    },
    [updateSettings],
  );

export const useRelayoutHideSubtitle = (hide: boolean) =>
  React.useEffect(() => {
    // There's a bug in plotly that keeps the subtitle rendering, hack to hide it
    const subtitleEl = document.getElementsByClassName("gtitle-subtitle")[0];
    if (!subtitleEl) {
      return;
    }

    if (hide) {
      subtitleEl.setAttribute("data-hidden", "true");
    } else {
      subtitleEl.removeAttribute("data-hidden");
    }
  }, [hide]);
