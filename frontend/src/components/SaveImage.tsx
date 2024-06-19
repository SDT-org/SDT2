import React from "react";
import Plotly from "../../vendor/plotly-custom.min.js";
import useAppState, { AppState } from "../appState";

export const SaveImage = ({ plotType }: { plotType?: string }) => {
  const { appState, setAppState } = useAppState();

  const handleSave = async () => {
    const element = document.getElementsByClassName("js-plotly-plot")[0];
    const config = {
      format: appState.client.saveFormat,
      width: 1000,
      height: 800,
    };
    // Run a warmup round to avoid the heatmap data missing from the render
    await Plotly.toImage(element, config);
    Plotly.toImage(element, config).then((data: string) => {
      window.pywebview.api.save_image({ data, format: config.format });
    });
  };

  return (
    <div className="group">
      <div className="field">
        <label htmlFor="save-format" className="header">
          Image Format
        </label>
        <select
          onChange={(e) =>
            setAppState((previous) => ({
              ...previous,
              client: {
                ...previous.client,
                saveFormat: e.target.value as AppState["client"]["saveFormat"],
              },
            }))
          }
        >
          <option value="png" selected={appState.client.saveFormat === "png"}>
            PNG
          </option>
          <option value="jpeg" selected={appState.client.saveFormat === "jpeg"}>
            JPEG
          </option>
          <option value="svg" selected={appState.client.saveFormat === "svg"}>
            SVG
          </option>
        </select>
      </div>
      <button
        type="button"
        aria-label="Save the plot as an image"
        onClick={handleSave}
      >
        Export{plotType ? ` ${plotType}` : null}...
      </button>
    </div>
  );
};
