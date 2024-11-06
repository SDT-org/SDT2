import React from "react";
import HistogramPlot from "./HistogramPlot";
import { NumberInput } from "./NumberInput";
import { DistributionData, DistributionSettings } from "../plotTypes";

enum ColorOption {
  White = "white",
  Black = "black",
  Red = "tomato",
  Blue = "lightblue",
  Green = "lightgreen",
  Purple = "plum",
  Pink = "lightcoral",
}

export const Distribution = ({
  data,
  settings,
  updateSettings,
  footer,
}: {
  data: DistributionData | undefined;
  settings: DistributionSettings;
  updateSettings: (_: Partial<DistributionSettings>) => void;
  footer?: React.ReactNode;
}) => {
  if (!data) {
    return (
      <div className="app-main centered">
        <h1>No data available</h1>
      </div>
    );
  }

  // Assign x and y explicitly
  const { x, y,raw_mat } = data;

  return (
    <>
      {/*Sidebar */}
      <div className="app-sidebar">
        <div className="app-sidebar-toolbar">
          <div className="form">
            <div className="group">
              <div className="field">
                <label className="header">Title</label>
                <input
                  type="text"
                  value={settings.plotTitle}
                  onChange={(e) =>
                    updateSettings({ plotTitle: e.target.value })
                  }
                />
              </div>
              <div className="row">
                <div className="col-2">
                  <div className="field">
                    <label htmlFor="showGrid">
                      <input
                        type="checkbox"
                        name="showGrid"
                        id="showGrid"
                        checked={settings.showGrid}
                        onChange={() =>
                          updateSettings({
                            showGrid: !settings.showGrid,
                          })
                        }
                      />
                      Grid
                    </label>
                  </div>
                  <div className="field">
                    <label htmlFor="showTickLabels">
                      <input
                        type="checkbox"
                        name="showTickLabels"
                        id="showTickLabels"
                        checked={settings.showTickLabels}
                        onChange={() =>
                          updateSettings({
                            showTickLabels: !settings.showTickLabels,
                          })
                        }
                      />
                      Tick Labels
                    </label>
                  </div>
                </div>
              </div>
              <div className="col-2">
                <div className="field">
                  <label htmlFor="showLine">
                    <input
                      type="checkbox"
                      name="showLine"
                      id="showLine"
                      checked={settings.showLine}
                      onChange={() =>
                        updateSettings({
                          showLine: !settings.showLine,
                        })
                      }
                    />
                    Axis Lines
                  </label>
                </div>
                <div className="field">
                  <label htmlFor="showAxisLabels">
                    <input
                      type="checkbox"
                      name="showAxisLabels"
                      id="showAxisLabels"
                      checked={settings.showAxisLabels}
                      onChange={() =>
                        updateSettings({
                          showAxisLabels: !settings.showAxisLabels,
                        })
                      }
                    />
                    Axis Title
                  </label>
                </div>
              </div>
            </div>
            <div className="group">
              <div className="field">
                <label className="header">
                  <input
                    type="checkbox"
                    checked={settings.showHistogram}
                    onChange={(e) =>
                      updateSettings({ showHistogram: e.target.checked })
                    }
                  />
                  Bar Plot
                </label>
              </div>
              <div className="field">
                <label htmlFor="bar-color">Color</label>
                <select
                  id="bar-color"
                  value={settings.color}
                  disabled={!settings.showHistogram}
                  onChange={(e) =>
                    updateSettings({ color: e.target.value })
                  }
                >
                  {Object.entries(ColorOption).map(([key, value]) => (
                    <option key={key} value={value}>
                      {key}
                    </option>
                  ))}
                </select>
              </div>
              <div className="col-2">
                <div className="field">
                  <label htmlFor="bar-line-color">Outline</label>
                  <select
                    id="bar-line-color"
                    value={settings.barlineColor}
                    disabled={!settings.showHistogram}
                    onChange={(e) =>
                      updateSettings({ barlineColor: e.target.value })
                    }
                  >
                    {Object.entries(ColorOption).map(([key, value]) => (
                      <option key={key} value={value}>
                        {key}
                      </option>
                    ))}
                  </select>
                </div>
                <NumberInput
                  label="Outline Width"
                  field="barOutlineWidth"
                  value={settings.barOutlineWidth}
                  updateValue={(field, value) =>
                    updateSettings({ [field]: value })
                  }
                  min={0}
                  max={5}
                  step={1}
                  disabled={!settings.showHistogram}
                />
              </div>
            </div>

            {/* Placeholder for Future Plot Settings */}
            {/* You can add more groups here for other plots like Line Plot, Scatter Plot, etc. */}
          </div>
        </div>

        {/* Footer area for additional content or controls */}
        <div className="app-sidebar-footer">{footer}</div>
      </div>

      {/* Main area where the Plotly charts are rendered */}
      <div className="app-main">
        {/* Render HistogramPlot conditionally */}
        {settings.showHistogram && (
          <HistogramPlot
            x={x}
            y={y}
            raw_mat={raw_mat}
            settings={{
              color: settings.barColor,
              outlineColor: settings.barlineColor,
              outlineWidth: settings.barOutlineWidth,
            }}
          />
        )}
        {/* Future Plots can be rendered here similarly */}
      </div>
    </>
  );
};
