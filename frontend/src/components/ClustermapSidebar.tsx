import React from "react";
import { Input, Label, TextField } from "react-aria-components";
import type { DocState } from "../appState";
import type { HeatmapSettings } from "../plotTypes";
import { NumberInput } from "./NumberInput";
import { Select, SelectItem } from "./Select";
import { Slider } from "./Slider";
import { Switch } from "./Switch";

export const ClustermapSidebar = ({
  settings,
  updateSettings,
  sequences_count,
}: {
  settings: DocState["clustermap"];
  updateSettings: (values: Partial<DocState["clustermap"]>) => void;
  sequences_count: number;
}) => {
  const maybeWarnPerformance = React.useCallback(
    (enabled: boolean, fn: () => void) => {
      if (
        enabled &&
        sequences_count > 99 &&
        !confirm(
          "Warning: Enabling this setting may significantly impact render performance.",
        )
      ) {
        return;
      }
      fn();
    },
    [sequences_count],
  );

  return (
    <div className="app-sidebar app-sidebar-right heatmap-sidebar">
      <div className="app-sidebar-toolbar">
        <div className="form">
          <div className="group">
            <label className="setting-header" htmlFor="colorscale">
              Clusters
            </label>
            <div className="drawer">
              <div className="col-2">
                <div className="field">
                  <NumberInput
                    label="Threshold 1"
                    field="threshold_one"
                    value={settings.threshold_one}
                    updateValue={updateSettings}
                    min={settings.threshold_two + 1}
                    max={100}
                    step={1}
                  />
                </div>
                <div className="field">
                  <NumberInput
                    label="Threshold 2"
                    field="threshold_two"
                    value={settings.threshold_two}
                    updateValue={updateSettings}
                    min={0}
                    max={settings.threshold_one - 1}
                    step={1}
                  />
                </div>
              </div>
              <hr className="compact" />
              <div className="field">
                <Slider
                  label="Cell Spacing"
                  labelClassName="sublabel"
                  id="cellspace"
                  onChange={(value) => updateSettings({ cellspace: value })}
                  minValue={0}
                  maxValue={20}
                  value={settings.cellspace}
                />
              </div>
            </div>
          </div>
          <div className="group">
            <Switch
              isSelected={settings.annotation}
              onChange={(value) =>
                maybeWarnPerformance(value, () =>
                  updateSettings({
                    annotation: value,
                  }),
                )
              }
            >
              Percent Identities
            </Switch>
            <div
              className="drawer"
              data-hidden={!settings.annotation}
              aria-hidden={!settings.annotation}
            >
              <div className="col-2">
                <div className="field">
                  <label htmlFor="round-vals">Precision</label>
                  <select
                    id="round-vals"
                    value={settings.annotation_rounding}
                    onChange={(event) =>
                      updateSettings({
                        annotation_rounding: Number.parseInt(
                          event.target.value,
                        ) as HeatmapSettings["annotation_rounding"],
                      })
                    }
                    disabled={!settings.annotation}
                  >
                    <option value={0}>0</option>
                    <option value={1}>1</option>
                    <option value={2}>2</option>
                  </select>
                </div>
              </div>
            </div>
          </div>
          <div className="group">
            <Switch
              isSelected={settings.axis_labels}
              onChange={(value) =>
                maybeWarnPerformance(value, () =>
                  updateSettings({
                    axis_labels: value,
                  }),
                )
              }
            >
              Axis Labels
            </Switch>
            <div
              className="drawer"
              data-hidden={!settings.axis_labels}
              aria-hidden={!settings.axis_labels}
            >
              <div className="col-2">
                <NumberInput
                  label="Font Size"
                  field="axlabel_xfontsize"
                  value={settings.axlabel_xfontsize}
                  updateValue={updateSettings}
                  min={1}
                  max={40}
                  step={1}
                />
                <NumberInput
                  label="X Rotation"
                  field="axlabel_xrotation"
                  value={settings.axlabel_xrotation}
                  updateValue={updateSettings}
                  min={-90}
                  max={90}
                  step={10}
                />
                <NumberInput
                  label="Y Rotation"
                  field="axlabel_yrotation"
                  value={settings.axlabel_yrotation}
                  updateValue={updateSettings}
                  min={-90}
                  max={90}
                  step={10}
                />
              </div>
            </div>
          </div>

          <div className="group">
            <Switch
              isSelected={settings.showTitles}
              onChange={(value) => {
                updateSettings({
                  showTitles: value,
                });
              }}
            >
              Plot Titles
            </Switch>
            <div
              className="drawer"
              data-hidden={!settings.showTitles}
              aria-hidden={!settings.showTitles}
            >
              <div className="col-2 auto-onefr align-items-center">
                <Label htmlFor="font">Font Type</Label>
                <Select
                  id="font"
                  data-compact
                  selectedKey={settings.titleFont}
                  onSelectionChange={(value) =>
                    updateSettings({
                      titleFont: value as typeof settings.titleFont,
                    })
                  }
                  items={["Sans Serif", "Monospace"].map((name) => ({
                    id: name,
                    name,
                  }))}
                >
                  {(item) => (
                    <SelectItem textValue={item.name}>{item.name}</SelectItem>
                  )}
                </Select>
              </div>

              <div className="field">
                <TextField
                  onChange={(value) => updateSettings({ title: value })}
                  value={settings.title}
                >
                  <Label>Title</Label>
                  <Input />
                </TextField>
              </div>
            </div>
          </div>
        </div>
      </div>
    </div>
  );
};
