import React from "react";
import { Input, Label, Text, TextField } from "react-aria-components";
import type { DocState } from "../appState";
import { reorderMethods } from "../constants";
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

  const [tempThreshold, setTempThreshold] = React.useState(settings.threshold);

  return (
    <div className="app-sidebar app-sidebar-right heatmap-sidebar">
      <div className="app-sidebar-toolbar">
        <div className="form">
          <div className="group">
            <label className="setting-header" htmlFor="colorscale">
              Clusters
            </label>
            <div className="drawer">
              <div className="col-2 auto-onefr align-items-center">
                <Label htmlFor="method">Linkage Method</Label>
                <Select
                  id="method"
                  data-compact
                  selectedKey={settings.method}
                  onSelectionChange={(value: React.Key) =>
                    updateSettings({
                      ...settings,
                      method: value as typeof settings.method,
                    })
                  }
                  items={Object.entries(reorderMethods).map(([key, value]) => ({
                    id: key,
                    name: value.name,
                    description: value.description,
                  }))}
                >
                  {(item) => (
                    <SelectItem textValue={item.name}>
                      <Text slot="label">{item.name}</Text>
                      <Text slot="description">{item.description}</Text>
                    </SelectItem>
                  )}
                </Select>
              </div>

              <hr className="compact" />

              <Slider
                label="Threshold"
                id="threshold"
                onChange={setTempThreshold}
                onChangeEnd={(value) => updateSettings({ threshold: value })}
                minValue={0}
                maxValue={100}
                step={1}
                value={tempThreshold}
              />

              <Slider
                label="Cell Spacing"
                id="cellspace"
                onChange={(value) => updateSettings({ cellspace: value })}
                minValue={0}
                maxValue={20}
                value={settings.cellspace}
              />
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
              <Slider
                label="Font Size"
                labelClassName="sublabel"
                onChange={(value) =>
                  updateSettings({ axlabel_fontsize: value })
                }
                value={settings.axlabel_fontsize}
                minValue={1}
                maxValue={20}
                step={1}
              />
              <Slider
                label="X Rotation"
                labelClassName="sublabel"
                onChange={(value) =>
                  updateSettings({ axlabel_xrotation: value })
                }
                value={settings.axlabel_xrotation}
                minValue={-90}
                maxValue={90}
                step={10}
              />
              <Slider
                label="Y Rotation"
                labelClassName="sublabel"
                onChange={(value) =>
                  updateSettings({ axlabel_yrotation: value })
                }
                value={settings.axlabel_yrotation}
                minValue={-90}
                maxValue={90}
                step={10}
              />
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
              Plot Title
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
                  onSelectionChange={(value: React.Key) =>
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

          <div className="group">
            <Switch
              isSelected={settings.showLegend}
              onChange={(value) => {
                updateSettings({
                  showLegend: value,
                });
              }}
            >
              Legend
            </Switch>
          </div>
        </div>
      </div>
    </div>
  );
};
