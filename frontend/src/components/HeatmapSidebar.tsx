import React from "react";
import { Input, Label, TextField, ToggleButton } from "react-aria-components";
import type { DocState } from "../appState";
import type { ColorScaleArray } from "../colorScales";
import { formatTitle } from "../helpers";
import type { ColorScaleKey, HeatmapSettings } from "../plotTypes";
import { NumberInput } from "./NumberInput";
import { Select, SelectItem } from "./Select";
import { Slider } from "./Slider";
import { Switch } from "./Switch";
import { Tooltip } from "./Tooltip";

export const HeatmapSidebar = ({
  settings,
  updateSettings,
  sequences_count,
  colorScales,
}: {
  settings: DocState["heatmap"];
  updateSettings: (values: Partial<DocState["heatmap"]>) => void;
  sequences_count: number;
  colorScales: {
    [K in HeatmapSettings["colorScaleKey"]]: ColorScaleArray;
  };
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
              Colorscale
            </label>

            <div className="drawer">
              <div className="field">
                <div className="col-2 onefr-atuo colorscale-setting">
                  <Select
                    id="colorscale"
                    items={Object.keys(colorScales).map((name) => ({
                      id: name,
                      name,
                    }))}
                    selectedKey={settings.colorScaleKey}
                    onSelectionChange={(value) => {
                      updateSettings({
                        colorScaleKey: value as ColorScaleKey,
                      });
                    }}
                  >
                    {(item) => (
                      <SelectItem key={item.name} textValue={item.name}>
                        <div className="color-scale-list-item">
                          <span
                            className="preview"
                            aria-hidden="true"
                            style={{
                              background: `linear-gradient(to right, ${colorScales[item.name as ColorScaleKey].map((v) => v[1]).join(", ")})`,
                            }}
                          />
                          <span>{formatTitle(item.name)}</span>
                        </div>
                      </SelectItem>
                    )}
                  </Select>
                  <Tooltip tooltip="Reverse colorscale" delay={600}>
                    <ToggleButton
                      aria-label="Toggle reverse colorscale"
                      id="reverse"
                      isSelected={settings.reverse}
                      onChange={(value) => {
                        updateSettings({
                          reverse: value,
                        });
                      }}
                    >
                      <svg
                        xmlns="http://www.w3.org/2000/svg"
                        viewBox="0 0 24 24"
                        aria-hidden="true"
                      >
                        <g
                          style={{
                            fill: "none",
                            stroke: "currentcolor",
                            strokeWidth: 2,
                            strokeLinecap: "round",
                            strokeLinejoin: "round",
                            strokeMiterlimit: 10,
                          }}
                        >
                          <path d="m15 1 5 5-5 5M9 23l-5-5 5-5" />
                          <path d="M4 18h13a6 6 0 0 0 6-6v-1M20 6H7a6 6 0 0 0-6 6v1" />
                        </g>
                      </svg>
                    </ToggleButton>
                  </Tooltip>
                </div>
              </div>

              <div className="field" style={{ paddingTop: "0.4rem" }}>
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
              {settings.colorScaleKey === "Discrete" ? (
                <div className="col-2">
                  <div className="field">
                    <NumberInput
                      label="Cutoff 1"
                      field="cutoff_1"
                      value={settings.cutoff_1}
                      updateValue={updateSettings}
                      min={settings.cutoff_2 + 1}
                      max={100}
                      step={1}
                    />
                  </div>
                  <div className="field">
                    <NumberInput
                      label="Cutoff 2"
                      field="cutoff_2"
                      value={settings.cutoff_2}
                      updateValue={updateSettings}
                      min={settings.vmin + 1}
                      max={settings.cutoff_1 - 1}
                      step={1}
                    />
                  </div>
                </div>
              ) : null}
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
              <Slider
                label="Font Size"
                labelClassName="sublabel"
                onChange={(value) =>
                  updateSettings({ axlabel_xfontsize: value })
                }
                value={settings.axlabel_xfontsize}
                minValue={1}
                maxValue={40}
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
              isSelected={settings.showscale}
              onChange={(value) =>
                updateSettings({
                  showscale: value,
                })
              }
            >
              Scale Bar
            </Switch>
            <div
              className="drawer"
              data-hidden={!settings.showscale}
              aria-hidden={!settings.showscale}
            >
              <div className="range-group">
                <Slider
                  label="Height"
                  minValue={1}
                  maxValue={10}
                  step={1}
                  value={settings.cbar_shrink}
                  onChange={(value) => updateSettings({ cbar_shrink: value })}
                />
                <Slider
                  label="Width"
                  id="cbar-aspect"
                  minValue={1}
                  maxValue={10}
                  step={1}
                  value={settings.cbar_aspect}
                  onChange={(value) => updateSettings({ cbar_aspect: value })}
                />
              </div>
              <div className="col-2">
                <NumberInput
                  label="Min"
                  field="vmin"
                  value={settings.vmin}
                  updateValue={updateSettings}
                  min={1}
                  max={settings.vmax}
                  step={1}
                  isDisabled={settings.colorScaleKey === "Discrete"}
                />
                <NumberInput
                  label="Max"
                  field="vmax"
                  value={settings.vmax}
                  updateValue={updateSettings}
                  min={settings.vmin}
                  max={100}
                  step={1}
                  isDisabled={settings.colorScaleKey === "Discrete"}
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
    </div>
  );
};
