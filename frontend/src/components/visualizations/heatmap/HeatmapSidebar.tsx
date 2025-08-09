import React from "react";
import { ToggleButton, ToggleButtonGroup } from "react-aria-components";
import type { DocState } from "../../../appState";
import type { ColorScaleArray } from "../../../colorScales";
import { formatTitle } from "../../../helpers";
import type { ColorScaleKey } from "../../../plotTypes";
import { Checkbox } from "../../primitives/Checkbox";
import { NumberInput } from "../../primitives/NumberInput";
import { Select, SelectItem } from "../../primitives/Select";
import { Slider } from "../../primitives/Slider";
import { Switch } from "../../primitives/Switch";
import { TitleField } from "../TitleField";

export const HeatmapSidebar = ({
  settings,
  updateSettings,
  sequences_count,
  continuousColorScales,
}: {
  settings: DocState["heatmap"];
  updateSettings: (values: Partial<DocState["heatmap"]>) => void;
  sequences_count: number;
  continuousColorScales: {
    [K in ColorScaleKey]: ColorScaleArray;
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

  const [colorScaleMode, setColorScaleMode] = React.useState<
    "continuous" | "discrete"
  >(settings.colorScaleKey === "Discrete" ? "discrete" : "continuous");
  const [continuousColorScaleKey, setContinuousColorScaleKey] =
    React.useState<ColorScaleKey>("Portland");

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
                <ToggleButtonGroup
                  selectionMode="single"
                  disallowEmptySelection={true}
                  selectedKeys={[colorScaleMode]}
                  onSelectionChange={(value) => {
                    const mode = value.values().next().value as
                      | "continuous"
                      | "discrete";
                    setColorScaleMode(mode);
                    updateSettings({
                      colorScaleKey:
                        mode === "continuous"
                          ? continuousColorScaleKey
                          : "Discrete",
                    });
                  }}
                >
                  <ToggleButton id="continuous">Continuous</ToggleButton>
                  <ToggleButton id="discrete">Discrete</ToggleButton>
                </ToggleButtonGroup>
              </div>
              {settings.colorScaleKey === "Discrete" ? (
                <>
                  <hr />
                  <div className="col-2">
                    <div className="field">
                      <NumberInput
                        label="Cutoff 1"
                        field="cutoff_1"
                        value={settings.cutoff_1}
                        updateValue={updateSettings}
                        min={settings.cutoff_2}
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
                        max={settings.cutoff_1}
                        step={1}
                      />
                    </div>
                  </div>
                </>
              ) : (
                <div className="field">
                  <Select
                    wide
                    id="colorscale"
                    items={Object.keys(continuousColorScales).map((name) => ({
                      id: name,
                      name,
                    }))}
                    selectedKey={settings.colorScaleKey}
                    onSelectionChange={(value) => {
                      setContinuousColorScaleKey(value as ColorScaleKey);
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
                              background: `linear-gradient(to right, ${continuousColorScales[item.name as ColorScaleKey].map((v) => v[1]).join(", ")})`,
                            }}
                          />
                          <span>{formatTitle(item.name)}</span>
                        </div>
                      </SelectItem>
                    )}
                  </Select>
                </div>
              )}
              <hr />
              <div className="field">
                <Checkbox
                  isSelected={settings.reverse}
                  onChange={(value) => {
                    updateSettings({
                      reverse: value,
                    });
                  }}
                >
                  Invert scale
                </Checkbox>
              </div>
              <hr />
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
              <div className="field">
                <label htmlFor="annotation-precision">
                  Precision
                  <ToggleButtonGroup
                    data-compact
                    selectionMode="single"
                    disallowEmptySelection={true}
                    selectedKeys={[settings.annotation_rounding]}
                    isDisabled={!settings.annotation}
                    onSelectionChange={(value) =>
                      updateSettings({
                        annotation_rounding: value.values().next()
                          .value as typeof settings.annotation_rounding,
                      })
                    }
                  >
                    <ToggleButton id={0}>0</ToggleButton>
                    <ToggleButton id={1}>1</ToggleButton>
                    <ToggleButton id={2}>2</ToggleButton>
                  </ToggleButtonGroup>
                </label>
              </div>
              <hr />
              <div className="field mask-values-field">
                <Checkbox
                  isSelected={settings.hideValuesBelowEnabled}
                  onChange={(value) => {
                    updateSettings({
                      hideValuesBelowEnabled: value,
                      hideValuesBelow: value ? 70 : 0,
                    });
                  }}
                >
                  Mask values
                </Checkbox>
                {settings.hideValuesBelowEnabled && (
                  <div className="mask-values-input-group">
                    <input
                      type="number"
                      className="mask-values-input"
                      value={settings.hideValuesBelow}
                      onChange={(e) =>
                        updateSettings({
                          hideValuesBelow: Number(e.target.value),
                        })
                      }
                      min={0}
                      max={100}
                      step={1}
                    />
                    <span className="mask-values-percent">%</span>
                  </div>
                )}
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
              Title
            </Switch>
            <div
              className="drawer"
              data-hidden={!settings.showTitles}
              aria-hidden={!settings.showTitles}
            >
              <TitleField
                textId="title"
                textValue={settings.title}
                onTextChange={(value) => updateSettings({ title: value })}
                fontValue={settings.titleFont}
                onFontChange={(value) => updateSettings({ titleFont: value })}
              />
            </div>
          </div>
        </div>
      </div>
    </div>
  );
};
