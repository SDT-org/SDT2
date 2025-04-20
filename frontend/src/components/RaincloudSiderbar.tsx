import {
  Input,
  Label,
  TextField,
  ToggleButton,
  ToggleButtonGroup,
} from "react-aria-components";
import type { ColorString } from "../colors";
import type { DistributionState } from "../distributionState";
import { ColorPicker } from "./ColorPicker";
import { Select, SelectItem } from "./Select";
import { Slider } from "./Slider";
import { Switch } from "./Switch";
import { Tooltip } from "./Tooltip";

export const RaincloudSidebar = ({
  sidebarComponent,
  settings,
  updateSettings,
}: {
  sidebarComponent?: React.ReactNode;
  settings: DistributionState["raincloud"];
  updateSettings: React.Dispatch<Partial<DistributionState["raincloud"]>>;
}) => (
  <div className="app-sidebar app-sidebar-right">
    <div className="app-sidebar-toolbar">
      <div className="form">
        {sidebarComponent}
        <div className="group">
          <div className="drawer">
            <ToggleButtonGroup
              data-icon-only
              selectionMode="multiple"
              selectedKeys={Object.keys(settings).filter(
                (key) =>
                  [
                    "showGrid",
                    "showTickLabels",
                    "showAxisLines",
                    "showAxisLabels",
                    "showMeanline",
                  ].includes(key) && settings[key as keyof typeof settings],
              )}
              onSelectionChange={(value) =>
                updateSettings({
                  showGrid: value.has("showGrid"),
                  showTickLabels: value.has("showTickLabels"),
                  showAxisLines: value.has("showAxisLines"),
                  showAxisLabels: value.has("showAxisLabels"),
                  showMeanline: value.has("showMeanline"),
                })
              }
            >
              <Tooltip tooltip="Toggle grid">
                <ToggleButton id="showGrid" aria-label="Toggle grid">
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
                      <path d="M9 3v18M15 3v18M3 15h18M21 9H3M19 21H5a2 2 0 0 1-2-2V5a2 2 0 0 1 2-2h14a2 2 0 0 1 2 2v14a2 2 0 0 1-2 2z" />
                    </g>
                  </svg>
                </ToggleButton>
              </Tooltip>
              <Tooltip tooltip="Toggle axis lines">
                <ToggleButton id="showAxisLines" aria-label="Toggle axis lines">
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
                      <path d="M3 5a2 2 0 0 1 2 -2h14a2 2 0 0 1 2 2v14a2 2 0 0 1-2 2h-14a2 2 0 0 1-2 -2v-14z" />
                      <path d="M3 10h18" />
                      <path d="M10 3v18" />
                    </g>
                  </svg>
                </ToggleButton>
              </Tooltip>
              <Tooltip tooltip="Toggle tick values">
                <ToggleButton
                  id="showTickLabels"
                  aria-label="Toggle tick values"
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
                      <path d="M18 22H6a4 4 0 0 1-4-4V6a4 4 0 0 1 4-4h12a4 4 0 0 1 4 4v12a4 4 0 0 1-4 4zM9 7v10M6 9l3-2" />
                      <path d="M15.5 17a2.5 2.5 0 0 1-2.5-2.5v-5a2.5 2.5 0 1 1 5 0v5a2.5 2.5 0 0 1-2.5 2.5z" />
                    </g>
                  </svg>
                </ToggleButton>
              </Tooltip>
              <Tooltip tooltip="Toggle mean line">
                <ToggleButton id="showMeanline" aria-label="Toggle mean line">
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
                      <path d="M5 12h2" />
                      <path d="M17 12h2" />
                      <path d="M11 12h2" />
                    </g>
                  </svg>
                </ToggleButton>
              </Tooltip>
            </ToggleButtonGroup>
            <hr className="compact" />
            <div className="col-2 auto-onefr align-items-center color-slider-gap">
              <ColorPicker
                value={settings.fillColor}
                onChange={(value) => {
                  updateSettings({
                    fillColor: value.toString() as ColorString,
                  });
                }}
              />
              <Slider
                label="Bandwidth"
                defaultValue={settings.bandwidth}
                onChange={(value) => updateSettings({ bandwidth: value })}
                minValue={1}
                maxValue={20}
                step={1}
              />
            </div>
            <hr className="compact" />
            <div className="col-2 auto-onefr align-items-center color-slider-gap">
              <ColorPicker
                value={settings.lineColor}
                onChange={(value) =>
                  updateSettings({
                    lineColor: value.toString() as ColorString,
                  })
                }
              />
              <Slider
                label="Width"
                value={settings.lineWidth}
                onChange={(value) => updateSettings({ lineWidth: value })}
                minValue={1}
                maxValue={20}
                step={1}
              />
            </div>
          </div>
        </div>
        <div className="group">
          <h4 className="setting-header">Points</h4>
          <div
            className="drawer"
            data-hidden={!settings.showPoints}
            aria-hidden={!settings.showPoints}
          >
            <div className="col-2 auto-onefr align-items-center">
              <Label htmlFor="points-type">Type</Label>
              <Select
                data-compact
                id="points-type"
                selectedKey={settings.points}
                onSelectionChange={(value) => {
                  updateSettings({
                    points: value as typeof settings.points,
                  });
                }}
                items={Object.entries({
                  all: "All",
                  outliers: "Outliers",
                  suspectedoutliers: "Suspected Outliers",
                }).map(([id, name]) => ({
                  id,
                  name,
                }))}
              >
                {(item) => (
                  <SelectItem textValue={item.name}>{item.name}</SelectItem>
                )}
              </Select>
            </div>
            <Slider
              label="Size"
              value={settings.markerSize}
              onChange={(value) => updateSettings({ markerSize: value })}
              minValue={1}
              maxValue={20}
              step={1}
            />
            <Slider
              label="Position"
              value={settings.pointPos}
              onChange={(value) => updateSettings({ pointPos: value })}
              minValue={-2}
              maxValue={-0.5}
              step={0.1}
            />
            <Slider
              label="Jitter"
              value={settings.jitter}
              onChange={(value) => updateSettings({ jitter: value })}
              minValue={0.1}
              maxValue={1}
              step={0.1}
            />
            <div className="col-2 onefr-auto small-color align-items-center">
              <Label>Color</Label>
              <ColorPicker
                value={settings.markerColor}
                onChange={(value) => {
                  updateSettings({
                    markerColor: value.toString() as ColorString,
                  });
                }}
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
);
