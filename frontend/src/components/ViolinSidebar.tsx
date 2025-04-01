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

export const ViolinSidebar = ({
  sidebarComponent,
  settings,
  updateSettings,
}: {
  sidebarComponent?: React.ReactNode;
  settings: DistributionState["violin"];
  updateSettings: React.Dispatch<Partial<DistributionState["violin"]>>;
}) => (
  <div className="app-sidebar app-sidebar-right">
    <div className="app-sidebar-toolbar form">
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
                    <path d="M3 5a2 2 0 0 1 2 -2h14a2 2 0 0 1 2 2v14a2 2 0 0 1 -2 2h-14a2 2 0 0 1-2 -2v-14z" />
                    <path d="M3 10h18" />
                    <path d="M10 3v18" />
                  </g>
                </svg>
              </ToggleButton>
            </Tooltip>
            <Tooltip tooltip="Toggle tick values">
              <ToggleButton id="showTickLabels" aria-label="Toggle tick values">
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

          <Label>Orientation</Label>
          <ToggleButtonGroup
            data-compact
            selectionMode="single"
            disallowEmptySelection={true}
            selectedKeys={[settings.plotOrientation]}
            onSelectionChange={(value) =>
              updateSettings({
                plotOrientation: value.values().next()
                  .value as typeof settings.plotOrientation,
              })
            }
          >
            <ToggleButton id="vertical">Vertical</ToggleButton>
            <ToggleButton id="horizontal">Horizontal</ToggleButton>
          </ToggleButtonGroup>
        </div>
      </div>
      <div className="group">
        <Switch
          isSelected={settings.showViolin}
          onChange={(value) =>
            updateSettings({
              showViolin: value,
              pointOrientation: !settings.showBox && value ? "Violin" : "Box",
            })
          }
        >
          Violin
        </Switch>
        <div
          className="drawer"
          data-hidden={!settings.showViolin}
          aria-hidden={!settings.showViolin}
        >
          <div className="col-2 auto-onefr align-items-center color-slider-gap">
            <ColorPicker
              className="medium-color"
              value={settings.fillColor}
              onChange={(value) => {
                updateSettings({
                  fillColor: value.toString() as ColorString,
                });
              }}
            />
            <Slider
              label="Band Width"
              defaultValue={settings.bandwidth}
              isDisabled={!settings.showViolin}
              onChange={(value) => updateSettings({ bandwidth: value })}
              minValue={1}
              maxValue={20}
              step={1}
            />
          </div>
          <hr className="compact" />
          <div className="col-2 auto-onefr align-items-center color-slider-gap">
            <ColorPicker
              className="medium-color"
              value={settings.lineColor}
              onChange={(value) =>
                updateSettings({
                  lineColor: value.toString() as ColorString,
                })
              }
            />
            <Slider
              label="Line Width"
              value={settings.lineWidth}
              isDisabled={!settings.showViolin}
              onChange={(value) => updateSettings({ lineWidth: value })}
              minValue={0}
              maxValue={20}
              step={1}
            />
          </div>
        </div>
      </div>
      <div className="group">
        <Switch
          isSelected={settings.showBox}
          onChange={(value) => {
            updateSettings({
              showBox: value,
              pointOrientation:
                !settings.showViolin && value ? "Box" : "Violin",
            });
          }}
        >
          Box
        </Switch>
        <div
          className="drawer"
          data-hidden={!settings.showBox}
          aria-hidden={!settings.showBox}
        >
          <div className="col-2 auto-onefr align-items-center color-slider-gap">
            <ColorPicker
              className="medium-color"
              value={settings.boxfillColor}
              onChange={(value) => {
                updateSettings({
                  boxfillColor: value.toString() as ColorString,
                });
              }}
            />
            <Slider
              label="Box Width"
              value={settings.boxWidth}
              isDisabled={!settings.showBox}
              onChange={(value) => updateSettings({ boxWidth: value })}
              minValue={0.05}
              maxValue={0.5}
              step={0.05}
            />
          </div>
          <hr className="compact" />
          <div className="col-2 auto-onefr align-items-center color-slider-gap">
            <ColorPicker
              className="medium-color"
              value={settings.boxlineColor}
              onChange={(value) => {
                updateSettings({
                  boxlineColor: value.toString() as ColorString,
                });
              }}
            />
            <Slider
              label="Line Width"
              value={settings.boxlineWidth}
              isDisabled={!settings.showBox}
              onChange={(value) => updateSettings({ boxlineWidth: value })}
              minValue={0}
              maxValue={20}
              step={1}
            />
          </div>
          <hr className="compact" />
          <Slider
            label="Whiskers"
            value={settings.whiskerWidth}
            isDisabled={!settings.showBox}
            onChange={(value) => updateSettings({ whiskerWidth: value })}
            minValue={0.1}
            maxValue={1}
            step={0.1}
          />
        </div>
      </div>
      <div className="group">
        <Switch
          isSelected={settings.showPoints}
          onChange={(value) => updateSettings({ showPoints: value })}
        >
          Points
        </Switch>
        <div
          className="drawer"
          data-hidden={!settings.showPoints}
          aria-hidden={!settings.showPoints}
        >
          <ToggleButtonGroup
            data-compact
            selectionMode="single"
            disallowEmptySelection={true}
            selectedKeys={
              settings.showViolin || settings.showBox
                ? [settings.pointOrientation]
                : []
            }
            onSelectionChange={(value) =>
              updateSettings({
                pointOrientation: value.values().next()
                  .value as typeof settings.pointOrientation,
              })
            }
          >
            <ToggleButton id="Violin" isDisabled={!settings.showViolin}>
              Violin
            </ToggleButton>
            <ToggleButton id="Box" isDisabled={!settings.showBox}>
              Box
            </ToggleButton>
          </ToggleButtonGroup>

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

          <div className="col-2 auto-onefr align-items-center color-slider-gap">
            <ColorPicker
              className="medium-color"
              value={settings.markerColor}
              onChange={(value) => {
                updateSettings({
                  markerColor: value.toString() as ColorString,
                });
              }}
            />
            <Slider
              label="Size"
              value={settings.markerSize}
              isDisabled={!settings.showPoints}
              onChange={(value) => updateSettings({ markerSize: value })}
              minValue={1}
              maxValue={20}
              step={1}
            />
          </div>

          <Slider
            label="Position"
            value={settings.pointPos}
            isDisabled={!settings.showViolin && !settings.showBox}
            onChange={(value) => updateSettings({ pointPos: value })}
            minValue={-2}
            maxValue={2}
            step={0.1}
          />

          <Slider
            label="Jitter"
            value={settings.jitter}
            isDisabled={!settings.showPoints}
            onChange={(value) => updateSettings({ jitter: value })}
            minValue={0.1}
            maxValue={1}
            step={0.1}
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
);
