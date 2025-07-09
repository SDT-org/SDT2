import { Label, ToggleButton, ToggleButtonGroup } from "react-aria-components";
import {
  TbGrid4X4,
  TbLetterA,
  TbNumber10,
  TbTableDashed,
} from "react-icons/tb";
import type { ColorString } from "../../../colors";
import type { DistributionState } from "../../../distributionState";
import { TitleField } from "../../fields/TitleField";
import { ColorPicker } from "../../ui/ColorPicker";
import { Slider } from "../../ui/Slider";
import { Switch } from "../../ui/Switch";
import { Tooltip } from "../../ui/Tooltip";

export const HistogramSidebar = ({
  sidebarComponent,
  settings,
  updateSettings,
}: {
  sidebarComponent?: React.ReactNode;
  settings: DistributionState["histogram"];
  updateSettings: React.Dispatch<Partial<DistributionState["histogram"]>>;
}) => {
  return (
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
                    <TbGrid4X4 size={18} />
                  </ToggleButton>
                </Tooltip>
                <Tooltip tooltip="Toggle axis lines">
                  <ToggleButton
                    id="showAxisLines"
                    aria-label="Toggle axis lines"
                  >
                    <TbTableDashed size={18} />
                  </ToggleButton>
                </Tooltip>
                <Tooltip tooltip="Toggle axis labels">
                  <ToggleButton
                    id="showAxisLabels"
                    aria-label="Toggle axis labels"
                  >
                    <TbLetterA size={18} />
                  </ToggleButton>
                </Tooltip>
                <Tooltip tooltip="Toggle tick values">
                  <ToggleButton
                    id="showTickLabels"
                    aria-label="Toggle axis tick values"
                  >
                    <TbNumber10 size={18} />
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
              <hr className="compact" />
              <div className="col-2 auto-onefr align-items-center color-slider-gap">
                <ColorPicker
                  value={settings.binColor}
                  onChange={(value) => {
                    updateSettings({
                      binColor: value.toString() as ColorString,
                    });
                  }}
                />
                <Slider
                  label="Bin Width"
                  defaultValue={settings.binSize}
                  onChange={(value) => updateSettings({ binSize: value })}
                  minValue={0.5}
                  maxValue={5}
                  step={0.5}
                />
              </div>
              <hr className="compact" />
              <div className="col-2 auto-onefr align-items-center color-slider-gap">
                <ColorPicker
                  value={settings.histlineColor}
                  onChange={(value) => {
                    updateSettings({
                      histlineColor: value.toString() as ColorString,
                    });
                  }}
                />
                <Slider
                  label="Outline Width"
                  defaultValue={settings.histOutlineWidth}
                  onChange={(value) =>
                    updateSettings({ histOutlineWidth: value })
                  }
                  minValue={0}
                  maxValue={5}
                  step={0.5}
                />
              </div>
              <hr className="compact" />
              <Slider
                label="Bar Gap"
                defaultValue={settings.barGap}
                onChange={(value) => updateSettings({ barGap: value })}
                minValue={0}
                maxValue={1}
                step={0.05}
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
