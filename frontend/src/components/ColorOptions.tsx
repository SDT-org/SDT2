import type React from "react";
import type { Color } from "react-aria-components";
import type { ColorString } from "../colors";
import { ColorPicker } from "./ColorPicker";
import { Select, SelectItem } from "./Select";

export interface DiscreteColorScheme {
  name: string;
  colors: [ColorString, ColorString, ColorString];
  description?: string;
}

// Colorblind-friendly discrete color schemes
export const discreteColorSchemes: DiscreteColorScheme[] = [
  {
    name: "Viridis",
    colors: ["hsl(68, 100%, 48%)", "hsl(184, 94%, 44%)", "hsl(253, 100%, 33%)"],
    description: "Colorblind-friendly scheme based on Viridis",
  },
  {
    name: "Cividis",
    colors: ["hsl(60, 100%, 65%)", "hsl(180, 17%, 46%)", "hsl(220, 100%, 15%)"],
    description: "Optimized for colorblind accessibility",
  },
  {
    name: "Okabe-Ito",
    colors: ["hsl(26, 100%, 50%)", "hsl(202, 100%, 50%)", "hsl(142, 53%, 40%)"],
    description: "Designed for maximum colorblind distinction",
  },
  {
    name: "Paul Tol",
    colors: ["hsl(201, 100%, 35%)", "hsl(329, 86%, 65%)", "hsl(55, 100%, 60%)"],
    description: "High contrast colorblind-safe palette",
  },
  {
    name: "Wong",
    colors: ["hsl(26, 100%, 50%)", "hsl(201, 100%, 42%)", "hsl(142, 70%, 35%)"],
    description: "Wong's colorblind-friendly palette",
  },
  {
    name: "Custom",
    colors: ["hsl(0, 100%, 50%)", "hsl(120, 100%, 50%)", "hsl(240, 100%, 50%)"],
    description: "Define your own colors",
  },
];

interface DiscreteColorOptionsProps {
  selectedScheme: string;
  customColors: [ColorString, ColorString, ColorString];
  onSchemeChange: (scheme: string) => void;
  onCustomColorChange: (index: 0 | 1 | 2, color: ColorString) => void;
}

export const DiscreteColorOptions: React.FC<DiscreteColorOptionsProps> = ({
  selectedScheme,
  customColors,
  onSchemeChange,
  onCustomColorChange,
}) => {
  const currentScheme =
    discreteColorSchemes.find((s) => s.name === selectedScheme) ||
    discreteColorSchemes[0];

  return (
    <div className="discrete-color-options">
      <div className="field">
        <Select
          wide
          id="discrete-color-scheme"
          label="Color Scheme"
          items={discreteColorSchemes.map((scheme) => ({
            id: scheme.name,
            name: scheme.name,
            description: scheme.description,
          }))}
          selectedKey={selectedScheme}
          onSelectionChange={(value) => onSchemeChange(value as string)}
        >
          {(item) => (
            <SelectItem key={item.name} textValue={item.name}>
              <div className="color-scheme-item">
                <div className="color-scheme-preview">
                  {discreteColorSchemes
                    .find((s) => s.name === item.name)
                    ?.colors.map((color, i) => (
                      <span
                        key={`${item.name}-color-${i}`}
                        className="color-swatch"
                        style={{ backgroundColor: color }}
                      />
                    ))}
                </div>
                <div className="color-scheme-info">
                  <span className="color-scheme-name">{item.name}</span>
                  {item.description && (
                    <span className="color-scheme-description">
                      {item.description}
                    </span>
                  )}
                </div>
              </div>
            </SelectItem>
          )}
        </Select>
      </div>

      {selectedScheme === "Custom" && (
        <>
          <hr />
          <div className="custom-colors-section">
            <div className="sublabel">Custom Colors</div>
            <div className="custom-color-pickers">
              <ColorPicker
                label="High"
                value={customColors[0]}
                onChange={(color: Color) =>
                  onCustomColorChange(0, color.toString() as ColorString)
                }
                aria-label="High value color"
              />
              <ColorPicker
                label="Medium"
                value={customColors[1]}
                onChange={(color: Color) =>
                  onCustomColorChange(1, color.toString() as ColorString)
                }
                aria-label="Medium value color"
              />
              <ColorPicker
                label="Low"
                value={customColors[2]}
                onChange={(color: Color) =>
                  onCustomColorChange(2, color.toString() as ColorString)
                }
                aria-label="Low value color"
              />
            </div>
          </div>
        </>
      )}

      <hr />
      <div className="color-preview-section">
        <div className="sublabel">Preview</div>
        <div className="color-preview-bar">
          {(selectedScheme === "Custom"
            ? customColors
            : currentScheme?.colors || customColors
          ).map((color, i) => (
            <div
              key={`preview-${i}-${color}`}
              className="color-preview-segment"
              style={{ backgroundColor: color }}
            >
              <span className="color-preview-label">
                {i === 0 ? "≥95%" : i === 1 ? "75-95%" : "<75%"}
              </span>
            </div>
          ))}
        </div>
      </div>
    </div>
  );
};
