import * as d3 from "d3";
import tinycolor from "tinycolor2";
import { z } from "zod";
import type { ColorScaleArray } from "./colorScales";

export type ColorString =
  | `hsl(${number}, ${number}%, ${number}%)`
  | `hsla(${number}, ${number}%, ${number}%, ${number})`;

export const ColorStringSchema = z
  .string()
  .refine(
    (color): color is ColorString =>
      tinycolor(color).isValid() && /^hsl/.test(color),
    {
      message: "Invalid HSL or HSLA color string",
    },
  );

export function createD3ColorScale(
  colorArray: ColorScaleArray,
  discrete: boolean,
  max: number,
  min: number,
): d3.ScaleLinear<string, string> {
  const domain = discrete
    ? colorArray.map(([stop]) => stop)
    : colorArray.map(([stop]) => stop * (max - min) + min);

  const range = colorArray.map(([_, color]) => color);

  return d3
    .scaleLinear<string>()
    .domain(domain)
    .range(range)
    .interpolate(d3.interpolateRgb)
    .clamp(true);
}
//https://stackoverflow.com/questions/10014271/generate-random-color-distinguishable-to-humans
export const distinctColor = (index: number) => {
  if (!index) {
    return "hsl(0, 0%, 96%)"; // Light grey for index 0
  }

  const hue = (index * 137.5) % 360; // Golden angle approximation for hue

  // Define saturation and alpha patterns
  const patterns = [
    { saturation: 90, alpha: 50 }, // Vibrant, Medium
    { saturation: 75, alpha: 65 }, // Less Vibrant, Lighter
    { saturation: 90, alpha: 35 }, // Dark and full
    { saturation: 60, alpha: 75 }, // Less Vibrant, Lighter
    { saturation: 75, alpha: 40 }, // Dark
    { saturation: 85, alpha: 50 }, // Default
  ];

  // cycle with modulo
  const cycle = ((index % patterns.length) + patterns.length) % patterns.length;

  // Default values
  let saturation = 85;
  let alpha = 50;

  if (patterns[cycle]) {
    saturation = patterns[cycle].saturation;
    alpha = patterns[cycle].alpha;
  }

  return `hsl(${hue.toFixed(1)}, ${saturation}%, ${alpha}%)`;
};

type ColorName =
  | "White"
  | "Black"
  | "Tomato"
  | "LightBlue"
  | "LightGreen"
  | "Purple"
  | "Orange"
  | "Yellow"
  | "Cyan"
  | "Teal"
  | "Magenta"
  | "Brown"
  | "Lime"
  | "Coral"
  | "Turquoise"
  | "Indigo"
  | "Violet"
  | "Lavender"
  | "Peach"
  | "SkyBlue"
  | "Olive"
  | "Tan"
  | "Salmon"
  | "Maroon"
  | "Navy"
  | "Khaki"
  | "Periwinkle"
  | "Mint"
  | "Azure"
  | "Chartreuse"
  | "Goldrod"
  | "SlateBlue"
  | "LightSeaGreen"
  | "DarkCyan"
  | "RosyBrown"
  | "PaleVioletRed"
  | "DeepPink"
  | "DarkOrange"
  | "Crimson"
  | "LightSalmon"
  | "Orchid"
  | "Thistle"
  | "DarkKhaki"
  | "LightCoral"
  | "MediumOrchid"
  | "None";

export const Colors: Record<ColorName, ColorString> = {
  White: "hsl(0, 0%, 100%)",
  Black: "hsl(0, 0%, 0%)",
  Tomato: "hsl(9, 100%, 64%)",
  LightBlue: "hsl(195, 53%, 79%)",
  LightGreen: "hsl(120, 73%, 75%)",
  Purple: "hsl(300, 47%, 75%)",
  Orange: "hsl(39, 100%, 50%)",
  Yellow: "hsl(51, 100%, 50%)",
  Cyan: "hsl(180, 100%, 50%)",
  Teal: "hsl(180, 100%, 25%)",
  Magenta: "hsl(300, 100%, 50%)",
  Brown: "hsl(25, 76%, 31%)",
  Lime: "hsl(120, 100%, 50%)",
  Coral: "hsl(16, 100%, 66%)",
  Turquoise: "hsl(174, 72%, 56%)",
  Indigo: "hsl(275, 100%, 25%)",
  Violet: "hsl(300, 76%, 72%)",
  Lavender: "hsl(240, 67%, 94%)",
  Peach: "hsl(28, 100%, 86%)",
  SkyBlue: "hsl(197, 71%, 73%)",
  Olive: "hsl(60, 100%, 25%)",
  Tan: "hsl(34, 44%, 69%)",
  Salmon: "hsl(6, 93%, 71%)",
  Maroon: "hsl(0, 100%, 25%)",
  Navy: "hsl(240, 100%, 25%)",
  Khaki: "hsl(54, 77%, 75%)",
  Periwinkle: "hsl(54, 77%, 75%)",
  Mint: "hsl(150, 100%, 98%)",
  Azure: "hsl(180, 100%, 97%)",
  Chartreuse: "hsl(90, 100%, 50%)",
  Goldrod: "hsl(43, 74%, 49%)",
  SlateBlue: "hsl(248, 53%, 58%)",
  LightSeaGreen: "hsl(177, 70%, 41%)",
  DarkCyan: "hsl(180, 100%, 27%)",
  RosyBrown: "hsl(0, 25%, 65%)",
  PaleVioletRed: "hsl(340, 60%, 65%)",
  DeepPink: "hsl(328, 100%, 54%)",
  DarkOrange: "hsl(33, 100%, 50%)",
  Crimson: "hsl(348, 83%, 47%)",
  LightSalmon: "hsl(17, 100%, 74%)",
  Orchid: "hsl(302, 59%, 65%)",
  Thistle: "hsl(300, 24%, 80%)",
  DarkKhaki: "hsl(56, 38%, 58%)",
  LightCoral: "hsl(0, 79%, 72%)",
  MediumOrchid: "hsl(288, 59%, 58%)",
  None: "hsla(0, 0%, 0%, 0)",
} as const;
