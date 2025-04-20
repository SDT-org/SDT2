declare module "d3-color-legend" {
  import type * as d3 from "d3";

  export interface LegendOptions {
    color: d3.ScaleLinear<string, string>;
    title?: string;
    width?: number;
    height?: number;
    marginTop?: number;
    marginRight?: number;
    marginBottom?: number;
    marginLeft?: number;
    ticks?: number;
    tickFormat?: (d: number) => string;
    tickValues?: number[];
    tickAdjust?: (
      selection: d3.Selection<SVGGElement, unknown, null, undefined>,
    ) => void;
    scheme?: string;
    [key: string]: unknown;
  }

  export function Legend(
    selection: d3.Selection<SVGGElement, unknown, null, undefined>,
    options: LegendOptions,
  ): void;
}
