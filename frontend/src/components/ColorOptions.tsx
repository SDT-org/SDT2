import React from "react";
import { ColorOption } from "../colors";
import { formatTitle } from "../helpers";

export const ColorOptions = () =>
  Object.entries(ColorOption).map(([key, value]) => (
    <option key={key} value={value}>
      {formatTitle(key)}
    </option>
  ));
