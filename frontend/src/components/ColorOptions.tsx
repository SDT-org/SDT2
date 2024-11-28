import { ColorOption } from "../colors";
import { formatTitle } from "../helpers";
import React from "react";

export const ColorOptions = () =>
  Object.entries(ColorOption).map(([key, value]) => (
    <option key={key} value={value}>
      {formatTitle(key)}
    </option>
  ));
