import React from "react";
import { Colors } from "../colors";
import { formatTitle } from "../helpers";

export const ColorOptions = () =>
  Object.entries(Colors).map(([key, value]) => (
    <option key={key} value={value}>
      {formatTitle(key)}
    </option>
  ));
