import React from "react";
import {
  Button,
  Group,
  Input,
  Label,
  NumberField,
  NumberFieldProps,
} from "react-aria-components";

export const NumberInput = ({
  field,
  updateValue,
  label,
  min,
  max,
  step,
  type = "int",
  ...props
}: {
  label: string;
  min: number;
  max: number;
  step?: number;
  type?: "int" | "float";
} & (
  | { updateValue: ({}) => void; field: string }
  | { updateValue?: ({}) => void; field?: string }
) &
  NumberFieldProps) => {
  return (
    <div className="field">
      <NumberField
        minValue={min}
        maxValue={max}
        onChange={(value) => updateValue?.({ [field || ""]: value })}
        formatOptions={
          type === "float"
            ? {
                minimumFractionDigits: 1,
                maximumFractionDigits: 2,
              }
            : {
                minimumFractionDigits: 0,
                maximumFractionDigits: 0,
              }
        }
        {...props}
      >
        <Label>{label}</Label>
        <Group className="react-aria-Group number-input-group">
          <Button slot="decrement">
            <svg
              height="8"
              width="8"
              viewBox="0 0 8 8"
              xmlns="http://www.w3.org/2000/svg"
            >
              <path
                d="M 1 3 L 4 6 L 7 3"
                fill="none"
                stroke="currentColor"
                stroke-width="2"
                stroke-linecap="round"
                stroke-linejoin="round"
              />
            </svg>
          </Button>
          <Input />
          <Button slot="increment">
            <svg
              height="8"
              width="8"
              viewBox="0 0 8 8"
              xmlns="http://www.w3.org/2000/svg"
            >
              <path
                d="M 1 5 L 4 2 L 7 5"
                fill="none"
                stroke="currentColor"
                stroke-width="2"
                stroke-linecap="round"
                stroke-linejoin="round"
              />
            </svg>
          </Button>
        </Group>
      </NumberField>
    </div>
  );
};
