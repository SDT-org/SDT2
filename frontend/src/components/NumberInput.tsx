import React from 'react';

const normalizeNumber = (
  value: any,
  type: "int" | "float",
  validMin: number,
  validMax: number,
  step?: number,
) => {
  try {
    value = type === "float" ? parseFloat(value) : parseInt(value, 10);
    if (isNaN(value)) {
      return undefined;
    }
  } catch {
    return undefined;
  }

  if (step) {
    value = Math.round(value / step) * step;
  }

  return Math.min(Math.max(value, validMin), validMax);
};

export const NumberInput = ({
  field,
  value,
  updateValue,
  label,
  min,
  max,
  step,
  type = "int",
  ...inputProps
}: {
  field: string;
  value: number;
  updateValue: ({}) => void;
  label: string;
  min: number;
  max: number;
  step?: number;
  type?: "int" | "float";
} & React.InputHTMLAttributes<HTMLInputElement>) => {
  const [displayValue, setDisplayValue] = React.useState(value?.toString());
  const [invalidValue, setInvalidValue] = React.useState(false);
  return (
    <div className="field" data-invalid-value={invalidValue}>
      <label htmlFor={field}>{label}</label>
      <input
        type="number"
        name={field}
        id={field}
        min={min}
        max={max}
        step={step}
        value={displayValue}
        onChange={(e) => {
          const newValue = normalizeNumber(e.target.value, type, min, max);

          setDisplayValue(newValue?.toString() || "");

          if (newValue !== undefined) {
            setInvalidValue(false);
            updateValue({
              [field]: newValue,
            });
          } else {
            setInvalidValue(true);
          }
        }}
        onBlur={() => {
          setDisplayValue(value.toString());
          setInvalidValue(false);
        }}
        {...inputProps}
      />
    </div>
  );
};
