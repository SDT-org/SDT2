// https://react-spectrum.adobe.com/react-aria/Checkbox.html
import type { CheckboxProps } from "react-aria-components";
import { Checkbox as RACCheckbox } from "react-aria-components";

export const Checkbox = ({
  children,
  ...props
}: Omit<CheckboxProps, "children"> & {
  children?: React.ReactNode;
}) => (
  <RACCheckbox {...props}>
    <div className="checkbox">
      <svg viewBox="0 0 18 18" aria-hidden="true">
        <polyline points="1 9 7 14 15 4" />
      </svg>
    </div>
    {children}
  </RACCheckbox>
);
