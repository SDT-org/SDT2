import React from "react";
import { Switch as RASwitch, SwitchProps } from "react-aria-components";

export const Switch = ({
  children,
  ...props
}: React.PropsWithChildren<SwitchProps>) => {
  return (
    <>
      <RASwitch {...props}>
        <div className="indicator" />
        {children}
      </RASwitch>
    </>
  );
};
