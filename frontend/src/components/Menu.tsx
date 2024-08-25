import React, { SVGProps } from "react";
import {
  Menu,
  MenuTrigger,
  type MenuItemProps,
  type MenuProps,
  type MenuTriggerProps,
  MenuItem,
  Button,
  Popover,
  Separator,
} from "react-aria-components";
import { AppState } from "../appState";

// AppMenuButton and AppMenuItem were derived from https://react-spectrum.adobe.com/react-aria/Menu.html#reusable-wrappers
interface MyMenuButtonProps<T>
  extends MenuProps<T>,
    Omit<MenuTriggerProps, "children"> {
  label?: string;
}

const AppMenuButton = <T extends object>({
  label,
  children,
  ...props
}: MyMenuButtonProps<T>) => {
  return (
    <MenuTrigger {...props}>
      <Button>{label}</Button>
      <Popover>
        <Menu {...props}>{children}</Menu>
      </Popover>
    </MenuTrigger>
  );
};

const AppMenuItem = (props: MenuItemProps) => {
  let textValue =
    props.textValue ||
    (typeof props.children === "string" ? props.children : "");
  return (
    <MenuItem
      {...props}
      textValue={textValue}
      className={({ isFocused, isOpen }) =>
        `app-item ${isFocused ? "focused" : ""} ${isOpen ? "open" : ""}`
      }
    >
      {({ hasSubmenu }) => (
        <>
          {props.children}
          {hasSubmenu && (
            <svg className="chevron" viewBox="0 0 24 24">
              <path d="m9 18 6-6-6-6" />
            </svg>
          )}
        </>
      )}
    </MenuItem>
  );
};

const Hamburger = () => (
  <svg
    xmlns="http://www.w3.org/2000/svg"
    xmlSpace="preserve"
    viewBox="0 0 122.88 95.95"
  >
    <path
      d="M8.94 0h105c4.92 0 8.94 4.02 8.94 8.94s-4.02 8.94-8.94 8.94h-105C4.02 17.88 0 13.86 0 8.94S4.02 0 8.94 0zm0 78.07h105c4.92 0 8.94 4.02 8.94 8.94s-4.02 8.94-8.94 8.94h-105C4.02 95.95 0 91.93 0 87.01s4.02-8.94 8.94-8.94zm0-39.04h105c4.92 0 8.94 4.02 8.94 8.94s-4.02 8.94-8.94 8.94h-105C4.02 56.91 0 52.89 0 47.97c0-4.91 4.02-8.94 8.94-8.94z"
      style={{
        fillRule: "evenodd",
        clipRule: "evenodd",
      }}
    />
  </svg>
);

export type MainMenuProps = {
  appState: AppState;
  onNew: () => void;
  onExport: () => void;
  onSettings: () => void;
  onManual: () => void;
  onAbout: () => void;
  onExit: () => void;
};

export const MainMenu = ({
  appState,
  onNew,
  onExport,
  onSettings,
  onManual,
  onAbout,
  onExit,
}: MainMenuProps) => (
  <AppMenuButton label="☰ Menu">
    <AppMenuItem onAction={onNew}>New...</AppMenuItem>
    {appState.view === "viewer" ? (
      <AppMenuItem onAction={onExport}>Export images and data...</AppMenuItem>
    ) : null}
    <Separator />
    <AppMenuItem onAction={onSettings}>Settings</AppMenuItem>
    <AppMenuItem onAction={onManual}>Manual</AppMenuItem>
    <AppMenuItem onAction={onAbout}>About</AppMenuItem>
    <Separator />
    <AppMenuItem onAction={onExit}>Exit</AppMenuItem>
  </AppMenuButton>
);
