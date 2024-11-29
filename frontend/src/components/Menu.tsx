import React from "react";
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

export type MainMenuProps = {
  appState: AppState;
  onNew: () => void;
  onOpen: () => void;
  onExport: () => void;
  onManual: () => void;
  onAbout: () => void;
  onExit: () => void;
};

export const MainMenu = ({
  appState,
  onNew,
  onOpen,
  onExport,
  onManual,
  onAbout,
  onExit,
}: MainMenuProps) => (
  <AppMenuButton label="☰">
    <AppMenuItem onAction={onNew}>New</AppMenuItem>
    <AppMenuItem onAction={onOpen}>Open...</AppMenuItem>
    {appState.view === "viewer" ? (
      <AppMenuItem onAction={onExport}>Export images and data...</AppMenuItem>
    ) : null}
    <Separator />
    <AppMenuItem onAction={onManual}>Manual</AppMenuItem>
    <AppMenuItem onAction={onAbout}>About</AppMenuItem>
    <Separator />
    <AppMenuItem onAction={onExit}>Exit</AppMenuItem>
  </AppMenuButton>
);
