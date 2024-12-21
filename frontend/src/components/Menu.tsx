import { createHideableComponent } from "@react-aria/collections";
import {
  Button,
  Menu,
  MenuItem,
  type MenuItemProps,
  type MenuProps,
  MenuTrigger,
  type MenuTriggerProps,
  Popover,
  Separator,
} from "react-aria-components";
import { type AppState, useAppState } from "../appState";
import useOpenFileDialog from "../hooks/useOpenFileDialog";

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
      <Button
        className="react-aria-Button main-menu-button"
        aria-label="Application Menu"
      >
        <svg
          height="16"
          width="16"
          aria-hidden="true"
          enableBackground="new 0 0 24 24"
          viewBox="0 0 24 24"
          xmlns="http://www.w3.org/2000/svg"
          color="currentcolor"
        >
          <g
            style={{
              fill: "none",
              stroke: "currentcolor",
              strokeWidth: 2,
              strokeLinecap: "round",
              strokeLinejoin: "round",
              strokeMiterlimit: "10",
            }}
          >
            <path d="m3 5h18" />
            <path d="m3 12h18" />
            <path d="m3 19h18" />
            <path d="m3 5h18" />
            <path d="m3 12h18" />
            <path d="m3 19h18" />
          </g>
        </svg>
      </Button>
      <Popover>
        <Menu {...props}>{children}</Menu>
      </Popover>
    </MenuTrigger>
  );
};

const AppMenuItem = (props: MenuItemProps) => {
  const textValue =
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
            <svg aria-hidden="true" className="chevron" viewBox="0 0 24 24">
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

export const MainMenu = createHideableComponent(() => {
  const { appState, setAppState } = useAppState();
  const openFileDialog = useOpenFileDialog(appState, setAppState);
  const onNew = () => {
    if (
      !(
        appState.view !== "runner" &&
        confirm("Are you sure? Current results will be cleared.")
      )
    ) {
      return;
    }

    window.pywebview.api.reset_state();
  };

  const onOpen = () => {
    if (
      appState.view !== "runner" &&
      !confirm("Are you sure? Current results will be cleared.")
    ) {
      return;
    }

    openFileDialog();
  };

  const onExport = () =>
    setAppState((previous) => ({
      ...previous,
      client: { ...previous.client, showExportModal: true },
    }));

  const onAbout = () => window.pywebview.api.show_about();
  const onManual = () => window.pywebview.api.show_manual();
  const onExit = () => {
    if (confirm("Are you sure you want to exit?")) {
      window.pywebview.api.close_app();
    }
  };

  return (
    <AppMenuButton label="☰">
      <AppMenuItem onAction={onNew}>New</AppMenuItem>
      <AppMenuItem onAction={onOpen}>Open...</AppMenuItem>
      {appState.view === "viewer" ? (
        <AppMenuItem onAction={onExport}>Export images and data...</AppMenuItem>
      ) : null}
      <AppMenuItem onAction={onManual}>Manual</AppMenuItem>
      <AppMenuItem onAction={onAbout}>About</AppMenuItem>
      <Separator />
      <AppMenuItem onAction={onExit}>Exit</AppMenuItem>
    </AppMenuButton>
  );
});
