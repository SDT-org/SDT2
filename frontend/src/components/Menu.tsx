import { createHideableComponent } from "@react-aria/collections";
import React from "react";
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
  SubmenuTrigger,
} from "react-aria-components";
import { TbMenu2 } from "react-icons/tb";
import { type AppState, findDoc, useAppState } from "../appState";
import { isSDTFile } from "../helpers";
import { useCloseDocument } from "../hooks/useCloseDocument";
import useNewDocument from "../hooks/useNewDocument";
import useOpenFileDialog from "../hooks/useOpenFileDialog";
import { useRecentFiles } from "../hooks/useRecentFiles";
import { useSaveActiveDocument } from "../hooks/useSaveActiveDocument";

// AppMenuButton and AppMenuItem were derived from https://react-spectrum.adobe.com/react-aria/Menu.html#reusable-wrappers
interface MyMenuButtonProps<T>
  extends MenuProps<T>,
    Omit<MenuTriggerProps, "children"> {
  label?: string;
}

const AppMenuButton = <T extends object>({
  children,
  ...props
}: MyMenuButtonProps<T>) => {
  return (
    <MenuTrigger {...props}>
      <Button
        className="react-aria-Button main-menu-button"
        aria-label="Application Menu"
      >
        <TbMenu2 size={18} />
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
  const activeDocState = React.useMemo(
    () => findDoc(appState.activeDocumentId, appState.documents),
    [appState],
  );

  const newDocument = useNewDocument(setAppState);
  const closeDocument = useCloseDocument(appState, setAppState);
  const saveDocument = useSaveActiveDocument(appState, setAppState);

  const openRecentFile = useRecentFiles(appState, setAppState);
  const sdtFile = activeDocState && isSDTFile(activeDocState.filetype);

  const openFileDialog = useOpenFileDialog(appState, setAppState);
  const onNew = newDocument;
  const onOpen = () => {
    openFileDialog();
  };
  const onSave = () => saveDocument(false);
  const onSaveAs = () => saveDocument(true);
  const onClose = () => closeDocument(appState.activeDocumentId);

  const onExport = () =>
    setAppState((previous) => ({
      ...previous,
      showExportModal: true,
    }));

  const onAbout = () => window.pywebview.api.show_about();
  const onManual = () => window.pywebview.api.show_manual();
  const onExit = () => {
    if (confirm("Are you sure you want to exit?")) {
      window.pywebview.api.close_app();
    }
  };

  return (
    <AppMenuButton>
      <AppMenuItem onAction={onNew}>New</AppMenuItem>
      <AppMenuItem onAction={onOpen}>Open...</AppMenuItem>
      {activeDocState ? (
        <SubmenuTrigger>
          <AppMenuItem>Open Recent</AppMenuItem>
          <Popover>
            <Menu>
              {appState.recentFiles.map((filePath) => (
                <AppMenuItem
                  key={filePath}
                  onAction={() => openRecentFile(filePath, activeDocState)}
                >
                  {filePath.split(/(\/|\\)/).pop()}
                </AppMenuItem>
              ))}
            </Menu>
          </Popover>
        </SubmenuTrigger>
      ) : null}
      <Separator />
      {activeDocState?.view === "viewer" ? (
        <>
          {
            <AppMenuItem onAction={onSave}>
              Save{!sdtFile ? "..." : null}
            </AppMenuItem>
          }
          {sdtFile ? (
            <AppMenuItem onAction={onSaveAs}>Save As...</AppMenuItem>
          ) : null}
          <AppMenuItem
            isDisabled={activeDocState?.view !== "viewer"}
            onAction={onExport}
          >
            Export...
          </AppMenuItem>
        </>
      ) : null}
      <AppMenuItem onAction={onClose}>Close</AppMenuItem>
      <Separator />
      <AppMenuItem onAction={onManual}>Manual</AppMenuItem>
      <AppMenuItem onAction={onAbout}>About</AppMenuItem>
      <Separator />
      <AppMenuItem onAction={onExit}>Exit</AppMenuItem>
    </AppMenuButton>
  );
});
