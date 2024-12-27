import type React from "react";
import { type DocState, useAppState } from "../appState";
import { useDocState } from "../hooks/useDocState";
import { Loader } from "./Loader";
import { Runner } from "./Runner";
import { Viewer } from "./Viewer";

export const Document = ({
  id,
  tabView,
  leftSidebarCollapsed,
}: {
  id: string;
  leftSidebarCollapsed: boolean;
  tabView: "tabs" | "select";
}) => {
  const { appState, setAppState } = useAppState();
  const { docState, setDocState, updateDocState } = useDocState(
    id,
    appState,
    setAppState,
  );
  const commonViewProps = {
    docState,
    setDocState,
    updateDocState,
    tabView,
    leftSidebarCollapsed,
  };

  const VIEWS: { [K in DocState["view"]]: React.ReactElement } = {
    runner: <Runner {...commonViewProps} />,
    loader: <Loader {...commonViewProps} />,
    viewer: <Viewer {...commonViewProps} />,
  };

  return VIEWS[docState?.view || "viewer"];
};
