import React from "react";
import { type DocState, useAppState } from "../appState";
import { useDocState } from "../hooks/useDocState";
import { ErrorPanel } from "./ErrorPanel";
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

  React.useEffect(() => {
    if (id !== appState.activeDocumentId) {
      return;
    }
    const handleKeyDown = (event: KeyboardEvent) => {
      if (
        event.shiftKey &&
        (event.metaKey || event.altKey) &&
        event.key === "i"
      ) {
        event.preventDefault();
        window.pywebview.api.documents.open_doc_folder(id);
      }
    };

    document.addEventListener("keydown", handleKeyDown);
    return () => document.removeEventListener("keydown", handleKeyDown);
  }, [id, appState.activeDocumentId]);

  if (docState.invalid) {
    return <ErrorPanel docState={docState} />;
  }

  return VIEWS[docState?.view || "viewer"];
};
