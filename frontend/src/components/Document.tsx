import type React from "react";
import { type DocState, useAppState } from "../appState";
import { useDocState } from "../hooks/useDocState";
import { Loader } from "./Loader";
import { MainMenu } from "./Menu";
import { Runner } from "./Runner";
import { Select, SelectItem } from "./Select";
import { Viewer } from "./Viewer";

export const Document = ({
  id,
}: {
  id: string;
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
    mainMenu: <MainMenu />,
  };

  const VIEWS: { [K in DocState["view"]]: React.ReactElement } = {
    runner: <Runner {...commonViewProps} />,
    loader: <Loader {...commonViewProps} />,
    viewer: <Viewer {...commonViewProps} />,
  };

  return VIEWS[docState?.view || "viewer"];
};
