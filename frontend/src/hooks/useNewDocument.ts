import type { SetAppState } from "../appState";
import { services } from "../services";

export const useNewDocument = (setAppState: SetAppState) => () =>
  services
    .newDocument()
    .then((id: string) =>
      setAppState((prev) => ({ ...prev, activeDocumentId: id })),
    );
