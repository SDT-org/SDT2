import type { SetAppState } from "../appState";
import { services } from "../services";

const useNewDocument = (setAppState: SetAppState) => () =>
  services
    .newDocument()
    .then((id: string) =>
      setAppState((prev) => ({ ...prev, activeDocumentId: id })),
    );

export default useNewDocument;
