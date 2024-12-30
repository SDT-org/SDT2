import type { DocState } from "../appState";
import { saveFileDialog } from "./files";

export const saveDocument = async (docState: DocState, saveAs?: boolean) => {
  let filename = docState.filename;

  if (saveAs) {
    const [success, result] = await saveFileDialog(
      `${docState.basename.split(".").slice(0, -1).join(".")}.sdt`,
    );
    if (success) {
      filename = result;
    } else {
      return false;
    }
  }

  return window.pywebview.api
    .save_doc_settings(docState)
    .then(() => window.pywebview.api.save_doc(docState.id, filename));
};

export const closeDocument = (docId: string) =>
  window.pywebview.api.close_doc(docId);

export const newDocument = () => window.pywebview.api.new_doc();
export const getDocument = (docId: string) =>
  window.pywebview.api.get_doc(docId);
