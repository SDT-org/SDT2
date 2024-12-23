import type { DocState } from "../appState";
import { saveFileDialog } from "./files";

export const saveDocument = async (docState: DocState) => {
  if (!docState.savepath) {
    return saveFileDialog(
      `${docState.basename.split(".").slice(0, -1).join(".")}.sdt`,
    ).then(([success, result]) => {
      console.log(success, result);
      if (success) {
        return window.pywebview.api.save_doc(docState.id, result);
      }
      return false;
    });
  }
  return window.pywebview.api.save_doc(docState.id, docState.savepath);
};

export const closeDocument = (docId: string) =>
  window.pywebview.api.close_doc(docId);
