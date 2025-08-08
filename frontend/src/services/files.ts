import messages from "../messages";

const handleOpenFileFailure = (e: unknown): [boolean, [string, string]] => {
  console.error(e);
  alert(
    e &&
      typeof e === "object" &&
      "message" in e &&
      typeof e.message === "string" &&
      Object.keys(messages).includes(e.message)
      ? `Error: ${messages[e.message as keyof typeof messages]}`
      : "An error occured while loading this file. Please ensure it is a valid, SDT-compatible file.",
  );
  return [false, ["", ""]];
};

export const openFile = async (path: string, docId?: string) => {
  try {
    const data = await window.pywebview.api.files.open_file(path, docId);
    return [!!data, data];
  } catch (e) {
    return handleOpenFileFailure(e);
  }
};

export const openFileDialog = async (
  docId?: string,
  path?: string,
): Promise<[boolean, [string, string]]> => {
  try {
    const data = await window.pywebview.api.files.open_file_dialog(docId, path);
    return [!!data, data || ["", ""]];
  } catch (e) {
    return handleOpenFileFailure(e);
  }
};

export const saveFileDialog = async (
  filename: string,
): Promise<[boolean, string]> => {
  try {
    const data = await window.pywebview.api.files.save_file_dialog(filename);
    return [!!data, data];
  } catch (e) {
    console.error(e);
    alert("An error occured while saving this file. Please try again.");
    return [false, ""];
  }
};
