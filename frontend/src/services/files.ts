export const openFileDialog = async (
  docId?: string,
  path?: string,
): Promise<[boolean, [string, string]]> => {
  try {
    const data = await window.pywebview.api.open_file_dialog(docId, path);
    return [!!data, data || ["", ""]];
  } catch (e) {
    console.error(e);
    alert(
      "An error occured while loading this file. Please ensure it is a valid, SDT-compatible file.",
    );
    return [false, ["", ""]];
  }
};

export const saveFileDialog = async (
  filename: string,
): Promise<[boolean, string]> => {
  try {
    const data = await window.pywebview.api.save_file_dialog(filename);
    return [!!data, data];
  } catch (e) {
    console.error(e);
    alert("An error occured while saving this file. Please try again.");
    return [false, ""];
  }
};
