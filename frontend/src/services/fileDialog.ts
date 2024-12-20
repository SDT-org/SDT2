export const openFileDialog = async (
  path?: string,
): Promise<[boolean, string]> => {
  try {
    const data = await window.pywebview.api.open_file_dialog(path);
    return [!!data, data || ""];
  } catch (e) {
    console.error(e);
    alert(
      "An error occured while loading this file. Please ensure it is a valid, SDT-compatible file.",
    );
    return [false, ""];
  }
};
