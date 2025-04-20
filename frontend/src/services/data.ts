export const getData = (docId: string) =>
  window.pywebview.api.get_data(docId).catch((e) => {
    console.error(e);
    alert(
      "An error occured while processing this file. Please ensure it is a valid, SDT-compatible file.",
    );
    window.pywebview.api.reset_state();
    return "";
  });
