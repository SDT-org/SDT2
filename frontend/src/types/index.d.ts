import { RunProcessDataArgs } from "../components/Runner";
import { AppState } from "../src/appState";

declare global {
  interface Window {
    pywebview: {
      api: {
        get_state: () => Promise<AppState>;
        reset_state: () => Promise<void>;
        open_file_dialog: () => Promise<void>;
        select_alignment_output_path: () => Promise<void>;
        select_export_path: () => Promise<void>;
        export_data: ({
          output_cluster: boolean,
          cluster_threshold_one: number,
          cluster_threshold_two: number,
        }) => Promise<boolean>;
        save_image: (args: {
          data: string;
          format: AppState["client"]["saveFormat"];
        }) => Promise<void>;
        run_process_data: (args: RunProcessDataArgs) => Promise<void>;
        cancel_run: () => Promise<void>;
        get_heatmap_data: () => Promise<string>;
        get_line_histo_data: () => Promise<string>;
      };
    };
    syncAppState: () => Promise<void>;
  }
}
