import type { SaveableImageFormat } from "../appState";
import type { RunProcessDataArgs } from "../components/Runner";
import type { AppState } from "../src/appState";

declare global {
  interface Window {
    pywebview: {
      api: {
        app_config: () => Promise<string>;
        get_state: () => Promise<AppState>;
        reset_state: () => Promise<void>;
        open_file_dialog: () => Promise<void>;
        select_alignment_output_path: () => Promise<void>;
        select_export_path: () => Promise<void>;
        export_data: (args: {
          output_cluster: boolean;
          cluster_threshold_one: number;
          cluster_threshold_two: number;
          heatmap_image_data: string;
          distribution_image_data: string;
          image_format: SaveableImageFormat;
        }) => Promise<boolean>;
        save_image: (args: {
          data: string;
          format: AppState["client"]["saveFormat"];
        }) => Promise<void>;
        run_process_data: (args: RunProcessDataArgs) => Promise<void>;
        processes_info: () => Promise<any>;
        get_available_memory: () => Promise<number>;
        cancel_run: () => Promise<void>;
        get_data: () => Promise<string>;
        show_about: () => Promise<void>;
        show_manual: () => Promise<void>;
        close_app: () => Promise<void>;
      };
    };
    syncAppState: (state: AppState) => void;
  }
}
