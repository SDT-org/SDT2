import type { SaveableImageFormat, SyncStateEvent } from "../appState";
import type { RunProcessDataArgs } from "../components/Runner";
import type { AppState } from "../src/appState";

declare global {
  interface Window {
    pywebview: {
      api: {
        app_config: () => Promise<string>;
        get_state: () => Promise<AppState>;
        reset_state: () => Promise<void>;
        open_file_dialog: (lastDataFilepath?: string) => Promise<string>;
        select_path_dialog: (defaultDirectory?: string) => Promise<string>;
        export_data: (args: {
          export_path: string;
          output_cluster: boolean;
          cluster_threshold_one: number;
          cluster_threshold_two: number;
          heatmap_image_data: string;
          histogram_image_data: string;
          violin_image_data: string;
          raincloud_image_data: string;
          image_format: SaveableImageFormat;
        }) => Promise<boolean>;
        save_image: (args: {
          export_path: string;
          data: string;
          format: AppState["client"]["saveFormat"];
        }) => Promise<void>;
        run_process_data: (args: RunProcessDataArgs) => Promise<void>;
        processes_info: () => Promise<string>;
        get_available_memory: () => Promise<number>;
        cancel_run: () => Promise<void>;
        get_data: () => Promise<string>;
        show_about: () => Promise<void>;
        show_manual: () => Promise<void>;
        close_app: () => Promise<void>;
      };
    };
    syncAppState: (state: AppState) => void;
    APP_STATE: AppState;
    LAST_ERROR: {
      error: Error;
      errorInfo: React.ErrorInfo;
    };
  }

  interface GlobalEventHandlersEventMap {
    "sync-state": SyncStateEvent;
  }
}
