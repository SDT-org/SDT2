import type {
  DocState,
  SaveableImageFormat,
  SyncProgressEvent,
  SyncStateEvent,
} from "../appState";
import type { RunProcessDataArgs } from "../components/Runner";
import { DistributionState } from "../distributionState";
import { type GetDataResponse, HeatmapSettings } from "../plotTypes";
import type { AppState } from "../src/appState";

declare global {
  interface Window {
    pywebview: {
      api: {
        app_config: () => Promise<AppState["config"]>;
        app_settings: () => Promise<{ recent_files: string[] }>;
        get_state: () => Promise<AppState>;
        reset_state: () => Promise<void>;
        new_doc: () => Promise<string>;
        get_doc: (doc_id: string) => Promise<DocState>;
        save_doc: (
          doc_id: string,
          path: string,
          save_as: boolean,
        ) => Promise<boolean>;
        close_doc: (doc_id: string) => Promise<void>;
        save_doc_settings: (docState: DocState) => Promise<void>;
        open_file: (
          filepath: string,
          doc_id?: string,
        ) => Promise<[string, string]>;
        open_file_dialog: (
          doc_id?: string,
          last_data_filepath?: string,
        ) => Promise<[string, string]> | Promise<void>;
        save_file_dialog: (
          filename: string,
          defaultDirectory?: string,
        ) => Promise<string>;
        select_path_dialog: (defaultDirectory?: string) => Promise<string>;
        generate_cluster_data: (
          doc_id: string,
          threshold_one: number,
          threshold_two: number,
        ) => Promise<{ id: string; group: number }[]>;
        export_data: (args: {
          doc_id: string;
          export_path: string;
          output_cluster: boolean;
          cluster_threshold_one: number;
          cluster_threshold_two: number;
          heatmap_image_data: string;
          clustermap_image_data: string;
          histogram_image_data: string;
          violin_image_data: string;
          image_format: SaveableImageFormat;
        }) => Promise<boolean>;
        start_run: (args: RunProcessDataArgs) => Promise<void>;
        set_window_title: (title: string) => Promise<void>;
        processes_info: () => Promise<string>;
        get_available_memory: () => Promise<number>;
        cancel_run: (
          doc_id: string,
          run_settings?: "preserve" | "clear",
        ) => Promise<void>;
        get_data: (docId: string) => Promise<string>;
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
    "sync-progress": SyncProgressEvent;
  }
}
