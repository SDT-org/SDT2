import type {
  DocState,
  SaveableImageFormat,
  SaveableImageKey,
  SyncProgressEvent,
  SyncStateEvent,
} from "../appState";
import type { RunSettings } from "../components/Runner";
import { DistributionState } from "../distributionState";
import type {
  ClusterDataItem,
  GetClustermapDataResponse,
  GetDataResponse,
  GetUMAPDataResponse,
  HeatmapData,
  HeatmapSettings,
} from "../plotTypes";
import type { AppState } from "../src/appState";

declare global {
  interface Window {
    pywebview: {
      api: {
        app_config: () => Promise<AppState["config"]>;
        app_settings: () => Promise<{
          recent_files: string[];
          export_path: string;
        }>;
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
        get_clustermap_data: (
          doc_id: string,
          threshold: number,
          method: string,
        ) => Promise<GetClustermapDataResponse>;
        get_umap_data: (
          doc_id: string,
          params: {
            n_neighbors: number;
            min_dist: number;
            threshold: number;
            methods: string[];
          },
        ) => Promise<GetUMAPDataResponse>;
        export: (args: {
          doc_id: string;
          export_path: string;
          output_cluster: boolean;
          cluster_threshold: number;
          cluster_method: string;
          image_format: SaveableImageFormat;
          open_folder: boolean;
          prefix: string;
        }) => Promise<boolean>;
        save_svg_element: (
          doc_id: string,
          selector: string,
          key: SaveableImageKey,
        ) => Promise<void>;
        save_svg_data: (
          doc_id: string,
          data: string,
          key: SaveableImageKey,
          format: SaveableImageFormat,
        ) => Promise<void>;
        save_raster_image: (
          doc_id: string,
          data: string,
          key: SaveableImageKey,
          format: SaveableImageFormat,
        ) => Promise<void>;
        start_workflow_run: (args: RunSettings) => Promise<void>;
        get_workflow_run_status: (
          doc_id: string,
        ) => Promise<{ stage: string; progress?: number }>;
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
        open_doc_folder: (doc_id: string) => Promise<void>;
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
