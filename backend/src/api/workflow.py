from multiprocessing import Manager, cpu_count as get_cpu_count
from document_paths import build_document_paths
from workflow.models import LzaniSettings, ParasailSettings, RunSettings, WorkflowRun
from workflow.runner import run_process
from app_state import get_state, set_state, get_document, update_document
from app_globals import get_workflow_runs, get_parsed_workflow_results, get_canceled, set_canceled


class Workflow:
    def start_workflow_run(self, args: dict):
        from app import get_lzani_exec_path
        
        app_state = get_state()
        if app_state.debug:
            print("\nAPI args:", args)

        if app_state.active_run_document_id:
            raise Exception("Multiple runs are not supported")

        doc = get_document(args["doc_id"])

        parsed_result = get_parsed_workflow_results().get(doc.id)
        if not parsed_result:
            raise Exception(
                "Parsed workflow result not found. Close the document and reselect the file."
            )

        if parsed_result.error:
            raise Exception("Workflow has an error that must be resolved.")

        settings = RunSettings(
            fasta_path=doc.filename,
            doc_paths=build_document_paths(doc.tempdir_path),
            output_path=doc.tempdir_path,
            cluster_method=args.get("cluster_method", ""),
            analysis_method=args.get("analysisMethod", "parasail"),
            lzani=LzaniSettings(
                exec_path=get_lzani_exec_path(),
                score_type=args.get("lzani_score_type", "ani"),
                aw=args.get("aw"),
                am=args.get("am"),
                mal=args.get("mal"),
                msl=args.get("msl"),
                mrd=args.get("mrd"),
                mqd=args.get("mqd"),
                reg=args.get("reg"),
                ar=args.get("ar"),
            ),
            parasail=ParasailSettings(
                process_count=max(
                    min(args.get("compute_cores", 1), get_cpu_count()), 1
                ),
                scoring_matrix=args.get("scoring_matrix"),
                open_penalty=args.get("open_penalty"),
                extend_penalty=args.get("extend_penalty"),
            ),
            export_alignments=args.get("export_alignments", False),
            alignment_export_path=args.get("alignment_export_path", ""),
        )

        workflow_run = WorkflowRun(
            result=parsed_result,
            settings=settings,
            progress=0,
        )

        get_workflow_runs()[doc.id] = workflow_run

        update_document(doc.id, view="loader")
        with Manager() as manager:
            canceled = manager.Value("b", False)
            set_canceled(canceled)
            result = run_process(workflow_run, canceled)

            if result.error:
                if result.error == "PROCESS_CANCELED":
                    return
                raise Exception(f"Workflow processing step failed: {result.error}")

        set_state(active_run_document_id=None)
        workflow_runs = get_workflow_runs()
        if doc.id in workflow_runs:
            del workflow_runs[doc.id]
        update_document(doc.id, view="viewer")

    def get_workflow_run_status(self, doc_id: str):
        workflow_run = get_workflow_runs().get(doc_id)
        if not workflow_run:
            raise Exception(f"No workflow run found for document ID: {doc_id}")
        return {"stage": workflow_run.stage, "progress": workflow_run.progress}

    def cancel_run(self, doc_id: str, run_settings: str):
        from app import do_cancel_run
        
        do_cancel_run()
        if run_settings == "preserve":
            update_document(
                doc_id,
                view="runner",
                progress=0,
                pair_progress=0,
                pair_count=0,
                estimated_time=None,
                stage="",
                sequences_count=0,
            )
        elif run_settings == "clear":
            from api.documents import Documents
            Documents().reset_state()