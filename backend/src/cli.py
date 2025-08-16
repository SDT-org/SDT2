#!/usr/bin/env python3

import argparse
from multiprocessing import Manager
import os
import sys
import tempfile

current_file_path = os.path.abspath(os.path.join(os.path.dirname(__file__)))
sys.path.append(os.path.join(current_file_path, "../../"))
sys.path.append(os.path.join(current_file_path, "."))

from mime_setup import register_mimetypes
from document_paths import build_document_paths
from workflow.models import LzaniSettings, ParasailSettings, RunSettings, WorkflowRun
from workflow.runner import run_parse, run_process
from export import build_source_target_pairs, do_export
from app_state import initialize_app_state
import platform
import psutil
from config import cpu_count


def parse_arguments():
    parser = argparse.ArgumentParser(description="SDT2 CLI")

    parser.add_argument("--input", "-i", required=True, help="Input FASTA file path")

    parser.add_argument("--output", "-o", required=True, help="Output directory path")

    parser.add_argument(
        "--analysis-method",
        choices=["parasail", "lzani"],
        default="parasail",
        help="Analysis method to use (default: parasail)",
    )

    parser.add_argument(
        "--cluster-method",
        default="",
        help="Clustering method (e.g., ward, complete, average)",
    )

    parser.add_argument("--scoring-matrix", help="Scoring matrix for parasail analysis")

    parser.add_argument(
        "--open-penalty", type=int, help="Gap open penalty for parasail analysis"
    )

    parser.add_argument(
        "--extend-penalty", type=int, help="Gap extend penalty for parasail analysis"
    )

    parser.add_argument(
        "--compute-cores",
        type=int,
        default=max(1, min(cpu_count, os.cpu_count() or 1)),
        help="Number of CPU cores to use for parasail analysis",
    )

    parser.add_argument(
        "--lzani-score-type",
        choices=["ani", "tani"],
        default="ani",
        help="Score type for LZANI analysis",
    )

    parser.add_argument("--aw", type=int, help="LZANI aw parameter")
    parser.add_argument("--am", type=int, help="LZANI am parameter")
    parser.add_argument("--mal", type=int, help="LZANI mal parameter")
    parser.add_argument("--msl", type=int, help="LZANI msl parameter")
    parser.add_argument("--mrd", type=float, help="LZANI mrd parameter")
    parser.add_argument("--mqd", type=float, help="LZANI mqd parameter")
    parser.add_argument("--reg", type=int, help="LZANI reg parameter")
    parser.add_argument("--ar", type=float, help="LZANI ar parameter")

    parser.add_argument(
        "--export-alignments", action="store_true", help="Export alignment data"
    )

    parser.add_argument(
        "--prefix",
        default="sdt_output",
        help="Prefix for output files (default: sdt_output)",
    )

    return parser.parse_args()


def get_lzani_exec_path():
    from api.workflow import get_lzani_exec_path

    return get_lzani_exec_path()


def validate_inputs(args):
    if not os.path.isfile(args.input):
        raise FileNotFoundError(f"Input FASTA file not found: {args.input}")

    if not os.path.exists(args.output):
        os.makedirs(args.output, exist_ok=True)
    elif not os.path.isdir(args.output):
        raise NotADirectoryError(f"Output path is not a directory: {args.output}")


def run_cli_workflow(args):
    print("SDT2 CLI - Starting workflow...")

    register_mimetypes()
    validate_inputs(args)

    initialize_app_state(
        debug=os.getenv("DEBUG", "false").lower() == "true",
        platform=dict(
            platform=platform.platform(),
            cores=cpu_count,
            memory=psutil.virtual_memory().total,
        ),
        on_update=lambda _: None,
    )

    with tempfile.TemporaryDirectory() as temp_dir:
        print(f"Processing in temporary directory: {temp_dir}")

        print("Step 1: Parsing FASTA file...")
        parsed_result = run_parse(args.input)

        if parsed_result.error:
            raise Exception(f"FASTA parsing failed: {parsed_result.error}")

        print(f"Successfully parsed {len(parsed_result.seq_dict)} sequences")

        doc_paths = build_document_paths(temp_dir)

        lzani_settings = LzaniSettings(
            exec_path=get_lzani_exec_path() if args.analysis_method == "lzani" else "",
            score_type=args.lzani_score_type,
            aw=args.aw,
            am=args.am,
            mal=args.mal,
            msl=args.msl,
            mrd=args.mrd,
            mqd=args.mqd,
            reg=args.reg,
            ar=args.ar,
        )

        parasail_settings = ParasailSettings(
            process_count=args.compute_cores,
            scoring_matrix=args.scoring_matrix,
            open_penalty=args.open_penalty,
            extend_penalty=args.extend_penalty,
        )

        run_settings = RunSettings(
            fasta_path=args.input,
            doc_paths=doc_paths,
            output_path=temp_dir,
            cluster_method=args.cluster_method,
            analysis_method=args.analysis_method,
            lzani=lzani_settings,
            parasail=parasail_settings,
            export_alignments=args.export_alignments,
            alignment_export_path="",
        )

        workflow_run = WorkflowRun(
            result=parsed_result,
            settings=run_settings,
            progress=0,
        )

        print(f"Step 2: Running {args.analysis_method} analysis...")

        class SimpleCanceled:
            value = False

        canceled = SimpleCanceled()

        final_result = run_process(workflow_run, canceled)  # type: ignore

        if final_result.error:
            raise Exception(f"Workflow processing failed: {final_result.error}")

        print("Step 3: Exporting results...")

        source_target_pairs = build_source_target_pairs(
            temp_dir,
            args.output,
            args.prefix,
            "svg",
            matrix_only=False,
        )

        existing_pairs = [
            (source, target)
            for source, target in source_target_pairs
            if os.path.exists(source)
        ]

        do_export(existing_pairs)

        print(f"Workflow completed successfully, results saved to: {args.output}")

        print("\nExported files:")
        for _, target in existing_pairs:
            print(f"  - {os.path.basename(target)}")


def main():
    try:
        args = parse_arguments()
        run_cli_workflow(args)
    except KeyboardInterrupt:
        print("\nWorkflow interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
