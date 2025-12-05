import argparse
import json
import sys
import traceback

from tcr_toolbox.utils.logger import stream_log_function_to_file


def main():
    parser = argparse.ArgumentParser(description="TCR Toolbox CLI")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Subcommand: tcr-assembly
    parser_assembly = subparsers.add_parser("run-tcr-assembly", help="Perform a TCR assembly run")
    parser_assembly.add_argument("config", help="Path to assembly config (json) file")

    # Subcommand: count-umi-cell
    parser_umi = subparsers.add_parser("count-umi-cell", help="Count UMIs per cell barcode")
    parser_umi.add_argument("config", help="Path to umi count config (json) file")

    # Subcommand: count-reads-bulk
    parser_reads = subparsers.add_parser("count-reads-bulk", help="Count reads from bulk gDNA")
    parser_reads.add_argument("config", help="Path to bulk read count config (json) file")

    # Subcommand: reconstruct-tcrs-simple
    parser_reconstruct = subparsers.add_parser("reconstruct-tcrs-simple", help="Reconstruct TCRs from a dataset file")
    parser_reconstruct.add_argument("config", help="Path to reconstruct TCRs config (json) file")

    args = parser.parse_args()

    with open(args.config, "r") as f:
        config = json.load(f)

    if args.command == "run-tcr-assembly":
        from tcr_toolbox.tcr_assembly.tcr_assembly_pipeline import run_tcr_assembly, init_run_dir
        numbered_run_name, run_path_local, log_file = init_run_dir(run_mode=config["run_mode"], run_name=config.get("run_name"), run_path=config.get("run_path"))
        config_wo_run_path_and_name = {k: v for k, v in config.items() if k not in ("run_path", "run_name")}
        try:
            stream_log_function_to_file(run_tcr_assembly, log_file, run_path=run_path_local, numbered_run_name=numbered_run_name, **config_wo_run_path_and_name)
        except Exception:
            traceback.print_exc(file=sys.stderr)
            sys.exit(1)

    elif args.command == "count-umi-cell":
        from tcr_toolbox.sequencing_analysis.align_count import count_umi_cell
        count_umi_cell(**config)

    elif args.command == "count-reads-bulk":
        from tcr_toolbox.sequencing_analysis.align_count import count_reads_bulk
        count_reads_bulk(**config)

    elif args.command == "reconstruct-tcrs-simple":
        from tcr_toolbox.tcr_reconstruction.reconstruct_tcrs_simple import reconstruct_tcrs_simple
        reconstruct_tcrs_simple(**config)
