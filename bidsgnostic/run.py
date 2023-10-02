#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path
from typing import Any

from snakebids.app import SnakeBidsApp

from bidsgnostic import group


def get_parser():
    """Exposes parser for sphinx doc generation, cwd is the docs dir"""
    app = SnakeBidsApp("../")
    return app.parser


def main():
    app = SnakeBidsApp(Path(__file__).resolve().parent.parent)  # to get repository root
    app.run_snakemake()


def group_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="""
        Group level diagnostic tool for BIDS datasets.
        Plots the number of files per particiapnt / session per:
        -  datatype,
        -  datatype / task,
        -  datatype / task and split by any other BIDS entity.
        """,
    )

    parser.add_argument(
        "bids_dir",
        help="""
        Fullpath to the directory with the input dataset
        formatted according to the BIDS standard.
        """,
        nargs=1,
    )
    parser.add_argument(
        "output_dir",
        help="""
        Fullpath to the directory where the output files will be stored.
        If you are running group level analysis this folder should be prepopulated
        with the results of the participant level analysis.
        """,
        nargs=1,
    )
    parser.add_argument(
        "analysis_level",
        help="""
        Level of the analysis that will be performed.
        Multiple participant level analyses can be run independently
        (in parallel) using the same ``output_dir``.
        """,
        choices=["group"],
        type=str,
        nargs=1,
    )
    parser.add_argument(
        "--participant_label",
        help="""
        The label(s) of the participant(s) that should be analyzed.
        The label corresponds to sub-<participant_label> from the BIDS spec
        (so it does not include "sub-").
        If this parameter is not provided all subjects should be analyzed.
        Multiple participants can be specified with a space separated list.
        """,
        type=str,
        nargs="+",
    )
    parser.add_argument(
        "--plot_by",
        help="""
        BIDS entity to split plots by.
        """,
        type=str,
        nargs="+",
        default="",
    )
    parser.add_argument(
        "--log_level",
        help="""
        The log_level level of the reporting that will be printed to the console.
        The default is "0", which means that only warnings and errors will be reported.
        If set to "1", all warnings, errors, and informational messages will be reported.
        If set to "2", all warnings, errors, informational and debug messages will be reported.
        """,
        choices=["0", "1", "2"],
        type=str,
        nargs=1,
        default="0",
    )

    return parser


def main_group(argv: Any = sys.argv) -> None:
    parser = group_parser()

    args = parser.parse_args(argv[1:])

    bids_dir = Path(args.bids_dir[0]).resolve()
    output_dir = Path(args.output_dir[0]).resolve()

    plot_by = args.plot_by
    if isinstance(plot_by, str):
        plot_by = [plot_by]
    if "" in plot_by:
        plot_by.pop(plot_by.index(""))
    if len(plot_by) == 0:
        plot_by = None

    filters = None
    if args.participant_label is not None:
        filters = {"subject": args.participant_label}

    group.main(
        bids_dir,
        output_dir,
        filters=filters,
        plot_by=plot_by,
        log_level=args.log_level[0],
    )


if __name__ == "__main__":
    main()
