from __future__ import annotations

import json
import shutil
from pathlib import Path

from bids import BIDSLayout
from loguru import logger

from bidsgnostic import _version
from bidsgnostic.visualize import LayoutPlotter

__version__ = _version.get_versions()["version"]


def set_logs(log_dir: Path) -> None:
    error_log_file = log_dir.joinpath("error.log")
    logger.add(error_log_file, filter="ERROR")

    error_log_file = log_dir.joinpath("warning.log")
    logger.add(error_log_file, filter="WARNING")


def add_readme(snakemake):
    output_file = Path(snakemake.output.readme)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    print("TODO: add readme", file=open(output_file, "w"))


def add_license(snakemake) -> None:
    """Copy CCO license to output directory."""
    input_file = str(Path(__file__).parent.joinpath("templates", "CCO"))
    output_file = Path(snakemake.output.license)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(input_file, output_file)


def create_ds_desc(snakemake) -> None:

    data = {
        "Name": "dataset name",
        "BIDSVersion": "1.7.0",
        "DatasetType": "derivatives",
        "License": "CCO",
        "Authors": ["", ""],
        "Acknowledgements": "Special thanks to ",
        "HowToAcknowledge": "",
        "Funding": ["", ""],
        "ReferencesAndLinks": [""],
        "DatasetDOI": "doi:",
        "GeneratedBy": [
            {
                "Name": "bidsgnostic",
                "Version": __version__,
                "Container": {"Type": "", "Tag": ""},
                "Description": "Generate visualizations of BIDS data content for easy review.",
                "CodeURL": "",
            }
        ],
        "SourceDatasets": [
            {
                "DOI": "doi:",
                "URL": "",
                "Version": "",
            }
        ],
    }

    output_file = Path(snakemake.output.desc)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, "w") as f:
        json.dump(data, f, indent=4)


def plot_layout(snakemake) -> None:

    logger.info(f"Indexing BIDS data: {snakemake.input.bids_dir}" "")
    layout = BIDSLayout(snakemake.input.bids_dir)
    this = LayoutPlotter(layout)
    output_file = Path(snakemake.output.datatype_html)
    this.plot(
        output_dir=output_file.parent,
        show=snakemake.params.log_level == "2",
    )
    datatype_html = output_file.parent.glob("*_splitby-datatype_summary.html")
    shutil.copy(next(datatype_html), snakemake.output.datatype_html)

    task_html = output_file.parent.glob("*_splitby-task_summary.html")
    shutil.copy(next(task_html), snakemake.output.task_html)
