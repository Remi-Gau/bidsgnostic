from __future__ import annotations

import json
from pathlib import Path

from loguru import logger

from bidsgnostic import _version

__version__ = _version.get_versions()["version"]


def set_logs(log_dir: Path) -> None:
    error_log_file = log_dir.joinpath("error.log")
    logger.add(error_log_file, filter="ERROR")

    error_log_file = log_dir.joinpath("warning.log")
    logger.add(error_log_file, filter="WARNING")


def create_dataset_description(snakemake) -> None:

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

    output_file = Path(snakemake.input.bids_dir).joinpath(snakemake.output.file)

    output_file.parent.mkdir(parents=True, exist_ok=True)

    with open(output_file, "w") as f:
        json.dump(data, f, indent=4)
