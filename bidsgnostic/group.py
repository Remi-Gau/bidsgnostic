from __future__ import annotations

from pathlib import Path

from bids import BIDSLayout
from loguru import logger

from bidsgnostic.visualize import LayoutPlotter


def main(
    bids_dir: Path,
    output_dir: Path,
    filters: dict[str, list[str]] | None = None,
    plot_by: list[str] | None = None,
    log_level: int = "1",
) -> None:

    if int(log_level) > 0:
        logger.info(f"indexing dataset {bids_dir}")
    layout = BIDSLayout(bids_dir)

    if int(log_level) > 0:
        logger.info(f"plotting events filtered by {filters} split by {plot_by}")
    LayoutPlotter(layout, filters=filters).plot(
        plot_by=plot_by, show=log_level == "2", output_dir=output_dir
    )
