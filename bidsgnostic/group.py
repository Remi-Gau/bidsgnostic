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
    log_level: str = "1",
) -> None:
    if log_level != "0":
        logger.info(f"indexing dataset {bids_dir}")
    layout = BIDSLayout(bids_dir)

    if log_level != "0":
        logger.info(f"plotting events filtered by {filters} split by {plot_by}")

    this = LayoutPlotter(layout, filters=filters)

    suffixes = ["datatype", "task"]
    if plot_by is not None:
        suffixes += plot_by
    for i in suffixes:
        if i in ("", None):
            continue
        with open(output_dir.joinpath(this._fig_name(suffix=i)), "w") as f:
            f.write(f"No data found in {bids_dir} for filters {filters}")

    with open(output_dir.joinpath(this._fig_name(suffix="datatype")), "w") as f:
        f.write(f"No data found in {bids_dir} for filters {filters}")

    if this.df_layout.size > 0:
        this.plot(plot_by=plot_by, show=log_level == "2", output_dir=output_dir)
