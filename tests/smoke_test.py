from __future__ import annotations

from pathlib import Path

from run import main


bids_dir = Path(__file__).parent.joinpath("data", "ds001")
main()
