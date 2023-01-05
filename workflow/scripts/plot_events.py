# pyright: reportUndefinedVariable=false
from pathlib import Path

from bids.visualize import EventPlotter

this = EventPlotter(snakemake.input.events)
this.plot()

if snakemake.params.log_level == "2":
    this.fig.show()

Path(snakemake.output.file).parent.mkdir(parents=True, exist_ok=True)
this.fig.write_html(snakemake.output.file)
