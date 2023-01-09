# pyright: reportUndefinedVariable=false
from pathlib import Path

from bidsgnostic.visualize import EventPlotter

this = EventPlotter(
    snakemake.input.events,
    include=snakemake.params.include_events,
    event_column=snakemake.params.event_column,
)

with open(snakemake.output.file, "w") as f:
    f.write(
        f"No events found in {snakemake.input.events} for event {snakemake.params.include_events}"
    )

if this.nb_trial_types > 0:
    this.plot()
    if snakemake.params.log_level == "2":
        this.fig.show()
    Path(snakemake.output.file).parent.mkdir(parents=True, exist_ok=True)
    this.fig.write_html(snakemake.output.file)
