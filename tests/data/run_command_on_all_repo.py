"""Script to run a command on all datasets

Assumes that you have cloned all the repositories in the same directory.

Then run this script from within the maintenance-tools dataset.
The script will apply the command to all the repos in `START_DIR`.
"""
import os
from pathlib import Path
from subprocess import run

from rich import print

COMMANDS = [
    "bidsgnostic . ./derivatives/ participant --cores all",
    "bidsgnostic_layout  . ./derivatives/ group --plot_by suffix",
]

DRY_RUN = False
VERBOSE = True

OUTPUT_FILE = Path(__file__).parent.joinpath("outputs", "log.md")
OUTPUT_FILE.unlink(missing_ok=True)

START_DIR = Path(__file__).parent.joinpath("bids-examples")

FOLDERS_TO_SKIP = [
    ".git",
    ".github",
    "maintenance-tools",
    ".vscode",
    "ds000001-fmriprep",
]

DATASET_TO_SKIP = ["7t_trt", "ds000246", "ds000247"]

OUTPUT_FILE.parent.mkdir(exist_ok=True, parents=True)


class DummyResult:
    def __init__(self, stdout="dummy output", stderr=None):
        self.stdout = stdout
        self.stderr = stderr


def run_cmd(cmd, dry_run=True, verbose=True, split=True):

    if split:
        cmd = cmd.split()

    if verbose:
        print(f"[green]{' '.join(cmd)}[/green]")

    if dry_run:
        return DummyResult()

    result = run(cmd, capture_output=True, text=True)

    if verbose and result.stderr:
        print(f"stdout: {result.stdout}")
        print(f"stderr: {result.stderr}")

    return result


def print_to_output(output_file, text):
    if output_file is not None:
        if not output_file.exists():
            with open(output_file, "w", encoding="utf8") as f:
                f.write(f"{text}\n")
        else:
            with open(output_file, "a", encoding="utf8") as f:
                f.write(f"{text}\n")
        return
    else:
        f.write(f"{text}\n")


def main():

    if OUTPUT_FILE is not None:
        with open(OUTPUT_FILE, "w") as log:
            print(f"# Output from '{COMMANDS}'\n", file=log)

    print(f"Applying to folders in:\n{START_DIR}")

    datasets_list = sorted([x for x in START_DIR.iterdir() if x.is_dir()])

    print(f"Running on the following datasets: {datasets_list}")

    for dataset in datasets_list:

        if dataset.is_file():
            continue

        if dataset.name in FOLDERS_TO_SKIP:
            continue

        os.chdir(dataset)

        print(f"\n[blue]{dataset.name}[/blue]")
        print_to_output(output_file=OUTPUT_FILE, text=f"## {dataset.name})")

        for i, cmd in enumerate(COMMANDS):

            if i == 0 and (
                dataset.name.startswith("pet")
                or dataset.name.startswith("asl")
                or dataset.name.startswith("qmri")
                or dataset.name.startswith("micr")
                or dataset.name.endswith("fmriprep")
                or dataset.name in DATASET_TO_SKIP
            ):
                continue

            result = run_cmd(cmd, verbose=VERBOSE, dry_run=DRY_RUN)
            print_to_output(output_file=OUTPUT_FILE, text="\n```")
            print_to_output(output_file=OUTPUT_FILE, text="OUTPUT:")
            print_to_output(output_file=OUTPUT_FILE, text=result.stdout)
            print_to_output(output_file=OUTPUT_FILE, text="\nERROR:")
            print_to_output(output_file=OUTPUT_FILE, text=result.stderr)
            print_to_output(output_file=OUTPUT_FILE, text="```\n")

            print(result.stderr)

    os.chdir(START_DIR)


if __name__ == "__main__":
    main()
