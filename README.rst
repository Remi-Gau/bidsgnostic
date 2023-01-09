

.. image:: https://img.shields.io/badge/license-GPL3-blue.svg
    :alt: license GPL3
.. image:: https://results.pre-commit.ci/badge/github/Remi-Gau/bidsgnostic/main.svg
    :target: https://results.pre-commit.ci/latest/github/Remi-Gau/bidsgnostic/main
    :alt: precommit status
.. image:: https://readthedocs.org/projects/bidsgnostic/badge/?version=latest
    :target: https://bidsgnostic.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation status
.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/psf/black
    :alt: black
.. image:: https://img.shields.io/badge/Sourcery-enabled-brightgreen.svg
    :target: https://sourcery.ai
    :alt: sourcery
.. image:: https://dl.circleci.com/status-badge/img/gh/Remi-Gau/bidsgnostic/tree/main.svg?style=svg
    :target: https://dl.circleci.com/status-badge/redirect/gh/Remi-Gau/bidsgnostic/tree/main
    :alt: circle ci

bidsgnostic
===========

BIDS app to view:

    - the dataset content
    - events file content


Install from github with pip
----------------------------

.. code-block:: bash

    pip install -e git+https://github.com/Remi-Gau/bidsgnostic#egg=bidsgnostic


Note: you can re-run this command to re-install with the latest version

Usage
-----

short form
**********

participant level
_________________

.. code-block:: bash

    bidsgnostic /path/to/bids/dir /path/to/output/dir participant --cores all

.. image:: https://raw.githubusercontent.com/Remi-Gau/bidsgnostic/main/docs/images/sub-01_task-balloonanalogrisktask_run-01_events.png
  :width: 1000
  :alt: example particicant output

See `here for participant level interactive figure as html <https://github.com/Remi-Gau/bidsgnostic/raw/main/docs/images/sub-01_task-balloonanalogrisktask_run-01_events.html>`_.

Figures will not be generated if your events have more than 14 type of events.
If that is the case use the ``include`` argument to select which events you want to see.

group level
___________

.. code-block:: bash

    bidsgnostic_layout /path/to/bids/dir /path/to/output/dir group

.. image:: https://raw.githubusercontent.com/Remi-Gau/bidsgnostic/main/docs/images/dataset-balloonanalogrisktakingtask_splitby-suffix_summary.png
  :width: 1000
  :alt: example group output

See `here for group level interactive figure as html <https://github.com/Remi-Gau/bidsgnostic/raw/main/docs/images/dataset-balloonanalogrisktakingtask_splitby-suffix_summary.html>`_.


all the gory details of the API
*******************************

subject Level
_____________

.. code-block:: text

    usage: bidsgnostic [-h]
                    [--pybidsdb-dir PYBIDSDB_DIR]
                    [--reset-db]
                    [--force-output]
                    [--help-snakemake]
                    [--participant_label PARTICIPANT_LABEL [PARTICIPANT_LABEL ...]]
                    [--exclude_participant_label EXCLUDE_PARTICIPANT_LABEL [EXCLUDE_PARTICIPANT_LABEL ...]]
                    [--derivatives DERIVATIVES [DERIVATIVES ...]]
                    [--log_level LOG_LEVEL]
                    [--filter-events FILTER_EVENTS [FILTER_EVENTS ...]]
                    [--wildcards-events WILDCARDS_EVENTS [WILDCARDS_EVENTS ...]]
                    [--path-events PATH_EVENTS]
                    bids_dir output_dir {participant}

    Snakebids helps build BIDS Apps with Snakemake

    optional arguments:
    -h, --help              show this help message and exit

    STANDARD:
    Standard options for all snakebids apps

    --pybidsdb-dir PYBIDSDB_DIR,
                            Optional path to directory of SQLite databasefile for PyBIDS.
                            If directory is passed and folder exists, indexing is skipped.
                            If reset_db is called, indexing will persist
    --reset-db, --reset_db
                            Reindex existing PyBIDS SQLite database
    --force-output, --force_output
                            Force output in a new directory that already has contents
    --help-snakemake, --help_snakemake
                            Options to Snakemake can also be passed directly at the command-line,
                            use this to print Snakemake usage

    SNAKEBIDS:
    Options for snakebids app

    bids_dir                The directory with the input dataset formatted according to the BIDS standard.

    output_dir              The directory where the output files should be stored.
                            If you are running group level analysis this folder should be prepopulated
                            with the results of the participant level analysis.

    {participant}           Level of the analysis that will be performed.

    --participant_label PARTICIPANT_LABEL [PARTICIPANT_LABEL ...],
                            The label(s) of the participant(s) that should be analyzed.
                            The label corresponds to sub-<participant_label>
                            from the BIDS spec (so it does not include "sub-").
                            If this parameter is not provided all subjects should be analyzed.
                            Multiple participants can be specified with a space separated list.

    --exclude_participant_label EXCLUDE_PARTICIPANT_LABEL [EXCLUDE_PARTICIPANT_LABEL ...],
                            The label(s) of the participant(s) that should be excluded.
                            The label corresponds to sub-<participant_label> from the BIDS spec
                            (so it does not include "sub-").
                            If this parameter is not provided all subjects should be analyzed.
                            Multiple participants can be specified with a space separated list.
    --derivatives DERIVATIVES [DERIVATIVES ...]
                            Path(s) to a derivatives dataset, for folder(s)
                            that contains multiple derivatives datasets (default: False)
    --log_level LOG_LEVEL, --log-level LOG_LEVEL
                            The log_level level of the reporting
                            that will be printed to the console.
                            The default is "0", which means that
                            only warnings and errors will be reported.
                            If set to "1", all warnings, errors,
                            and informational messages will be reported.
                            If set to "2", all warnings, errors,
                            informational and debug messages will be reported.
    --include-events INCLUDE_EVENTS [INCLUDE_EVENTS ...],
                        A list of events name to include in the plot.
                        If none is passed then all events will be included.
                        Must be less than 15 events.

    BIDS FILTERS:
    Filters to customize PyBIDS get() as key=value pairs

    --filter-events FILTER_EVENTS [FILTER_EVENTS ...],
                            (default: suffix=events extension=.tsv)

    INPUT WILDCARDS:
    File path entities to use as wildcards in snakemake

    --wildcards-events WILDCARDS_EVENTS [WILDCARDS_EVENTS ...],
                            (default: subject session acquisition task run)

    PATH OVERRIDE:
    Options for overriding BIDS by specifying absolute paths that include wildcards,
    e.g.: /path/to/my_data/{subject}/t1.nii.gz

    --path-events PATH_EVENTS, --path_events PATH_EVENTS


group Level
___________

.. code-block:: text

    usage: bidsgnostic_layout [-h]
                              [--participant_label PARTICIPANT_LABEL [PARTICIPANT_LABEL ...]]
                              [--plot_by PLOT_BY [PLOT_BY ...]]
                              [--log_level {0,1,2}]
                            bids_dir output_dir {group}

    Group level diagnostic tool for BIDS datasets.
    Plots the number of files per participant / session per:
        - datatype,
        - datatype / task,
        - datatype / task and split by any other BIDS entity.

    positional arguments:

    bids_dir
    output_dir
    {group}

    options:

    -h, --help
    --participant_label PARTICIPANT_LABEL [PARTICIPANT_LABEL ...]
    --log_level {0,1,2}
    --plot_by PLOT_BY [PLOT_BY ...]
                            BIDS entity to split plots by.
