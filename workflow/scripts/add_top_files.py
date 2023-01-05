# pyright: reportUndefinedVariable=false
from bidsgnostic.utils import add_license
from bidsgnostic.utils import add_readme
from bidsgnostic.utils import create_dataset_description

create_dataset_description(snakemake)
add_license(snakemake)
add_readme(snakemake)
